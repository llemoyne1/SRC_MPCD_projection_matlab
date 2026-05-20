function [stateOut, diag] = mpcd_step_projection_piston(state, params, stepIndex)
%MPCD_STEP_PROJECTION_PISTON Moving-piston classic step plus Q6/Q9 projection.
%
%   [stateOut, diag] = mpcd_step_projection_piston(state, params, stepIndex)
%
% This wrapper keeps the projection core used by Taylor-Green and Poiseuille:
%   classic SRC/MPCD step -> optional Q6 velocity projection -> optional Q9
%   low-k mass-flux correction -> optional final cleanup -> thermostat.
%
% Supported method names:
%   'classic'        : moving piston + wallVP-v2, no projection.
%   'q6_projection'  : Q6 velocity projection only.
%   'q9_full' / 'q9' : Q6 + Q9 low-k mass-flux correction, i.e. the full
%                     Q9 model used for TG/Poiseuille validations.
%   'q9_historical'  : same as q9_full but defaults to the legacy
%                     periodic-x/Neumann-y mass-flux operator and lowpass_fft.
%   'q9_mass_only'   : piston-specific diagnostic branch kept for comparison:
%                     Q9 mass-flux correction without preliminary Q6.
%
% The active projection domain follows the moving piston height Ly=yTop.
% The fixed Nx-by-Ny grid means that a homogeneous cell occupancy gamma
% corresponds to an increasing physical density as yTop decreases.

if nargin < 3 || isempty(stepIndex)
    stepIndex = get_param(params, 'pistonStepIndex', 1);
end

method = normalize_method(get_param(params, 'method', 'classic'));

[stepParams, pistonInfo] = piston_runtime_params(params, state, stepIndex);

% For projection methods, follow the TG/Poiseuille ordering: the classic SRC
% substep is not thermostatted; the thermostat is applied once after the
% Q6/Q9 projection block.  For the pure classic control case, keep the
% user/default thermostatAfterStep behavior.
if ~strcmp(method, 'classic')
    stepParams.thermostatAfterStep = false;
end

[stateClassic, classicDiag] = mpcd_step_classic_piston(state, stepParams);
if isfield(stateClassic, 'piston') && isstruct(stateClassic.piston)
    pistonInfo = merge_struct_fields(pistonInfo, stateClassic.piston);
end

if strcmp(method, 'classic')
    stateOut = stateClassic;
    stateOut.piston = pistonInfo;
    diag = struct();
    diag.method = 'classic';
    diag.classic = classicDiag;
    diag.wallInfo = classicDiag.wallInfo;
    diag.piston = pistonInfo;
    diag.projectionApplied = false;
    diag.meanVyAfterProjection = mean(stateOut.v(:,2), 'omitnan');
    diag.kBTCellAfterProjection = classicDiag.kBTCell;
    diag.populationAfterClassic = classicDiag.population;
    diag.populationAfterProjection = classicDiag.population;
    diag.rmsMassFluxDivResidual = NaN;
    diag.rmsMassFluxDivBefore = NaN;
    diag.rmsMassFluxDivProjectedAfter = NaN;
    return;
end

projParams = configure_projection_params(stepParams, params, pistonInfo, stateClassic, method);

% Optional virial EOS closure is applied after the validated Q6/Q9 block.
% If the kick is active, postpone the thermostat until after the kick so the
% final state remains at the target temperature without hiding the virial
% impulse in the projection subroutine.
virialDiagnosticsEnable = logical(get_param(params, 'virialDiagnosticsEnable', false));
virialKickEnable = logical(get_param(params, 'virialKickEnable', false));
thermostatAfterVirial = virialKickEnable && logical(get_param(params, 'thermostatAfterProjection', true));
if thermostatAfterVirial
    projParams.thermostatAfterProjection = false;
end

[stateOut, projectionDiag] = mpcd_apply_q9_projection_channel(stateClassic, projParams);
stateOut.piston = pistonInfo;

virialInfo = struct();
if virialDiagnosticsEnable || virialKickEnable
    virialParams = projParams;
    virialParams.Ly0 = pistonInfo.Ly0;
    virialParams.Ly = pistonInfo.yTop;
    virialParams.useMovingPiston = true;
    virialParams.boundary_top = 'piston';
    virialParams.virialDiagnosticsEnable = virialDiagnosticsEnable || virialKickEnable;
    virialParams.virialKickEnable = virialKickEnable;
    [stateOut, virialInfo] = mpcd_apply_virial_pressure_kick_channel(stateOut, virialParams);
    stateOut.piston = pistonInfo;

    if thermostatAfterVirial
        [stateOut.v, virialThermostatInfo] = projection_apply_cell_thermostat(stateOut.x, stateOut.v, virialParams, ...
            'periodicX', true, 'periodicY', false);
        virialInfo.thermostatAfterVirial = true;
        virialInfo.thermostatInfo = virialThermostatInfo;
        virialInfo.kBTAfterVirialThermostat = getf(virialThermostatInfo, 'meanKBTAfter', NaN);
    else
        virialInfo.thermostatAfterVirial = false;
    end
end

diag = projectionDiag;
diag.method = method;
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.piston = pistonInfo;
diag.projectionApplied = true;
diag.virial = virialInfo;
diag.virialDiagnosticsEnable = virialDiagnosticsEnable;
diag.virialKickEnable = virialKickEnable;
diag.thermostatAfterVirial = thermostatAfterVirial;
if isstruct(virialInfo) && isfield(virialInfo, 'enabled')
    diag.virialKvirial = getf(virialInfo, 'Kvirial', NaN);
    diag.virialBeta = getf(virialInfo, 'betaVirial', NaN);
    diag.virialPvirMean = getf(virialInfo, 'PvirMean', NaN);
    diag.virialPtotMean = getf(virialInfo, 'PtotMean', NaN);
    diag.virialDuRms = getf(virialInfo, 'duVirialAppliedRms', NaN);
    diag.virialResidualMomentumNorm = getf(virialInfo, 'residualMomentumKickNorm', NaN);
end
diag.activeProjectionParams = struct('Lx', projParams.Lx, 'Ly', projParams.Ly, ...
    'Nx', projParams.Nx, 'Ny', projParams.Ny, 'gamma', projParams.gamma, ...
    'method', method, ...
    'projectionStrength', get_param(projParams, 'projectionStrength', NaN), ...
    'massFluxProjectionMode', char(string(get_param(projParams, 'massFluxProjectionMode', 'off'))), ...
    'massFluxApplyAfterVelocityProjection', logical(get_param(projParams, 'massFluxApplyAfterVelocityProjection', false)), ...
    'massFluxProjectionOperator', char(string(get_param(projParams, 'massFluxProjectionOperator', 'legacy_periodic_x_neumann_y'))), ...
    'massFluxTargetFilter', char(string(get_param(projParams, 'massFluxTargetFilter', 'none'))), ...
    'massFluxFinalVelocityProjectionCleanup', logical(get_param(projParams, 'massFluxFinalVelocityProjectionCleanup', false)), ...
    'massFluxFinalVelocityProjectionStrength', get_param(projParams, 'massFluxFinalVelocityProjectionStrength', NaN), ...
    'pistonQ9TargetMode', char(string(get_param(params, 'pistonQ9TargetMode', 'reference_gamma'))), ...
    'freezeEllipticMetric', logical(get_param(params, 'pistonQ9FreezeEllipticMetric', false)), ...
    'massFluxEllipticLxOverride', getf(projParams, 'massFluxEllipticLxOverride', NaN), ...
    'massFluxEllipticLyOverride', getf(projParams, 'massFluxEllipticLyOverride', NaN));
end

function projParams = configure_projection_params(stepParams, params, pistonInfo, stateClassic, method)
projParams = stepParams;
projParams.Ly0 = pistonInfo.Ly0;
projParams.Ly = pistonInfo.yTop;
projParams.useMovingPiston = true;
projParams.boundary_top = 'piston';

% The general elliptic operator can optionally keep a fixed metric for
% benchmarking speed.  Default false here: q9_full should reproduce the
% mathematical model on the active moving-height domain.  Set true only for
% a performance sensitivity test.
if logical(get_param(params, 'pistonQ9FreezeEllipticMetric', false))
    projParams.massFluxEllipticLxOverride = get_param(params, 'Lx', projParams.Lx);
    projParams.massFluxEllipticLyOverride = pistonInfo.Ly0;
else
    if isfield(projParams, 'massFluxEllipticLxOverride'), projParams = rmfield(projParams, 'massFluxEllipticLxOverride'); end
    if isfield(projParams, 'massFluxEllipticLyOverride'), projParams = rmfield(projParams, 'massFluxEllipticLyOverride'); end
end

projParams.projectionEnable = true;
projParams.thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', true));
projParams.projectionInterpolationMethod = get_param(params, 'projectionInterpolationMethod', 'nearest');
projParams.projectionMomentumCorrectionEnable = get_param(params, 'projectionMomentumCorrectionEnable', true);
projParams.projectionMomentumCorrectionMode = get_param(params, 'projectionMomentumCorrectionMode', 'particle_global_exact');
projParams.computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));
projParams.projectionTransportDiagnosticsEnable = logical(get_param(params, 'projectionTransportDiagnosticsEnable', true));
projParams.densityTransportDiagnosticsEnable = logical(get_param(params, 'densityTransportDiagnosticsEnable', true));
projParams.massFluxLowKMaxIndex = get_param(params, 'massFluxLowKMaxIndex', get_param(params, 'lowKMaxIndex', 2));
projParams.lowKMaxIndex = get_param(params, 'lowKMaxIndex', projParams.massFluxLowKMaxIndex);
projParams.massFluxMinCellCount = get_param(params, 'massFluxMinCellCount', 1.0);
projParams.massFluxProjectionRegularization = get_param(params, 'massFluxProjectionRegularization', 1e-12);
projParams.massFluxBCLeft = 'periodic';
projParams.massFluxBCRight = 'periodic';
projParams.massFluxBCBottom = 'wall';
projParams.massFluxBCTop = 'wall';
projParams.massFluxEllipticUseSolveCache = logical(get_param(params, 'massFluxEllipticUseSolveCache', true));
projParams.massFluxEllipticSolver = get_param(params, 'massFluxEllipticSolver', 'cached_backslash');
projParams.massFluxEllipticFactorization = get_param(params, 'massFluxEllipticFactorization', 'auto');

% Homogeneous-density target.  With a fixed active grid, reference gamma is
% the consistent quasi-incompressible target; physical density still changes
% as Ly=yTop changes.  The instant_mean option is kept for robustness if the
% particle count is not exactly gamma*Nx*Ny.
targetMode = lower(strrep(char(string(get_param(params, 'pistonQ9TargetMode', 'reference_gamma'))), '-', '_'));
switch targetMode
    case {'reference_gamma','fixed_gamma','initial_gamma'}
        projParams.gamma = get_param(params, 'gamma', size(stateClassic.x, 1) / max(projParams.Nx * projParams.Ny, 1));
    case {'instant_mean_occupancy','instant_mean','mean_occupancy','preserve_global_compression'}
        projParams.gamma = size(stateClassic.x, 1) / max(projParams.Nx * projParams.Ny, 1);
    otherwise
        error('Unknown pistonQ9TargetMode: %s', targetMode);
end

switch method
    case 'q6_projection'
        projParams.projectionStrength = get_param(params, 'projectionStrength', 1.0);
        projParams.massFluxProjectionMode = 'off';
        projParams.massFluxProjectionStrength = 0.0;
        projParams.massFluxDensityRelaxationBeta = 0.0;
        projParams.massFluxApplyAfterVelocityProjection = false;
        projParams.massFluxFinalVelocityProjectionCleanup = false;
        projParams.massFluxFinalVelocityProjectionStrength = 0.0;
        projParams.massFluxTargetFilter = 'none';

    case 'q9_full'
        % Full model used by TG/Poiseuille: Q6 first, then Q9 mass-flux
        % correction.  The operator/filter remain parametric so the same
        % wrapper can test the legacy and elliptic paths.
        projParams.projectionStrength = get_param(params, 'projectionStrength', 1.0);
        projParams.massFluxProjectionMode = get_param(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk');
        projParams.massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', 1.0);
        projParams.massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 5.0e-4);
        projParams.massFluxApplyAfterVelocityProjection = true;
        projParams.massFluxFinalVelocityProjectionCleanup = logical(get_param(params, 'massFluxFinalVelocityProjectionCleanup', false));
        projParams.massFluxFinalVelocityProjectionStrength = get_param(params, 'massFluxFinalVelocityProjectionStrength', 0.5);
        projParams.massFluxProjectionOperator = get_param(params, 'massFluxProjectionOperator', 'general_bc');
        projParams.massFluxTargetFilter = get_param(params, 'massFluxTargetFilter', 'elliptic_lowpass');

    case 'q9_historical'
        % Historical piston Q9 defaults from the q6/q9 snapshot.
        projParams.projectionStrength = get_param(params, 'projectionStrength', 1.0);
        projParams.massFluxProjectionMode = get_param(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk');
        projParams.massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', 1.0);
        projParams.massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 0.002);
        projParams.massFluxApplyAfterVelocityProjection = true;
        projParams.massFluxFinalVelocityProjectionCleanup = logical(get_param(params, 'massFluxFinalVelocityProjectionCleanup', true));
        projParams.massFluxFinalVelocityProjectionStrength = get_param(params, 'massFluxFinalVelocityProjectionStrength', 0.5);
        projParams.massFluxProjectionOperator = get_param(params, 'massFluxProjectionOperator', 'periodic_x_neumann_y');
        projParams.massFluxTargetFilter = get_param(params, 'massFluxTargetFilter', 'lowpass_fft');

    case 'q9_mass_only'
        % Diagnostic branch only: no preliminary Q6 projection.  This is not
        % the validated full Q9 model, but is useful for sensitivity tests.
        projParams.projectionStrength = 0.0;
        projParams.massFluxProjectionMode = get_param(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk');
        projParams.massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', 1.0);
        projParams.massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 5.0e-4);
        projParams.massFluxApplyAfterVelocityProjection = false;
        projParams.massFluxFinalVelocityProjectionCleanup = false;
        projParams.massFluxFinalVelocityProjectionStrength = 0.0;
        projParams.massFluxProjectionOperator = get_param(params, 'massFluxProjectionOperator', 'general_bc');
        projParams.massFluxTargetFilter = get_param(params, 'massFluxTargetFilter', 'elliptic_lowpass');

    otherwise
        error('Unsupported piston projection method: %s', method);
end
end

function method = normalize_method(value)
method = lower(strrep(char(string(value)), '-', '_'));
switch method
    case {'classic','none','src'}
        method = 'classic';
    case {'q6','q6_projection','projection','projection_only','q6_full'}
        method = 'q6_projection';
    case {'q9','q9_full','q9_projection','q6_q9','full_q9','q9_validated'}
        method = 'q9_full';
    case {'q9_historical','historical_q9','q9_legacy','legacy_q9'}
        method = 'q9_historical';
    case {'q9_mass_only','mass_only','q9_piston_mass_only'}
        method = 'q9_mass_only';
    otherwise
        error('Unknown piston method: %s', method);
end
end

function [p, info] = piston_runtime_params(params, state, stepIndex)
p = params;
Ly0 = get_param(params, 'Ly0', get_param(params, 'LyReference', params.Ly));
y0 = get_param(params, 'pistonY0', Ly0);
yMin = get_param(params, 'pistonYMin', 0.95 * Ly0);
vy = get_param(params, 'pistonVy', 0.0);
dt = params.dt;
stopOnMin = logical(get_param(params, 'pistonStopOnMin', true));
motionMode = lower(char(string(get_param(params, 'pistonMotionMode', 'linear'))));

if strcmp(motionMode, 'cycle') || strcmp(motionMode, 'compression_hold_decompression')
    [yPrev, yTop, up, phaseName, phaseIndex] = piston_cycle_kinematics(params, stepIndex, Ly0, y0, yMin, dt);
elseif strcmp(motionMode, 'staircase') || strcmp(motionMode, 'staircase_eos')
    [yPrev, yTop, up, phaseName, phaseIndex] = piston_staircase_kinematics(params, stepIndex, Ly0, y0, dt);
elseif any(strcmp(motionMode, {'smooth_ramp','smoothstep_ramp','quasistatic_ramp'}))
    [yPrev, yTop, up, phaseName, phaseIndex] = piston_smooth_ramp_kinematics(params, stepIndex, Ly0, y0, dt);
else
    if isfield(state, 'piston') && isfield(state.piston, 'yTop') && ~isempty(state.piston.yTop)
        yPrev = state.piston.yTop;
    else
        yPrev = y0 + max(stepIndex - 1, 0) * dt * vy;
    end

    yTop = y0 + stepIndex * dt * vy;
    if stopOnMin
        if vy < 0
            yTop = max(yMin, yTop);
        elseif vy > 0
            yTop = min(yMin, yTop);
        end
    end

    up = vy;
    if stopOnMin && ((vy < 0 && yTop <= yMin + eps(Ly0)) || (vy > 0 && yTop >= yMin - eps(Ly0)))
        up = 0.0;
    end
    phaseName = 'linear';
    phaseIndex = 1;
end

yTop = min(max(yTop, eps(Ly0)), Ly0);
yPrev = min(max(yPrev, eps(Ly0)), Ly0);

p.Ly0 = Ly0;
p.LyReference = Ly0;
p.pistonYPrev = yPrev;
p.pistonYCurrent = yTop;
p.pistonVyCurrent = up;
p.Ly = yTop;
p.useMovingPiston = true;
p.boundary_top = 'piston';

info = struct();
info.yTop = yTop;
info.yPrev = yPrev;
info.Up = up;
info.Ly0 = Ly0;
info.activeHeight = yTop;
info.compression = 1.0 - yTop / Ly0;
info.stepIndex = stepIndex;
info.t = stepIndex * dt;
info.motionMode = motionMode;
info.phase = phaseName;
info.phaseIndex = phaseIndex;
info.dV = p.Lx * (yTop - yPrev);
info.dVdt = info.dV / max(dt, eps);
info.activeArea = p.Lx * yTop;
info.rhoPhysicalMean = size(state.x, 1) / max(info.activeArea, eps);
info.gammaMeanGeometric = size(state.x, 1) / max(p.Nx * p.Ny, 1);
end


function [yPrev, yTop, up, phaseName, phaseIndex] = piston_smooth_ramp_kinematics(params, stepIndex, Ly0, y0, dt)
initialHold = max(0, round(get_param(params, 'pistonSmoothRampInitialHoldSteps', 0)));
rampSteps = max(1, round(get_param(params, 'pistonSmoothRampSteps', get_param(params, 'nSteps', 1))));
finalHold = max(0, round(get_param(params, 'pistonSmoothRampFinalHoldSteps', 0))); %#ok<NASGU>
cFinal = get_param(params, 'pistonSmoothRampCompressionFinal', get_param(params, 'pistonCompressionTarget', 0.0));
cFinal = min(max(double(cFinal), 0.0), 0.95);
profile = lower(char(string(get_param(params, 'pistonSmoothRampProfile', 'smoothstep'))));

k0 = max(stepIndex - 1, 0);
k1 = max(stepIndex, 0);
[cPrev, ~, ~] = smooth_ramp_compression_at(k0, cFinal, initialHold, rampSteps, profile);
[cTop, phaseName, phaseIndex] = smooth_ramp_compression_at(k1, cFinal, initialHold, rampSteps, profile);
yPrev = y0 * (1.0 - cPrev);
yTop = y0 * (1.0 - cTop);
if stepIndex <= 0
    up = 0.0;
else
    up = (yTop - yPrev) / max(dt, eps);
end
yPrev = min(max(yPrev, eps(Ly0)), Ly0);
yTop = min(max(yTop, eps(Ly0)), Ly0);
end

function [c, phaseName, phaseIndex] = smooth_ramp_compression_at(k, cFinal, initialHold, rampSteps, profile)
if k <= initialHold
    c = 0.0;
    phaseName = 'smooth_ramp_initial_hold';
    phaseIndex = 1;
    return;
end
q = (k - initialHold) / max(rampSteps, 1);
if q < 1.0
    q = min(max(q, 0.0), 1.0);
    switch profile
        case {'smoothstep','s_curve','scurve'}
            s = q*q*(3.0 - 2.0*q);
        case {'smootherstep'}
            s = q*q*q*(q*(q*6.0 - 15.0) + 10.0);
        case {'linear'}
            s = q;
        otherwise
            error('Unknown pistonSmoothRampProfile: %s', profile);
    end
    c = cFinal * s;
    phaseName = 'smooth_ramp';
    phaseIndex = 2;
else
    c = cFinal;
    phaseName = 'smooth_ramp_final_hold';
    phaseIndex = 3;
end
end

function [yPrev, yTop, up, phaseName, phaseIndex] = piston_staircase_kinematics(params, stepIndex, Ly0, y0, dt)
levels = get_param(params, 'pistonStaircaseCompressionList', [0 0.0025 0.005 0.01]);
levels = double(levels(:).');
if isempty(levels)
    levels = 0;
end
levels(~isfinite(levels)) = 0;
levels = min(max(levels, 0), 0.95);
if levels(1) ~= 0
    levels = [0, levels];
end
initialHold = max(0, round(get_param(params, 'pistonStaircaseInitialHoldSteps', 0)));
moveSteps = max(1, round(get_param(params, 'pistonStaircaseMoveSteps', 100)));
holdSteps = max(0, round(get_param(params, 'pistonStaircaseHoldSteps', 400)));
[yPrev, ~, ~] = staircase_position(max(stepIndex - 1, 0), levels, y0, initialHold, moveSteps, holdSteps);
[yTop, phaseName, phaseIndex] = staircase_position(max(stepIndex, 0), levels, y0, initialHold, moveSteps, holdSteps);
if stepIndex <= 0
    up = 0.0;
else
    up = (yTop - yPrev) / max(dt, eps);
end
yPrev = min(max(yPrev, eps(Ly0)), Ly0);
yTop = min(max(yTop, eps(Ly0)), Ly0);
end

function [y, phaseName, phaseIndex] = staircase_position(k, levels, y0, initialHold, moveSteps, holdSteps)
if k <= initialHold
    phaseName = 'staircase_initial_hold';
    phaseIndex = 1;
    y = y0 * (1.0 - levels(1));
    return;
end
rem = k - initialHold;
prev = levels(1);
for j = 2:numel(levels)
    target = levels(j);
    if rem <= moveSteps
        q = rem / max(moveSteps, 1);
        c = prev + (target - prev) * q;
        phaseName = 'staircase_move';
        phaseIndex = j;
        y = y0 * (1.0 - c);
        return;
    end
    rem = rem - moveSteps;
    if rem <= holdSteps
        phaseName = 'staircase_hold';
        phaseIndex = j;
        y = y0 * (1.0 - target);
        return;
    end
    rem = rem - holdSteps;
    prev = target;
end
phaseName = 'staircase_done';
phaseIndex = numel(levels);
y = y0 * (1.0 - levels(end));
end

function [yPrev, yTop, up, phaseName, phaseIndex] = piston_cycle_kinematics(params, stepIndex, Ly0, y0, yMin, dt)
nC = max(1, round(get_param(params, 'pistonCycleCompressSteps', get_param(params, 'pistonCompressSteps', 5000))));
nH = max(0, round(get_param(params, 'pistonCycleHoldSteps', get_param(params, 'pistonHoldSteps', 0))));
nD = max(1, round(get_param(params, 'pistonCycleDecompressSteps', get_param(params, 'pistonDecompressSteps', nC))));
nF = max(0, round(get_param(params, 'pistonCycleFinalHoldSteps', get_param(params, 'pistonFinalHoldSteps', 0))));

k0 = max(stepIndex - 1, 0);
k1 = max(stepIndex, 0);
yPrev = cycle_position(k0, y0, yMin, nC, nH, nD, nF);
yTop = cycle_position(k1, y0, yMin, nC, nH, nD, nF);
up = (yTop - yPrev) / max(dt, eps);
[phaseName, phaseIndex] = cycle_phase(k1, nC, nH, nD, nF);
yPrev = min(max(yPrev, eps(Ly0)), Ly0);
yTop = min(max(yTop, eps(Ly0)), Ly0);
end

function y = cycle_position(k, y0, yMin, nC, nH, nD, nF) %#ok<INUSD>
if k <= nC
    y = y0 + (yMin - y0) * (k / max(nC, 1));
elseif k <= nC + nH
    y = yMin;
elseif k <= nC + nH + nD
    q = (k - nC - nH) / max(nD, 1);
    y = yMin + (y0 - yMin) * q;
else
    y = y0;
end
end

function [name, idx] = cycle_phase(k, nC, nH, nD, nF)
if k <= nC
    name = 'compression'; idx = 1;
elseif k <= nC + nH
    name = 'hold_compressed'; idx = 2;
elseif k <= nC + nH + nD
    name = 'decompression'; idx = 3;
elseif k <= nC + nH + nD + nF
    name = 'hold_relaxed'; idx = 4;
else
    name = 'done'; idx = 5;
end
end


function dst = merge_struct_fields(dst, src)
if ~isstruct(src)
    return;
end
names = fieldnames(src);
for i = 1:numel(names)
    dst.(names{i}) = src.(names{i});
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end

function value = getf(s, name, defaultValue)
if nargin < 3
    defaultValue = NaN;
end
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
