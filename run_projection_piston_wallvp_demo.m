function out = run_projection_piston_wallvp_demo(varargin)
%RUN_PROJECTION_PISTON_WALLVP_DEMO Moving-piston wallVP-v2 validation demo.
%
%   out = run_projection_piston_wallvp_demo('method','classic')
%   out = run_projection_piston_wallvp_demo('method','q6_projection')
%   out = run_projection_piston_wallvp_demo('method','q9_full')
%   out = run_projection_piston_wallvp_demo('method','q9_mass_only')
%
% This runner is a piston application wrapper.  The core Q6/Q9 model remains
% mpcd_apply_q9_projection_channel, as used by TG and Poiseuille.  The
% q9_mass_only branch is kept only as a diagnostic sensitivity mode; the
% validated Q9 method is method='q9_full' or method='q9'.

params = parse_inputs(varargin{:});
params = set_default_params(params);

rng(params.seed);
initParams = params;
initParams.Ly = params.pistonY0;
[state, initInfo] = projection_initialize_particles_poiseuille(initParams, ...
    'initialPopulationMode', params.initialPopulationMode, ...
    'initialVelocityZeroGlobalMean', params.initialVelocityZeroGlobalMean, ...
    'initialMeanVelocityX', params.initialMeanVelocityX, ...
    'initialMeanVelocityY', params.initialMeanVelocityY);
state.piston = struct('yTop', params.pistonY0, 'yPrev', params.pistonY0, ...
    'Up', params.pistonVy, 'Ly0', params.Ly0, 'compression', 0.0, ...
    'activeHeight', params.pistonY0);
Np0 = size(state.x, 1);

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
rows = cell(nSamples, 1);
UxProfiles = nan(params.Ny, nSamples);
UyProfiles = nan(params.Ny, nSamples);
NProfiles = nan(params.Ny, nSamples);
if params.storeDensityMaps
    NMaps = nan(params.Nx, params.Ny, nSamples);
else
    NMaps = [];
end

isamp = 0;
lastStepDiag = [];
actualStep = 0;
wallClockTic = tic;
fig = [];
if params.visualEnable
    fig = figure(params.visualFigureId);
    clf(fig);
    setappdata(fig, 'stopRequested', false);
end

for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        sampleDiag = piston_sample_diagnostics(state, params, it, lastStepDiag);
        rows{isamp} = sampleDiag.row;
        UxProfiles(:, isamp) = sampleDiag.UxProfile;
        UyProfiles(:, isamp) = sampleDiag.UyProfile;
        NProfiles(:, isamp) = sampleDiag.NProfile;
        if params.storeDensityMaps
            NMaps(:, :, isamp) = sampleDiag.G.N;
        end
    end

    if params.visualEnable && (mod(it, params.visualEvery) == 0 || it == params.nSteps)
        diagTableTmp = rows_to_table(rows(1:isamp));
        projection_piston_visualize_frame(state, params, diagTableTmp, it, sampleDiag, fig);
        drawnow limitrate;
        pause(params.visualPause);
        if isgraphics(fig) && isappdata(fig, 'stopRequested') && getappdata(fig, 'stopRequested')
            fprintf('Stop requested from piston visualization at step %d.\n', it);
            break;
        end
    end

    if it == params.nSteps
        actualStep = it;
        break;
    end

    nextStep = it + 1;
    stepParams = params;
    stepParams.computeDiagnostics = params.computeFullDiagnosticsEveryStep || ...
        mod(nextStep, params.sampleEvery) == 0 || ...
        (params.progressEvery > 0 && mod(nextStep, params.progressEvery) == 0) || ...
        nextStep == params.nSteps;

    [state, lastStepDiag] = mpcd_step_projection_piston(state, stepParams, nextStep);
    actualStep = nextStep;

    if params.progressEvery > 0 && mod(nextStep, params.progressEvery) == 0
        sd = piston_sample_diagnostics(state, params, nextStep, lastStepDiag);
        fprintf('  piston %s step %d/%d, t=%.4g, yTop=%.6g, comp=%.4g, rho=%.4g, stdN=%.3g, lowK=%.3e, Pkin=%.4g, Ptop=%.4g, kBT=%.4g, elapsed=%.1fs\n', ...
            upper(char(string(params.method))), nextStep, params.nSteps, nextStep*params.dt, ...
            sd.row.yTop, sd.row.compression, sd.row.rhoPhysicalMean, sd.row.stdN, ...
            sd.row.lowKDensity, sd.row.PkinMean, sd.row.pressureTopWallTotal, sd.row.kBTCell, toc(wallClockTic));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_piston_wallvp_demo:TimeLimit', ...
            'Stopping early at step %d/%d after %.1f s because maxWallClockSeconds=%.1f.', ...
            actualStep, params.nSteps, toc(wallClockTic), params.maxWallClockSeconds);
        break;
    end
end

elapsed = toc(wallClockTic);
validRows = rows(~cellfun(@isempty, rows));
diagTable = rows_to_table(validRows);
validSamples = height(diagTable);
UxProfiles = UxProfiles(:, 1:validSamples);
UyProfiles = UyProfiles(:, 1:validSamples);
NProfiles = NProfiles(:, 1:validSamples);
if params.storeDensityMaps
    NMaps = NMaps(:, :, 1:validSamples);
end

out = struct();
out.params = params;
out.initialization = initInfo;
out.state = state;
out.diagTable = diagTable;
out.UxProfiles = UxProfiles;
out.UyProfiles = UyProfiles;
out.NProfiles = NProfiles;
out.NMaps = NMaps;
out.actualStep = actualStep;
out.elapsedWallClock = elapsed;
out.summary = summarize_piston(out, Np0);
out.compressibility = analyze_piston_compressibility(out, ...
    'fitStartFraction', get_param(params, 'pistonCompressibilityFitStartFraction', 0.2), ...
    'fitEndFraction', get_param(params, 'pistonCompressibilityFitEndFraction', 1.0));

fprintf('\n=== run_projection_piston_wallvp_demo (%s) ===\n', upper(char(string(params.method))));
fprintf('Np/grid/gamma                : %d / %d x %d / %.6g\n', Np0, params.Nx, params.Ny, params.gamma);
fprintf('dt, kBT                     : %.6g / %.6g\n', params.dt, params.kBT);
fprintf('piston y0/yMin/final        : %.12g / %.12g / %.12g\n', params.pistonY0, params.pistonYMin, out.summary.finalYTop);
fprintf('compression final           : %.12g\n', out.summary.finalCompression);
fprintf('wallVP geometry/density     : %s / %.6g\n', char(string(params.wallVirtualParticlesGeometryMode)), params.wallVirtualParticlesDensityFactor);
fprintf('steps completed             : %d / %d in %.2f s\n', out.actualStep, params.nSteps, out.elapsedWallClock);
fprintf('mean/final rho physical     : %.12g / %.12g\n', out.summary.meanRhoPhysicalMean, out.summary.finalRhoPhysicalMean);
fprintf('mean/final rho rel. error   : %.12g / %.12g\n', out.summary.meanRhoPhysicalRelError, out.summary.finalRhoPhysicalRelError);
fprintf('mean/final std(N)           : %.12g / %.12g\n', out.summary.meanStdN, out.summary.finalStdN);
fprintf('mean/final low-k density    : %.12g / %.12g\n', out.summary.meanLowKDensity, out.summary.finalLowKDensity);
fprintf('mean/final Pkin             : %.12g / %.12g\n', out.summary.meanPkinMean, out.summary.finalPkinMean);
fprintf('mean/final Pkin/(rho kBT)   : %.12g / %.12g\n', out.summary.meanPkinIdealRatio, out.summary.finalPkinIdealRatio);
if isfield(out.summary, 'finalPtotMean') && isfinite(out.summary.finalPtotMean)
    fprintf('mean/final Pvir             : %.12g / %.12g\n', out.summary.meanPvirMean, out.summary.finalPvirMean);
    fprintf('mean/final Ptot             : %.12g / %.12g\n', out.summary.meanPtotMean, out.summary.finalPtotMean);
    fprintf('mean/final virial du rms    : %.12g / %.12g\n', out.summary.meanVirialDuRms, out.summary.finalVirialDuRms);
    fprintf('mean/final virial limiter   : %.12g / %.12g  (duMax %.4g)\n', ...
        out.summary.meanVirialLimitedCellFraction, out.summary.finalVirialLimitedCellFraction, ...
        out.summary.finalVirialLimiterDuMax);
end
fprintf('mean/final kBT cell         : %.12g / %.12g\n', out.summary.meanKBTCell, out.summary.finalKBTCell);
fprintf('mean/final top wall pressure: %.12g / %.12g\n', out.summary.meanTopWallPressureTotal, out.summary.finalTopWallPressureTotal);
fprintf('  impact / wallVP top mean  : %.12g / %.12g\n', out.summary.meanTopWallPressureImpact, out.summary.meanTopWallPressureVP);
fprintf('mean/final piston work      : %.12g / %.12g\n', out.summary.meanPistonWorkOnFluidCumulative, out.summary.finalPistonWorkOnFluidCumulative);
fprintf('mean/final Pkin top layer   : %.12g / %.12g\n', out.summary.meanPkinTopLayerMean, out.summary.finalPkinTopLayerMean);
if isstruct(out.compressibility) && isfield(out.compressibility, 'KeffPkin')
    fprintf('Keff Pkin / topWall total   : %.12g / %.12g  (fit R2 %.4g / %.4g)\n', ...
        out.compressibility.KeffPkin, out.compressibility.KeffTopWallTotal, ...
        out.compressibility.R2Pkin, out.compressibility.R2TopWallTotal);
    if isfield(out.compressibility, 'KeffPtot') && isfinite(out.compressibility.KeffPtot)
        fprintf('Keff Pvir / Ptot            : %.12g / %.12g  (fit R2 %.4g / %.4g)\n', ...
            out.compressibility.KeffPvir, out.compressibility.KeffPtot, ...
            out.compressibility.R2Pvir, out.compressibility.R2Ptot);
    end
end
methodLower = lower(strrep(char(string(params.method)), '-', '_'));
if any(strcmp(methodLower, {'q9','q9_full','q9_historical','q9_mass_only','mass_only'}))
    fprintf('mean/final mass-flux before : %.12g / %.12g\n', out.summary.meanMassFluxBefore, out.summary.finalMassFluxBefore);
    fprintf('mean/final mass-flux after  : %.12g / %.12g\n', out.summary.meanMassFluxAfter, out.summary.finalMassFluxAfter);
    fprintf('mean/final mass-flux resid  : %.12g / %.12g\n', out.summary.meanMassFluxResidual, out.summary.finalMassFluxResidual);
end
end

function params = parse_inputs(varargin)
if nargin == 1 && isstruct(varargin{1})
    params = varargin{1};
    return;
end
params = struct();
if mod(nargin, 2) ~= 0
    error('Use name-value pairs or a single struct.');
end
for k = 1:2:nargin
    params.(char(varargin{k})) = varargin{k+1};
end
end

function params = set_default_params(params)
params = set_default(params, 'method', 'classic');
params = set_default(params, 'Lx', 1.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Ly0', params.Ly);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 32);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 0.01);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.0);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'nSteps', 5000);
params = set_default(params, 'sampleEvery', 100);
params = set_default(params, 'progressEvery', 500);
params = set_default(params, 'seed', 11);
params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'initialMeanVelocityX', 0.0);
params = set_default(params, 'initialMeanVelocityY', 0.0);
params = set_default(params, 'thermostatAfterStep', true);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'thermostatMinParticlesPerCell', 2);
params = set_default(params, 'thermostatMaxScale', 5.0);
params = set_default(params, 'wallModeY', 'bounceback');
params = set_default(params, 'pistonTangentialMode', 'bounceback');
params = set_default(params, 'wallSigma', sqrt(params.kBT));
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'Vbottom', 0.0);
params = set_default(params, 'pistonY0', params.Ly);

% Compression aliases and precedence.
% Historical wrappers use `pistonCompressionTarget`, while some direct calls
% use `compressionTarget`.  In earlier versions, direct calls with only
% `compressionTarget` were ignored here, leaving the default 5% compression.
hasCompressionTarget = isfield(params, 'compressionTarget') && ~isempty(params.compressionTarget);
hasPistonCompressionTarget = isfield(params, 'pistonCompressionTarget') && ~isempty(params.pistonCompressionTarget);
hasPistonYMin = isfield(params, 'pistonYMin') && ~isempty(params.pistonYMin);

if hasCompressionTarget && ~hasPistonCompressionTarget
    params.pistonCompressionTarget = params.compressionTarget;
end
params = set_default(params, 'pistonCompressionTarget', 0.05);

if hasPistonYMin
    % Explicit pistonYMin is the most concrete geometric input.
    params.pistonCompressionTarget = 1.0 - params.pistonYMin / max(params.pistonY0, eps);
else
    params.pistonYMin = params.pistonY0 * (1.0 - params.pistonCompressionTarget);
end
params.compressionTarget = params.pistonCompressionTarget;

if ~isfield(params, 'pistonVy') || isempty(params.pistonVy)
    params.pistonVy = -(params.pistonY0 - params.pistonYMin) / max(params.nSteps * params.dt, eps);
end
params = set_default(params, 'pistonStopOnMin', true);
params = set_default(params, 'pistonAffineReposition', true);
params = set_default(params, 'pistonMotionMode', 'linear');
params = set_default(params, 'pistonStaircaseCompressionList', [0 0.0025 0.005 0.01]);
params = set_default(params, 'pistonStaircaseInitialHoldSteps', 0);
params = set_default(params, 'pistonStaircaseMoveSteps', 100);
params = set_default(params, 'pistonStaircaseHoldSteps', 400);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'wallVirtualParticlesEnable', true);
if isfield(params, 'densityFactor') && ~isfield(params, 'wallVirtualParticlesDensityFactor')
    params.wallVirtualParticlesDensityFactor = params.densityFactor;
end
params = set_default(params, 'wallVirtualParticlesGeometryMode', 'shifted_solid_fraction');
params = set_default(params, 'wallVirtualParticlesDensityFactor', 0.8);
params = set_default(params, 'wallVirtualParticlesThermal', true);
params = set_default(params, 'wallVirtualParticlesKBT', params.kBT);
params = set_default(params, 'wallVirtualParticlesForceRandomShiftY', true);
params = set_default(params, 'wallVirtualParticlesStochasticCount', true);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'projectionMomentumCorrectionEnable', true);
params = set_default(params, 'projectionMomentumCorrectionMode', 'particle_global_exact');
params = set_default(params, 'massFluxProjectionOperator', 'general_bc');
params = set_default(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk');
params = set_default(params, 'massFluxProjectionStrength', 1.0);
params = set_default(params, 'massFluxDensityRelaxationBeta', 5.0e-4);
params = set_default(params, 'massFluxTargetFilter', 'elliptic_lowpass');
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', true);
params = set_default(params, 'massFluxFinalVelocityProjectionStrength', 0.5);
params = set_default(params, 'massFluxLowKMaxIndex', 2);
params = set_default(params, 'lowKMaxIndex', 2);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', false);
params = set_default(params, 'massFluxEllipticUseSolveCache', true);
params = set_default(params, 'pistonQ9TargetMode', 'reference_gamma');
params = set_default(params, 'pistonQ9FreezeEllipticMetric', false);
params = set_default(params, 'pistonQ9Level2NoQ6Projection', false);
params = set_default(params, 'pistonDiagnosticsTopLayerCells', 2);
% Virial EOS/kick closure.  Disabled by default: diagnostics can be enabled
% without applying any velocity kick.  Kvirial and betaEOS are historical
% aliases kept for compatibility with the former liquid-closure scripts.
params = set_default(params, 'virialDiagnosticsEnable', false);
params = set_default(params, 'virialKickEnable', false);
params = set_default(params, 'Kvirial', 0.0);
params = set_default(params, 'virialBeta', get_param(params, 'betaEOS', 0.0));
params = set_default(params, 'virialDriveTargetMode', 'current_uniform');
params = set_default(params, 'virialRhoEOSRefMode', 'initial_physical_density');
params = set_default(params, 'virialRhoUniformMode', 'reference_gamma_current_volume');
params = set_default(params, 'virialRhoKickMode', 'uniform_now');
params = set_default(params, 'virialLimiterEnable', true);
params = set_default(params, 'virialMaxDuFractionThermal', 0.1);
params = set_default(params, 'virialStoreMaps', false);
params = set_default(params, 'pistonCompressibilityFitStartFraction', 0.2);
params = set_default(params, 'pistonCompressibilityFitEndFraction', 1.0);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'densityBandFraction', 0.20);
params = set_default(params, 'visualEnable', true);
params = set_default(params, 'visualEvery', 500);
params = set_default(params, 'visualPause', 0.001);
params = set_default(params, 'visualFigureId', 620);
params = set_default(params, 'visualMaxParticles', 5000);
params = set_default(params, 'maxWallClockSeconds', Inf);
end

function sd = piston_sample_diagnostics(state, params0, it, lastStepDiag)
params = active_params_from_state(params0, state);
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
thermal = projection_thermal_diagnostics(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
pop = projection_population_diagnostics(state.x, params, ...
    'periodicX', true, 'periodicY', false, 'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction);

Tmap = local_temperature_map(state.x, state.v, params, G);
PkinMap = G.rho .* max(Tmap, 0);
PkinMean = mean(PkinMap(isfinite(PkinMap)), 'omitnan');
if ~isfinite(PkinMean)
    PkinMean = (size(state.x,1) / max(params.Lx * params.Ly, eps)) * thermal.kBTCellRelative;
end
lowK = lowk_density_energy(G.N, params.lowKMaxIndex, params.gamma);

if isfield(state, 'piston')
    piston = state.piston;
else
    piston = struct('yTop', params.Ly, 'compression', 0, 'Up', 0, 'Ly0', params.Ly);
end
rhoPhysicalExpected = size(state.x, 1) / max(params.Lx * piston.yTop, eps);
rhoPhysicalMean = rhoPhysicalExpected;
rhoPhysicalRelError = (rhoPhysicalMean - rhoPhysicalExpected) / max(abs(rhoPhysicalExpected), eps);
gammaMeanExpected = size(state.x, 1) / max(params.Nx * params.Ny, 1);
gammaMeanRelError = (pop.meanN - gammaMeanExpected) / max(abs(gammaMeanExpected), eps);

layerCells = max(1, round(get_param(params, 'pistonDiagnosticsTopLayerCells', 2)));
botIds = 1:min(params.Ny, layerCells);
topIds = max(1, params.Ny-layerCells+1):params.Ny;
PkinBottomLayerMean = mean(PkinMap(:, botIds), 'all', 'omitnan');
PkinTopLayerMean = mean(PkinMap(:, topIds), 'all', 'omitnan');
PkinIdealMean = rhoPhysicalMean * thermal.kBTCellRelative;
PkinIdealRatio = PkinMean / max(PkinIdealMean, eps);

virialInfo = struct();
if isstruct(lastStepDiag)
    virialInfo = getf(lastStepDiag, 'virial', struct());
end
PvirMean = getf(virialInfo, 'PvirMean', NaN);
PtotMeanFromVirialInfo = getf(virialInfo, 'PtotMean', NaN);
PdriveMeanFromVirialInfo = getf(virialInfo, 'PdriveMean', NaN);

% The virial closure is applied inside mpcd_step_projection_piston after the
% Q6/Q9 block. When virialKickEnable=true, the kick and the final thermostat
% can modify the velocities after the virial pressure field has been built.
% Therefore virialInfo.PtotMean is a pre-kick/pre-final-thermostat diagnostic,
% whereas PkinMean above is recomputed here from the final state used for all
% piston diagnostics. The EOS total pressure reported in diagTable/summary
% must be evaluated in one consistent frame:
%
%     Ptot(final) = Pkin(final) + Pvir(rho final)
%
% Positions are unchanged by the velocity kick, so PvirMean is still valid at
% the end of the step. The pressure field used for the kick is only defined up
% to a spatially uniform offset: the historical/current_uniform drive subtracts
% rhoUniformNow, whereas the EOS pressure uses rhoEOSRef. These two choices
% have the same gradient and therefore the same velocity kick. For
% pressure-density diagnostics and compressibility fits, use the EOS gauge
% PdriveMean = PtotMean, and preserve the raw centered kick-field mean as
% PdriveMeanFromVirialInfo.
if isfinite(PvirMean)
    PtotMean = PkinMean + PvirMean;
else
    PtotMean = PtotMeanFromVirialInfo;
end
PdriveMean = PtotMean;
PtotIdealRatio = PtotMean / max(PkinIdealMean, eps);

massFluxResidual = NaN;
massFluxBefore = NaN;
massFluxAfter = NaN;
if isstruct(lastStepDiag)
    massFluxResidual = getf(lastStepDiag, 'rmsMassFluxDivResidual', NaN);
    massFluxBefore = getf(lastStepDiag, 'rmsMassFluxDivBefore', NaN);
    massFluxAfter = getf(lastStepDiag, 'rmsMassFluxDivProjectedAfter', NaN);
end
wallInfo = struct();
if isstruct(lastStepDiag)
    wallInfo = getf(lastStepDiag, 'wallInfo', struct());
end

row = table();
row.step = it;
row.t = it * params.dt;
row.yTop = piston.yTop;
row.Up = getf(piston, 'Up', NaN);
row.compression = getf(piston, 'compression', 1 - piston.yTop / params.Ly0);
row.pistonPhaseIndex = getf(piston, 'phaseIndex', NaN);
row.pistonPhaseName = string(getf(piston, 'phase', 'unknown'));
row.pistonUp = getf(piston, 'Up', NaN);
row.rhoPhysicalMean = rhoPhysicalMean;
row.rhoPhysicalExpected = rhoPhysicalExpected;
row.rhoPhysicalRelError = rhoPhysicalRelError;
row.gammaMeanCell = size(state.x, 1) / max(params.Nx * params.Ny, 1);
row.gammaMeanExpected = gammaMeanExpected;
row.gammaMeanRelError = gammaMeanRelError;
row.meanN = pop.meanN;
row.stdN = pop.stdN;
row.relStdN = pop.stdN / max(pop.meanN, eps);
row.outBandFraction = pop.outBandFraction;
row.lowKDensity = lowK;
row.PkinMean = PkinMean;
row.PkinIdealMean = PkinIdealMean;
row.PkinIdealRatio = PkinIdealRatio;
row.PkinBottomLayerMean = PkinBottomLayerMean;
row.PkinTopLayerMean = PkinTopLayerMean;
row.PvirMean = PvirMean;
row.PtotMean = PtotMean;
row.PtotMeanFromVirialInfo = PtotMeanFromVirialInfo;
row.PdriveMean = PdriveMean;
row.PdriveMeanFromVirialInfo = PdriveMeanFromVirialInfo;
row.PtotIdealRatio = PtotIdealRatio;
row.virialKvirial = getf(virialInfo, 'Kvirial', NaN);
row.virialBeta = getf(virialInfo, 'betaVirial', NaN);
row.virialRhoEOSRef = getf(virialInfo, 'rhoEOSRef', NaN);
row.virialRhoUniformNow = getf(virialInfo, 'rhoUniformNow', NaN);
row.virialRhoDefectRms = getf(virialInfo, 'rhoDefectRms', NaN);
row.virialRhoDefectRelRms = getf(virialInfo, 'rhoDefectRelRms', NaN);
row.virialGradPdriveRms = getf(virialInfo, 'gradPdriveRms', NaN);
row.virialDuRms = getf(virialInfo, 'duVirialAppliedRms', NaN);
row.virialDuMaxAbs = getf(virialInfo, 'duVirialAppliedMaxAbs', NaN);
row.virialDuOverThermalRms = getf(virialInfo, 'duVirialOverThermalRms', NaN);
row.virialLimitedCellFraction = getf(virialInfo, 'limitedCellFraction', NaN);
row.virialLimitedCellCount = getf(virialInfo, 'limitedCellCount', NaN);
row.virialLimiterDuMax = getf(virialInfo, 'duLimiterMax', NaN);
row.virialDuRmsBeforeLimiter = getf(virialInfo, 'duRmsBeforeLimiter', NaN);
row.virialDuRmsAfterLimiter = getf(virialInfo, 'duRmsAfterLimiter', NaN);
row.virialMaxDuFractionThermal = getf(virialInfo, 'duMaxFractionThermal', NaN);
row.virialRawMomentumKickNorm = getf(virialInfo, 'rawMomentumKickNorm', NaN);
row.virialResidualMomentumKickNorm = getf(virialInfo, 'residualMomentumKickNorm', NaN);
row.kBTCell = thermal.kBTCellRelative;
row.meanVx = mean(state.v(:,1), 'omitnan');
row.meanVy = mean(state.v(:,2), 'omitnan');
row.nTopHits = getf(wallInfo, 'nTop', NaN);
row.nBottomHits = getf(wallInfo, 'nBot', NaN);
row.pressureTopWall = getf(wallInfo, 'pressureTopWall', NaN);
row.pressureBottomWall = getf(wallInfo, 'pressureBottomWall', NaN);
row.pressureTopWallImpact = getf(wallInfo, 'pressureTopWallImpact', row.pressureTopWall);
row.pressureTopWallVP = getf(wallInfo, 'pressureTopWallVP', 0.0);
row.pressureTopWallTotal = getf(wallInfo, 'pressureTopWallTotal', row.pressureTopWallImpact + row.pressureTopWallVP);
row.pressureBottomWallImpact = getf(wallInfo, 'pressureBottomWallImpact', row.pressureBottomWall);
row.pressureBottomWallVP = getf(wallInfo, 'pressureBottomWallVP', 0.0);
row.pressureBottomWallTotal = getf(wallInfo, 'pressureBottomWallTotal', row.pressureBottomWallImpact + row.pressureBottomWallVP);
row.pistonPowerOnFluid = getf(wallInfo, 'pistonPowerOnFluid', NaN);
row.pistonWorkIncrement = getf(wallInfo, 'pistonWorkIncrement', NaN);
if isfield(piston, 'workOnFluidCumulative')
    row.pistonWorkOnFluidCumulative = piston.workOnFluidCumulative;
else
    row.pistonWorkOnFluidCumulative = getf(wallInfo, 'pistonWorkOnFluidCumulative', NaN);
end
row.massFluxResidual = massFluxResidual;
row.massFluxBefore = massFluxBefore;
row.massFluxAfter = massFluxAfter;

sd = struct();
sd.row = row;
sd.G = G;
sd.Tmap = Tmap;
sd.PkinMap = PkinMap;
sd.virial = virialInfo;
sd.UxProfile = mean(G.Ux, 1, 'omitnan').';
sd.UyProfile = mean(G.Uy, 1, 'omitnan').';
sd.NProfile = mean(G.N, 1, 'omitnan').';
end

function params = active_params_from_state(params0, state)
params = params0;
params.Ly0 = get_param(params0, 'Ly0', get_param(params0, 'LyReference', params0.Ly));
if isfield(state, 'piston') && isfield(state.piston, 'yTop') && ~isempty(state.piston.yTop)
    params.Ly = state.piston.yTop;
else
    params.Ly = params.pistonY0;
end
end

function T = local_temperature_map(x, v, params, G)
Nx = params.Nx;
Ny = params.Ny;
dx = params.Lx / Nx;
dy = params.Ly / Ny;
ix = floor(mod(x(:,1), params.Lx) / dx) + 1;
iy = floor(min(max(x(:,2), 0), params.Ly - eps(params.Ly)) / dy) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
cellId = iy + Ny*(ix-1);
UxVec = reshape(G.Ux.', [], 1);
UyVec = reshape(G.Uy.', [], 1);
NVec = reshape(G.N.', [], 1);
ux = UxVec(cellId);
uy = UyVec(cellId);
rel2 = (v(:,1)-ux).^2 + (v(:,2)-uy).^2;
sumRel2 = accumarray(cellId, rel2, [Nx*Ny, 1], @sum, 0);
dof = 2*max(double(NVec)-1, 0);
Tvec = nan(Nx*Ny,1);
mask = dof > 0;
Tvec(mask) = sumRel2(mask)./dof(mask);
T = reshape(Tvec, [Ny, Nx]).';
end

function e = lowk_density_energy(N, kmax, gamma)
D = double(N) - mean(double(N(:)), 'omitnan');
if ~isfinite(gamma) || gamma <= 0
    gamma = max(mean(double(N(:)), 'omitnan'), eps);
end
D = D ./ gamma;
[Nx, Ny] = size(D);
H = fft2(D);
ix = [0:floor(Nx/2), -ceil(Nx/2)+1:-1];
iy = [0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix) > Nx, ix = ix(1:Nx); end
if numel(iy) > Ny, iy = iy(1:Ny); end
[KX, KY] = ndgrid(ix, iy);
mask = sqrt(double(KX).^2 + double(KY).^2) <= double(kmax);
mask(1,1) = false;
H(~mask) = 0;
e = mean(abs(H(:)).^2, 'omitnan') / max(numel(H)^2, 1);
end

function T = rows_to_table(rows)
rows = rows(~cellfun(@isempty, rows));
if isempty(rows)
    T = table();
else
    T = vertcat(rows{:});
end
end

function summary = summarize_piston(out, Np0)
T = out.diagTable;
summary = struct();
summary.Np0 = Np0;
summary.NpFinal = size(out.state.x, 1);
summary.finalYTop = last_or_nan(T.yTop);
summary.finalCompression = last_or_nan(T.compression);
summary.finalRhoPhysicalMean = last_or_nan(T.rhoPhysicalMean);
summary.meanRhoPhysicalMean = mean(T.rhoPhysicalMean, 'omitnan');
summary.finalRhoPhysicalRelError = last_or_nan(T.rhoPhysicalRelError);
summary.meanRhoPhysicalRelError = mean(T.rhoPhysicalRelError, 'omitnan');
summary.finalStdN = last_or_nan(T.stdN);
summary.meanStdN = mean(T.stdN, 'omitnan');
summary.finalLowKDensity = last_or_nan(T.lowKDensity);
summary.meanLowKDensity = mean(T.lowKDensity, 'omitnan');
summary.finalPkinMean = last_or_nan(T.PkinMean);
summary.meanPkinMean = mean(T.PkinMean, 'omitnan');
summary.finalPkinIdealMean = last_or_nan(T.PkinIdealMean);
summary.meanPkinIdealMean = mean(T.PkinIdealMean, 'omitnan');
summary.finalPkinIdealRatio = last_or_nan(T.PkinIdealRatio);
summary.meanPkinIdealRatio = mean(T.PkinIdealRatio, 'omitnan');
summary.finalPkinBottomLayerMean = last_or_nan(T.PkinBottomLayerMean);
summary.meanPkinBottomLayerMean = mean(T.PkinBottomLayerMean, 'omitnan');
summary.finalPkinTopLayerMean = last_or_nan(T.PkinTopLayerMean);
summary.meanPkinTopLayerMean = mean(T.PkinTopLayerMean, 'omitnan');
summary.finalPvirMean = last_or_nan(column_or_nan(T, 'PvirMean'));
summary.meanPvirMean = mean(column_or_nan(T, 'PvirMean'), 'omitnan');
summary.finalPtotMean = last_or_nan(column_or_nan(T, 'PtotMean'));
summary.meanPtotMean = mean(column_or_nan(T, 'PtotMean'), 'omitnan');
summary.finalPtotMeanFromVirialInfo = last_or_nan(column_or_nan(T, 'PtotMeanFromVirialInfo'));
summary.meanPtotMeanFromVirialInfo = mean(column_or_nan(T, 'PtotMeanFromVirialInfo'), 'omitnan');
summary.finalPdriveMean = last_or_nan(column_or_nan(T, 'PdriveMean'));
summary.meanPdriveMean = mean(column_or_nan(T, 'PdriveMean'), 'omitnan');
summary.finalPdriveMeanFromVirialInfo = last_or_nan(column_or_nan(T, 'PdriveMeanFromVirialInfo'));
summary.meanPdriveMeanFromVirialInfo = mean(column_or_nan(T, 'PdriveMeanFromVirialInfo'), 'omitnan');
summary.finalPtotIdealRatio = last_or_nan(column_or_nan(T, 'PtotIdealRatio'));
summary.meanPtotIdealRatio = mean(column_or_nan(T, 'PtotIdealRatio'), 'omitnan');
summary.finalVirialDuRms = last_or_nan(column_or_nan(T, 'virialDuRms'));
summary.meanVirialDuRms = mean(column_or_nan(T, 'virialDuRms'), 'omitnan');
summary.finalVirialDuOverThermalRms = last_or_nan(column_or_nan(T, 'virialDuOverThermalRms'));
summary.meanVirialDuOverThermalRms = mean(column_or_nan(T, 'virialDuOverThermalRms'), 'omitnan');
summary.finalVirialLimitedCellFraction = last_or_nan(column_or_nan(T, 'virialLimitedCellFraction'));
summary.meanVirialLimitedCellFraction = mean(column_or_nan(T, 'virialLimitedCellFraction'), 'omitnan');
summary.maxVirialLimitedCellFraction = max(column_or_nan(T, 'virialLimitedCellFraction'), [], 'omitnan');
summary.finalVirialLimitedCellCount = last_or_nan(column_or_nan(T, 'virialLimitedCellCount'));
summary.maxVirialLimitedCellCount = max(column_or_nan(T, 'virialLimitedCellCount'), [], 'omitnan');
summary.finalVirialLimiterDuMax = last_or_nan(column_or_nan(T, 'virialLimiterDuMax'));
summary.meanVirialDuRmsBeforeLimiter = mean(column_or_nan(T, 'virialDuRmsBeforeLimiter'), 'omitnan');
summary.meanVirialDuRmsAfterLimiter = mean(column_or_nan(T, 'virialDuRmsAfterLimiter'), 'omitnan');
summary.finalVirialMaxDuFractionThermal = last_or_nan(column_or_nan(T, 'virialMaxDuFractionThermal'));
summary.finalVirialResidualMomentumKickNorm = last_or_nan(column_or_nan(T, 'virialResidualMomentumKickNorm'));
summary.meanVirialResidualMomentumKickNorm = mean(column_or_nan(T, 'virialResidualMomentumKickNorm'), 'omitnan');
summary.finalKBTCell = last_or_nan(T.kBTCell);
summary.meanKBTCell = mean(T.kBTCell, 'omitnan');
summary.meanTopWallPressure = mean(T.pressureTopWall, 'omitnan');
summary.finalTopWallPressureImpact = last_or_nan(T.pressureTopWallImpact);
summary.meanTopWallPressureImpact = mean(T.pressureTopWallImpact, 'omitnan');
summary.finalTopWallPressureVP = last_or_nan(T.pressureTopWallVP);
summary.meanTopWallPressureVP = mean(T.pressureTopWallVP, 'omitnan');
summary.finalTopWallPressureTotal = last_or_nan(T.pressureTopWallTotal);
summary.meanTopWallPressureTotal = mean(T.pressureTopWallTotal, 'omitnan');
summary.finalPistonPowerOnFluid = last_or_nan(T.pistonPowerOnFluid);
summary.meanPistonPowerOnFluid = mean(T.pistonPowerOnFluid, 'omitnan');
summary.finalPistonWorkOnFluidCumulative = last_or_nan(T.pistonWorkOnFluidCumulative);
summary.meanPistonWorkOnFluidCumulative = mean(T.pistonWorkOnFluidCumulative, 'omitnan');
summary.finalMassFluxBefore = last_or_nan(T.massFluxBefore);
summary.meanMassFluxBefore = mean(T.massFluxBefore, 'omitnan');
summary.finalMassFluxAfter = last_or_nan(T.massFluxAfter);
summary.meanMassFluxAfter = mean(T.massFluxAfter, 'omitnan');
summary.finalMassFluxResidual = last_or_nan(T.massFluxResidual);
summary.meanMassFluxResidual = mean(T.massFluxResidual, 'omitnan');
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
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

function x = column_or_nan(T, name)
if istable(T) && ismember(name, T.Properties.VariableNames)
    x = T.(name);
else
    x = nan(height(T), 1);
end
end

function v = last_or_nan(x)
if isempty(x)
    v = NaN;
else
    idx = find(isfinite(x), 1, 'last');
    if isempty(idx)
        v = NaN;
    else
        v = x(idx);
    end
end
end
