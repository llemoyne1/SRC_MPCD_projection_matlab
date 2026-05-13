function [stateOut, diag] = mpcd_step_projection_piston(state, params, stepIndex)
%MPCD_STEP_PROJECTION_PISTON Classic piston MPCD step plus Q6/Q9 projection.
%
%   [stateOut, diag] = mpcd_step_projection_piston(state, params, stepIndex)
%
% The piston kinematics are updated from stepIndex, then the classic piston
% step is performed on the active height.  The Q6/Q9 projection is applied
% on the same active domain by temporarily setting params.Ly = pistonYCurrent.

if nargin < 3 || isempty(stepIndex)
    stepIndex = get_param(params, 'pistonStepIndex', 1);
end

[stepParams, pistonInfo] = piston_runtime_params(params, state, stepIndex);

[stateClassic, classicDiag] = mpcd_step_classic_piston(state, stepParams);

projParams = stepParams;
projParams.Ly0 = pistonInfo.Ly0;
projParams.Ly = pistonInfo.yTop;

[stateOut, projectionDiag] = mpcd_apply_q9_projection_channel(stateClassic, projParams);
stateOut.piston = pistonInfo;

diag = projectionDiag;
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.piston = pistonInfo;
diag.activeProjectionParams = struct('Lx', projParams.Lx, 'Ly', projParams.Ly, ...
    'Nx', projParams.Nx, 'Ny', projParams.Ny, 'gamma', get_param(projParams, 'gamma', NaN));
end

function [p, info] = piston_runtime_params(params, state, stepIndex)
p = params;
Ly0 = get_param(params, 'Ly0', get_param(params, 'LyReference', params.Ly));
y0 = get_param(params, 'pistonY0', Ly0);
yMin = get_param(params, 'pistonYMin', 0.9 * Ly0);
vy = get_param(params, 'pistonVy', 0.0);
dt = params.dt;
stopOnMin = logical(get_param(params, 'pistonStopOnMin', true));
motionMode = lower(char(string(get_param(params, 'pistonMotionMode', 'linear'))));

if strcmp(motionMode, 'cycle') || strcmp(motionMode, 'compression_hold_decompression')
    [yPrev, yTop, up, phaseName, phaseIndex] = piston_cycle_kinematics(params, stepIndex, Ly0, y0, yMin, dt);
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


function [yPrev, yTop, up, phaseName, phaseIndex] = piston_cycle_kinematics(params, stepIndex, Ly0, y0, yMin, dt)
%PISTON_CYCLE_KINEMATICS Slow compression-hold-decompression cycle.
% The profile is defined at integer step indices. yPrev corresponds to
% stepIndex-1 and yTop to stepIndex, so Up=(yTop-yPrev)/dt is exact for the
% piecewise-linear ramp used during compression/decompression.

nC = max(1, round(get_param(params, 'pistonCycleCompressSteps', get_param(params, 'pistonCompressSteps', 10000))));
nH = max(0, round(get_param(params, 'pistonCycleHoldSteps', get_param(params, 'pistonHoldSteps', 0))));
nD = max(1, round(get_param(params, 'pistonCycleDecompressSteps', get_param(params, 'pistonDecompressSteps', nC))));
nF = max(0, round(get_param(params, 'pistonCycleFinalHoldSteps', get_param(params, 'pistonFinalHoldSteps', 0))));

k0 = max(stepIndex - 1, 0);
k1 = max(stepIndex, 0);
yPrev = cycle_position(k0, y0, yMin, nC, nH, nD, nF);
yTop = cycle_position(k1, y0, yMin, nC, nH, nD, nF);
up = (yTop - yPrev) / max(dt, eps);
[phaseName, phaseIndex] = cycle_phase(k1, nC, nH, nD, nF);

% Numerical guard.
yPrev = min(max(yPrev, eps(Ly0)), Ly0);
yTop = min(max(yTop, eps(Ly0)), Ly0);
end

function y = cycle_position(k, y0, yMin, nC, nH, nD, nF)
if k <= nC
    y = y0 + (yMin - y0) * (k / max(nC, 1));
elseif k <= nC + nH
    y = yMin;
elseif k <= nC + nH + nD
    q = (k - nC - nH) / max(nD, 1);
    y = yMin + (y0 - yMin) * q;
else
    %#ok<NASGU> nF is kept to document the intended profile length.
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

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
