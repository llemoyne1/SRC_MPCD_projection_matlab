function out = run_projection_poiseuille_long_demo(params)
%RUN_PROJECTION_POISEUILLE_LONG_DEMO Longer Poiseuille validation run.
%
%   out = run_projection_poiseuille_long_demo(params)
%
% Uses the Q4 diagnostics: cell-relative thermal energy, optional local
% thermostat after projection, continuous density-transport proxy and a
% late-time Poiseuille fit. This wrapper intentionally keeps the numerical
% core in run_projection_poiseuille_demo.

if nargin < 1 || isempty(params)
    params = struct();
end

params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 48);
params = set_default(params, 'Ny', 24);
params = set_default(params, 'gamma', 16);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 1.0);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 2.0e-2);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'nSteps', 2500);
params = set_default(params, 'sampleEvery', 25);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'projectionTransportDiagnosticsEnable', true);
params = set_default(params, 'densityTransportDiagnosticsEnable', true);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'thermostatMinParticlesPerCell', 2);
params = set_default(params, 'thermostatMaxScale', 5.0);
params = set_default(params, 'wallModeY', 'thermalize');
params = set_default(params, 'wallSigma', sqrt(params.kBT));
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'seed', 4);
params = set_default(params, 'excludeWallCellsFit', 2);
params = set_default(params, 'fitStartFraction', 0.5);
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'progressEvery', 250);
params = set_default(params, 'maxWallClockSeconds', 180);

out = run_projection_poiseuille_demo(params);

fprintf('\n=== physical Poiseuille checks ===\n');
fprintf('fit window t              : [%.6g, %.6g]\n', out.viscosity.tMin, out.viscosity.tMax);
fprintf('curvature a2              : %.12e\n', out.viscosity.a2);
fprintf('curvature sign ok         : %d\n', out.viscosity.hasCorrectCurvature);
fprintf('center minus wall U       : %.12e\n', out.viscosity.centerMinusWall);
fprintf('profile noise rms         : %.12e\n', out.viscosity.profileNoiseRms);
fprintf('signal/noise estimate     : %.6g\n', out.viscosity.signalToNoise);
fprintf('physical candidate        : %d\n', out.viscosity.physicalCandidate);
fprintf('mean cell kBT after       : %.12e\n', out.summary.meanKBTCellAfter);
fprintf('mean rho transport classic/projected/diff: %.6g / %.6g / %.6g\n', ...
    out.summary.meanDensityTransportClassicRms, ...
    out.summary.meanDensityTransportProjectedRms, ...
    out.summary.meanDensityTransportProjectedMinusClassicRms);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end
