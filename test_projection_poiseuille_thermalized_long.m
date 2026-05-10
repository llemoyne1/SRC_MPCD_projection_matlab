function result = test_projection_poiseuille_thermalized_long()
%TEST_PROJECTION_POISEUILLE_THERMALIZED_LONG Longer sanity test for Q4 diagnostics.
%
% This is not meant to certify final physical viscosity. It checks that a
% longer run with thermostat-after-projection stays finite, projected,
% thermally controlled, and population-preserving.

params = struct();
params.Nx = 32;
params.Ny = 16;
params.gamma = 10;
params.nSteps = 800;
params.sampleEvery = 20;
params.bodyForceX = 2e-2;
params.progressEvery = 0;
params.maxWallClockSeconds = 90;
params.makeFigures = false;
params.thermostatAfterProjection = true;
params.wallModeY = 'thermalize';
params.fitStartFraction = 0.5;
params.seed = 7;

out = run_projection_poiseuille_long_demo(params);

lastReduction = out.summary.lastRmsDivParticle / max(out.summary.lastRmsDivBefore, eps);
kBTTarget = out.params.kBT;
kBTCell = out.summary.lastKBTCellAfter;

result = struct();
result.lastReduction = lastReduction;
result.lastKBTCellAfter = kBTCell;
result.nuEff = out.viscosity.nuEff;
result.R2 = out.viscosity.R2;
result.hasCorrectCurvature = out.viscosity.hasCorrectCurvature;
result.signalToNoise = out.viscosity.signalToNoise;
result.physicalCandidate = out.viscosity.physicalCandidate;
result.maxPopDeltaProjection = out.summary.maxPopDeltaProjection;
result.passFinite = all(isfinite(out.diagHistory(:))) && isfinite(out.viscosity.nuEff);
result.passProjection = lastReduction < 1e-8;
result.passPopulationNoChange = out.summary.maxPopDeltaProjection == 0;
result.passKBTCell = isfinite(kBTCell) && kBTCell > 0.25*kBTTarget && kBTCell < 4.0*kBTTarget;
result.passDensityTransportFinite = isfinite(out.summary.lastDensityTransportProjectedRms);
result.passed = result.passFinite && result.passProjection && result.passPopulationNoChange && ...
    result.passKBTCell && result.passDensityTransportFinite;

fprintf('\n=== test_projection_poiseuille_thermalized_long ===\n');
fprintf('last reduction             : %.12e\n', result.lastReduction);
fprintf('last kBT cell after        : %.12e\n', result.lastKBTCellAfter);
fprintf('max pop delta projection   : %.12e\n', result.maxPopDeltaProjection);
fprintf('nu_eff                     : %.12e\n', result.nuEff);
fprintf('R2                         : %.6f\n', result.R2);
fprintf('curvature sign ok          : %d\n', result.hasCorrectCurvature);
fprintf('signal/noise estimate      : %.6g\n', result.signalToNoise);
fprintf('physical candidate         : %d\n', result.physicalCandidate);
fprintf('passFinite                 : %d\n', result.passFinite);
fprintf('passProjection             : %d\n', result.passProjection);
fprintf('passPopulationNoChange     : %d\n', result.passPopulationNoChange);
fprintf('passKBTCell                : %d\n', result.passKBTCell);
fprintf('passDensityTransportFinite : %d\n', result.passDensityTransportFinite);
fprintf('passed                     : %d\n', result.passed);

if ~result.passed
    error('test_projection_poiseuille_thermalized_long failed.');
end
end
