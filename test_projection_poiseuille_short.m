function result = test_projection_poiseuille_short()
%TEST_PROJECTION_POISEUILLE_SHORT Short integration test for Q3 Poiseuille.

params = struct();
params.Lx = 1.5;
params.Ly = 1.0;
params.Nx = 24;
params.Ny = 12;
params.gamma = 8;
params.dt = 1e-3;
params.kBT = 1.0;
params.alphaDeg = 90;
params.bodyForceX = 0.02;
params.nSteps = 80;
params.sampleEvery = 10;
params.projectionEnable = true;
params.projectionStrength = 1.0;
params.projectionInterpolationMethod = 'nearest';
params.wallModeY = 'bounceback';
params.makeFigures = false;
params.seed = 4;

out = run_projection_poiseuille_demo(params);
H = out.diagHistory;
state = out.state;

lastDivBefore = H(end, 3);
lastDivParticle = H(end, 5);
lastReduction = lastDivParticle / max(lastDivBefore, eps);
maxPopDeltaProjection = max(H(:, 14));

result = struct();
result.lastDivBefore = lastDivBefore;
result.lastDivParticle = lastDivParticle;
result.lastReduction = lastReduction;
result.maxPopDeltaProjection = maxPopDeltaProjection;
result.nuEff = out.viscosity.nuEff;
result.R2 = out.viscosity.R2;
result.passFinite = all(isfinite(state.x(:))) && all(isfinite(state.v(:))) && all(isfinite(H(:)));
result.passDomain = all(state.x(:,1) >= 0 & state.x(:,1) < params.Lx) && all(state.x(:,2) >= 0 & state.x(:,2) <= params.Ly);
result.passProjection = lastReduction < 1e-10;
result.passPopulationUnchangedByProjection = maxPopDeltaProjection == 0;
result.passKBT = H(end, 7) > 0 && H(end, 7) < 10 * params.kBT;
result.passViscosityFinite = isfinite(out.viscosity.nuEff);
result.passed = result.passFinite && result.passDomain && result.passProjection && ...
    result.passPopulationUnchangedByProjection && result.passKBT && result.passViscosityFinite;

fprintf('\n=== test_projection_poiseuille_short ===\n');
fprintf('last rms div before        : %.12e\n', result.lastDivBefore);
fprintf('last rms div particle      : %.12e\n', result.lastDivParticle);
fprintf('last reduction             : %.12e\n', result.lastReduction);
fprintf('max pop delta projection   : %.12e\n', result.maxPopDeltaProjection);
fprintf('nu_eff                     : %.12e\n', result.nuEff);
fprintf('R2                         : %.6f\n', result.R2);
fprintf('passFinite                 : %d\n', result.passFinite);
fprintf('passDomain                 : %d\n', result.passDomain);
fprintf('passProjection             : %d\n', result.passProjection);
fprintf('passPopulationNoChange     : %d\n', result.passPopulationUnchangedByProjection);
fprintf('passKBT                    : %d\n', result.passKBT);
fprintf('passViscosityFinite        : %d\n', result.passViscosityFinite);
fprintf('passed                     : %d\n', result.passed);

if ~result.passed
    error('test_projection_poiseuille_short failed.');
end
end
