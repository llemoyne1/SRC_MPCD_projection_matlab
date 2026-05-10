function result = test_projection_periodic_particles()
%TEST_PROJECTION_PERIODIC_PARTICLES Test particle-grid pressure-projection coupling.
%
% This validates the second prototype layer after the grid-only projection:
% a classic periodic MPCD step, particle-to-grid deposit, FFT projection, and
% grid-to-particle velocity correction.

params = struct();
params.Lx = 1.0;
params.Ly = 1.0;
params.Nx = 32;
params.Ny = 32;
params.gamma = 20;
params.dt = 1.0e-3;
params.kBT = 1.0;
params.alphaDeg = 90;
params.useRandomGridShift = true;
params.bodyForceX = 0.0;
params.bodyForceY = 0.0;
params.projectionEnable = true;
params.projectionStrength = 1.0;
params.projectionInterpolationMethod = 'nearest';

rng(42);
Np = round(params.gamma * params.Nx * params.Ny);
state0 = struct();
state0.x = [params.Lx * rand(Np, 1), params.Ly * rand(Np, 1)];
state0.v = sqrt(params.kBT) * randn(Np, 2);
state0.v = state0.v - mean(state0.v, 1);
state0.v(:, 1) = state0.v(:, 1) + 0.25 * sin(2*pi*state0.x(:, 1)/params.Lx);
state0.v(:, 2) = state0.v(:, 2) + 0.25 * sin(2*pi*state0.x(:, 2)/params.Ly);

[state1, diag1] = mpcd_step_projection_periodic(state0, params);

% With nearest correction and no empty cells, the particle re-deposited field
% should match the projected grid field to round-off. We keep a tolerant bound
% to allow rare empty cells in random initializations.
passProjected = diag1.rmsDivProjectedAfter < 1e-10 * max(diag1.rmsDivBefore, 1);
passParticle = diag1.rmsDivParticleAfter < 1e-8 * max(diag1.rmsDivBefore, 1);
passFinite = all(isfinite(state1.x(:))) && all(isfinite(state1.v(:))) && isfinite(diag1.kBTAfterProjection);
passDomain = all(state1.x(:, 1) >= 0 & state1.x(:, 1) < params.Lx) && ...
             all(state1.x(:, 2) >= 0 & state1.x(:, 2) < params.Ly);

% projectionStrength = 0 must reduce exactly to the classic step if the RNG
% is reset, because no grid correction is applied.
params0 = params;
params0.projectionStrength = 0.0;
rng(123);
[stateClassic, ~] = mpcd_step_classic_periodic(state0, params0);
rng(123);
[stateProj0, diag0] = mpcd_step_projection_periodic(state0, params0);
passStrengthZero = max(abs(stateClassic.x(:) - stateProj0.x(:))) == 0 && ...
                   max(abs(stateClassic.v(:) - stateProj0.v(:))) == 0 && ...
                   diag0.dvRms == 0;

result = struct();
result.rmsDivBefore = diag1.rmsDivBefore;
result.rmsDivProjectedAfter = diag1.rmsDivProjectedAfter;
result.rmsDivParticleAfter = diag1.rmsDivParticleAfter;
result.divReductionProjected = diag1.divReductionProjected;
result.divReductionParticle = diag1.divReductionParticle;
result.nEmptyCellsAfter = diag1.nEmptyCellsAfter;
result.kBTBeforeProjection = diag1.kBTBeforeProjection;
result.kBTAfterProjection = diag1.kBTAfterProjection;
result.passProjected = passProjected;
result.passParticle = passParticle;
result.passFinite = passFinite;
result.passDomain = passDomain;
result.passStrengthZero = passStrengthZero;
result.passed = passProjected && passParticle && passFinite && passDomain && passStrengthZero;

fprintf('\n=== test_projection_periodic_particles ===\n');
fprintf('rms div before         : %.12e\n', result.rmsDivBefore);
fprintf('rms div projected      : %.12e\n', result.rmsDivProjectedAfter);
fprintf('rms div particle grid  : %.12e\n', result.rmsDivParticleAfter);
fprintf('reduction projected    : %.12e\n', result.divReductionProjected);
fprintf('reduction particle     : %.12e\n', result.divReductionParticle);
fprintf('nEmptyCellsAfter       : %d\n', result.nEmptyCellsAfter);
fprintf('kBT before projection  : %.12e\n', result.kBTBeforeProjection);
fprintf('kBT after projection   : %.12e\n', result.kBTAfterProjection);
fprintf('passProjected          : %d\n', passProjected);
fprintf('passParticle           : %d\n', passParticle);
fprintf('passFinite             : %d\n', passFinite);
fprintf('passDomain             : %d\n', passDomain);
fprintf('passStrengthZero       : %d\n', passStrengthZero);
fprintf('passed                 : %d\n', result.passed);

if ~result.passed
    error('test_projection_periodic_particles failed.');
end
end
