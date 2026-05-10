function result = test_projection_exact_cell_initialization()
%TEST_PROJECTION_EXACT_CELL_INITIALIZATION Verify exact-per-cell initialization.

params = struct();
params.Lx = 1.5;
params.Ly = 1.0;
params.Nx = 12;
params.Ny = 8;
params.gamma = 6;
params.kBT = 1.0;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;

rng(123);
[state, info] = projection_initialize_particles(params);
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 0);

passCount = size(state.x, 1) == params.Nx * params.Ny * params.gamma;
passExact = all(G.N == params.gamma, 'all');
passStd = info.initialStdN == 0;
passDomain = all(state.x(:,1) >= 0 & state.x(:,1) < params.Lx) && ...
    all(state.x(:,2) >= 0 & state.x(:,2) <= params.Ly);
passVelocityFinite = all(isfinite(state.v(:))) && ...
    abs(mean(state.v(:,1))) < 1e-12 && abs(mean(state.v(:,2))) < 1e-12;
passed = passCount && passExact && passStd && passDomain && passVelocityFinite;

result = struct();
result.info = info;
result.passCount = passCount;
result.passExact = passExact;
result.passStd = passStd;
result.passDomain = passDomain;
result.passVelocityFinite = passVelocityFinite;
result.passed = passed;

fprintf('\n=== test_projection_exact_cell_initialization ===\n');
fprintf('mode                 : %s\n', info.mode);
fprintf('Np                   : %d\n', info.Np);
fprintf('initial std(N)       : %.12g\n', info.initialStdN);
fprintf('initial max |N-gamma|: %.12g\n', info.initialMaxAbsDeltaFromGamma);
fprintf('initial out-band20   : %.12g\n', info.initialOutBand20);
fprintf('passCount            : %d\n', passCount);
fprintf('passExact            : %d\n', passExact);
fprintf('passStd              : %d\n', passStd);
fprintf('passDomain           : %d\n', passDomain);
fprintf('passVelocityFinite   : %d\n', passVelocityFinite);
fprintf('passed               : %d\n', passed);

if ~passed
    error('test_projection_exact_cell_initialization failed.');
end
end
