function result = test_projection_population_transport_diagnostics()
%TEST_PROJECTION_POPULATION_TRANSPORT_DIAGNOSTICS Validate density transport diagnostics.
%
% This test checks the distinction between:
%   1. instantaneous population, unchanged by pressure projection, and
%   2. one-step forecast population, changed by the velocity correction.

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
params.bodyForceY = 0.0;
params.projectionEnable = true;
params.projectionStrength = 1.0;
params.projectionInterpolationMethod = 'nearest';
params.projectionTransportDiagnosticsEnable = true;
params.wallModeY = 'bounceback';
params.useRandomGridShiftX = true;
params.useRandomGridShiftY = false;

rng(11);
Np = round(params.gamma * params.Nx * params.Ny);
state = struct();
state.x = [params.Lx * rand(Np, 1), params.Ly * rand(Np, 1)];
state.v = sqrt(params.kBT) * randn(Np, 2);
state.v = state.v - mean(state.v, 1);

[stateOut, diag] = mpcd_step_projection_poiseuille(state, params); %#ok<ASGLU>
tr = diag.populationTransport;

result = struct();
result.maxPopDeltaProjection = diag.populationProjectionDeltaMaxAbs;
result.transportClassicRms = tr.classic.rms;
result.transportProjectedRms = tr.projected.rms;
result.transportProjectedMinusClassicRms = tr.projectedMinusClassic.rms;
result.transportProjectedMinusClassicMax = tr.projectedMinusClassic.maxAbs;
result.massErrorClassic = tr.massErrorClassic;
result.massErrorProjected = tr.massErrorProjected;
result.massErrorProjectedMinusClassic = tr.massErrorProjectedMinusClassic;
result.passPopulationInstant = diag.populationProjectionDeltaMaxAbs == 0;
result.passTransportFinite = all(isfinite([result.transportClassicRms, ...
    result.transportProjectedRms, result.transportProjectedMinusClassicRms, ...
    result.transportProjectedMinusClassicMax]));
result.passMass = result.massErrorClassic == 0 && result.massErrorProjected == 0 && ...
    result.massErrorProjectedMinusClassic == 0;
result.passNonNegative = result.transportClassicRms >= 0 && result.transportProjectedRms >= 0 && ...
    result.transportProjectedMinusClassicRms >= 0;
result.passed = result.passPopulationInstant && result.passTransportFinite && ...
    result.passMass && result.passNonNegative;

fprintf('\n=== test_projection_population_transport_diagnostics ===\n');
fprintf('instant pop delta max      : %.12e\n', result.maxPopDeltaProjection);
fprintf('forecast rms classic       : %.12e\n', result.transportClassicRms);
fprintf('forecast rms projected     : %.12e\n', result.transportProjectedRms);
fprintf('forecast rms proj-classic  : %.12e\n', result.transportProjectedMinusClassicRms);
fprintf('forecast max proj-classic  : %.12e\n', result.transportProjectedMinusClassicMax);
fprintf('mass errors                : %.12e / %.12e / %.12e\n', ...
    result.massErrorClassic, result.massErrorProjected, result.massErrorProjectedMinusClassic);
fprintf('passPopulationInstant      : %d\n', result.passPopulationInstant);
fprintf('passTransportFinite        : %d\n', result.passTransportFinite);
fprintf('passMass                   : %d\n', result.passMass);
fprintf('passNonNegative            : %d\n', result.passNonNegative);
fprintf('passed                     : %d\n', result.passed);

if ~result.passed
    error('test_projection_population_transport_diagnostics failed.');
end
end
