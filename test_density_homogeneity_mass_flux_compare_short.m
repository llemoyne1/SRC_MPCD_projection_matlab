function result = test_density_homogeneity_mass_flux_compare_short()
%TEST_DENSITY_HOMOGENEITY_MASS_FLUX_COMPARE_SHORT Short Q8 integration test.

params = struct();
params.Nx = 12;
params.Ny = 8;
params.gamma = 6;
params.nSteps = 60;
params.sampleEvery = 20;
params.progressEvery = 0;
params.bodyForceX = 0.02;
params.wallModeY = 'thermalize';
params.thermostatAfterProjection = true;
params.seed = 321;
params.initialPopulationMode = 'exact_per_cell';
params.projectedStrength = 1.0;
params.massFluxProjectionMode = 'conservative';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.0;
params.makeFigures = false;
params.saveFigures = false;
params.writeSummary = false;
params.maxWallClockSeconds = Inf;

cmp = run_compare_density_projection_mass_flux_poiseuille(params);
S = cmp.summary;

result = struct();
result.initialStdClassic = S.initialStdClassic;
result.initialStdProjected = S.initialStdProjected;
result.initialStdMassFlux = S.initialStdMassFlux;
result.meanStdClassic = S.meanStdClassic;
result.meanStdProjected = S.meanStdProjected;
result.meanStdMassFlux = S.meanStdMassFlux;
result.rhoTransportClassic = S.meanRhoTransportClassic;
result.rhoTransportProjected = S.meanRhoTransportProjected;
result.rhoTransportMassFlux = S.meanRhoTransportMassFlux;
result.massFluxResidual = S.meanMassFluxDivResidualMassFlux;
result.passFinite = all(isfinite([S.meanStdClassic, S.meanStdProjected, S.meanStdMassFlux, ...
    S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.meanRhoTransportMassFlux]));
result.passMaps = ~isempty(cmp.classic.NMaps) && ~isempty(cmp.projected.NMaps) && ~isempty(cmp.massFlux.NMaps);
result.passSameInitial = S.initialStdClassic == S.initialStdProjected && S.initialStdProjected == S.initialStdMassFlux;
result.passExactInitial = S.initialStdClassic == 0;
result.passMassFluxRan = isfinite(S.meanMassFluxDivResidualMassFlux);
result.passTransportFinite = isfinite(S.meanRhoTransportMassFlux);
result.passed = result.passFinite && result.passMaps && result.passSameInitial && ...
    result.passExactInitial && result.passMassFluxRan && result.passTransportFinite;

fprintf('\n=== test_density_homogeneity_mass_flux_compare_short ===\n');
fprintf('initial std C/P/M      : %.6g / %.6g / %.6g\n', ...
    result.initialStdClassic, result.initialStdProjected, result.initialStdMassFlux);
fprintf('mean std C/P/M         : %.6g / %.6g / %.6g\n', ...
    result.meanStdClassic, result.meanStdProjected, result.meanStdMassFlux);
fprintf('rho transport C/P/M    : %.6g / %.6g / %.6g\n', ...
    result.rhoTransportClassic, result.rhoTransportProjected, result.rhoTransportMassFlux);
fprintf('mass-flux residual M   : %.6g\n', result.massFluxResidual);
fprintf('passed                 : %d\n', result.passed);

if ~result.passed
    error('test_density_homogeneity_mass_flux_compare_short failed.');
end
end
