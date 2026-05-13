function result = test_density_homogeneity_lowk_mass_flux_compare_short()
%TEST_DENSITY_HOMOGENEITY_LOWK_MASS_FLUX_COMPARE_SHORT Short Q9 integration test.

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
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
params.makeFigures = false;
params.saveFigures = false;
params.writeSummary = false;
params.maxWallClockSeconds = Inf;

cmp = run_compare_density_projection_lowk_mass_flux_poiseuille(params);
S = cmp.summary;

result = struct();
result.initialStdClassic = S.initialStdClassic;
result.initialStdProjected = S.initialStdProjected;
result.initialStdHybrid = S.initialStdMassFlux;
result.meanStdClassic = S.meanStdClassic;
result.meanStdProjected = S.meanStdProjected;
result.meanStdHybrid = S.meanStdMassFlux;
result.lowKClassic = S.meanLowKClassic;
result.lowKProjected = S.meanLowKProjected;
result.lowKHybrid = S.meanLowKMassFlux;
result.massFluxResidual = S.meanMassFluxDivResidualMassFlux;
result.passFinite = all(isfinite([S.meanStdClassic, S.meanStdProjected, S.meanStdMassFlux, ...
    S.meanLowKClassic, S.meanLowKProjected, S.meanLowKMassFlux, S.meanMassFluxDivResidualMassFlux]));
result.passMaps = ~isempty(cmp.classic.NMaps) && ~isempty(cmp.projected.NMaps) && ~isempty(cmp.massFlux.NMaps);
result.passSameInitial = S.initialStdClassic == S.initialStdProjected && S.initialStdProjected == S.initialStdMassFlux;
result.passExactInitial = S.initialStdClassic == 0;
result.passHybridRan = isfinite(S.meanMassFluxDivResidualMassFlux);
result.passLowKFinite = isfinite(S.meanLowKMassFlux);
result.passed = result.passFinite && result.passMaps && result.passSameInitial && ...
    result.passExactInitial && result.passHybridRan && result.passLowKFinite;

fprintf('\n=== test_density_homogeneity_lowk_mass_flux_compare_short ===\n');
fprintf('initial std C/P/H      : %.6g / %.6g / %.6g\n', ...
    result.initialStdClassic, result.initialStdProjected, result.initialStdHybrid);
fprintf('mean std C/P/H         : %.6g / %.6g / %.6g\n', ...
    result.meanStdClassic, result.meanStdProjected, result.meanStdHybrid);
fprintf('low-k C/P/H            : %.6g / %.6g / %.6g\n', ...
    result.lowKClassic, result.lowKProjected, result.lowKHybrid);
fprintf('mass-flux residual H   : %.6g\n', result.massFluxResidual);
fprintf('passed                 : %d\n', result.passed);

if ~result.passed
    error('test_density_homogeneity_lowk_mass_flux_compare_short failed.');
end
end
