function result = test_density_homogeneity_virial_compare_short()
%TEST_DENSITY_HOMOGENEITY_VIRIAL_COMPARE_SHORT Smoke test for Q7 comparison.

params = struct();
params.Nx = 12;
params.Ny = 8;
params.Lx = 1.5;
params.Ly = 1.0;
params.gamma = 6;
params.nSteps = 60;
params.sampleEvery = 20;
params.progressEvery = 0;
params.bodyForceX = 0.01;
params.wallModeY = 'thermalize';
params.thermostatAfterProjection = true;
params.projectionTransportDiagnosticsEnable = true;
params.densityTransportDiagnosticsEnable = true;
params.computeFullDiagnosticsEveryStep = false;
params.seed = 123;
params.makeFigures = false;
params.writeSummary = false;
params.projectedStrength = 1.0;
params.initialPopulationMode = 'exact_per_cell';
params.virialDensityKickStrength = 0.05;
params.virialK = 1.0;
params.virialSmoothPasses = 0;

cmp = run_compare_density_projection_virial_poiseuille(params);
S = cmp.summary;

passFinite = all(isfinite([S.meanStdClassic, S.meanStdProjected, S.meanStdVirial, ...
    S.meanLowKClassic, S.meanLowKProjected, S.meanLowKVirial, ...
    S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.meanRhoTransportVirial, ...
    S.meanVirialDVRms, S.meanRmsDivAfterVirialRaw]));
passMaps = ~isempty(cmp.classic.NMaps) && ~isempty(cmp.projected.NMaps) && ~isempty(cmp.virial.NMaps) && ...
    isequal(size(cmp.classic.NMaps), size(cmp.projected.NMaps)) && ...
    isequal(size(cmp.classic.NMaps), size(cmp.virial.NMaps));
passSameInitial = isequal(cmp.classic.NMaps(:, :, 1), cmp.projected.NMaps(:, :, 1)) && ...
    isequal(cmp.classic.NMaps(:, :, 1), cmp.virial.NMaps(:, :, 1));
passExactInitial = cmp.metricsClassic.summary.initialStdN == 0 && ...
    cmp.metricsProjected.summary.initialStdN == 0 && cmp.metricsVirial.summary.initialStdN == 0;
passVirialRan = cmp.virial.params.useVirialDensityKick && cmp.virial.summary.meanVirialDVRms > 0;
passNoDirectPopulationChange = cmp.virial.summary.maxPopDeltaProjection == 0;
passTransportFinite = isfinite(S.ratioMeanRhoTransportVirialProjected);
passed = passFinite && passMaps && passSameInitial && passExactInitial && ...
    passVirialRan && passNoDirectPopulationChange && passTransportFinite;

result = struct();
result.cmp = cmp;
result.passFinite = passFinite;
result.passMaps = passMaps;
result.passSameInitial = passSameInitial;
result.passExactInitial = passExactInitial;
result.passVirialRan = passVirialRan;
result.passNoDirectPopulationChange = passNoDirectPopulationChange;
result.passTransportFinite = passTransportFinite;
result.passed = passed;

fprintf('\n=== test_density_homogeneity_virial_compare_short ===\n');
fprintf('initial std C/P/V      : %.6g / %.6g / %.6g\n', S.initialStdClassic, S.initialStdProjected, S.initialStdVirial);
fprintf('mean std C/P/V         : %.6g / %.6g / %.6g\n', S.meanStdClassic, S.meanStdProjected, S.meanStdVirial);
fprintf('mean lowK C/P/V        : %.6g / %.6g / %.6g\n', S.meanLowKClassic, S.meanLowKProjected, S.meanLowKVirial);
fprintf('rho transport C/P/V    : %.6g / %.6g / %.6g\n', S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.meanRhoTransportVirial);
fprintf('virial dv rms          : %.6g\n', S.meanVirialDVRms);
fprintf('passFinite             : %d\n', passFinite);
fprintf('passMaps               : %d\n', passMaps);
fprintf('passSameInitial        : %d\n', passSameInitial);
fprintf('passExactInitial       : %d\n', passExactInitial);
fprintf('passVirialRan          : %d\n', passVirialRan);
fprintf('passNoDirectPopulation : %d\n', passNoDirectPopulationChange);
fprintf('passTransportFinite    : %d\n', passTransportFinite);
fprintf('passed                 : %d\n', passed);

if ~passed
    error('test_density_homogeneity_virial_compare_short failed.');
end
end
