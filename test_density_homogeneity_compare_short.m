function result = test_density_homogeneity_compare_short()
%TEST_DENSITY_HOMOGENEITY_COMPARE_SHORT Smoke test for Q5 density comparison.

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

cmp = run_compare_density_projection_poiseuille(params);
S = cmp.summary;

passFinite = all(isfinite([S.meanStdClassic, S.meanStdProjected, S.meanRmsRelClassic, ...
    S.meanRmsRelProjected, S.meanLowKClassic, S.meanLowKProjected, ...
    S.meanRhoTransportClassic, S.meanRhoTransportProjected]));
passMaps = ~isempty(cmp.classic.NMaps) && ~isempty(cmp.projected.NMaps) && ...
    isequal(size(cmp.classic.NMaps), size(cmp.projected.NMaps));
passSameInitial = isequal(cmp.classic.NMaps(:, :, 1), cmp.projected.NMaps(:, :, 1));
passExactInitial = cmp.metricsClassic.summary.initialStdN == 0 && ...
    cmp.metricsProjected.summary.initialStdN == 0 && ...
    all(cmp.classic.NMaps(:, :, 1) == params.gamma, 'all') && ...
    all(cmp.projected.NMaps(:, :, 1) == params.gamma, 'all');
passProjectionRan = cmp.projected.params.projectionStrength > cmp.classic.params.projectionStrength;
passNoDirectPopulationChange = cmp.projected.summary.maxPopDeltaProjection == 0;
passTransportFinite = isfinite(S.ratioMeanRhoTransport);
passed = passFinite && passMaps && passSameInitial && passExactInitial && passProjectionRan && ...
    passNoDirectPopulationChange && passTransportFinite;

result = struct();
result.cmp = cmp;
result.passFinite = passFinite;
result.passMaps = passMaps;
result.passSameInitial = passSameInitial;
result.passExactInitial = passExactInitial;
result.passProjectionRan = passProjectionRan;
result.passNoDirectPopulationChange = passNoDirectPopulationChange;
result.passTransportFinite = passTransportFinite;
result.passed = passed;

fprintf('\n=== test_density_homogeneity_compare_short ===\n');
fprintf('initial std classic/projected: %.6g / %.6g\n', S.initialStdClassic, S.initialStdProjected);
fprintf('mean std classic/projected  : %.6g / %.6g\n', S.meanStdClassic, S.meanStdProjected);
fprintf('mean rmsRel classic/projected: %.6g / %.6g\n', S.meanRmsRelClassic, S.meanRmsRelProjected);
fprintf('mean lowK classic/projected : %.6g / %.6g\n', S.meanLowKClassic, S.meanLowKProjected);
fprintf('mean rho transport C/P      : %.6g / %.6g\n', S.meanRhoTransportClassic, S.meanRhoTransportProjected);
fprintf('passFinite                  : %d\n', passFinite);
fprintf('passMaps                    : %d\n', passMaps);
fprintf('passSameInitial             : %d\n', passSameInitial);
fprintf('passExactInitial            : %d\n', passExactInitial);
fprintf('passProjectionRan           : %d\n', passProjectionRan);
fprintf('passNoDirectPopulationChange: %d\n', passNoDirectPopulationChange);
fprintf('passTransportFinite         : %d\n', passTransportFinite);
fprintf('passed                      : %d\n', passed);

if ~passed
    error('test_density_homogeneity_compare_short failed.');
end
end
