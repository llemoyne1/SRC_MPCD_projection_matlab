function result = test_density_homogeneity_repair_compare_short()
%TEST_DENSITY_HOMOGENEITY_REPAIR_COMPARE_SHORT Short Q7c integration smoke test.

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
params.seed = 123;
params.initialPopulationMode = 'exact_per_cell';
params.projectedStrength = 1.0;
params.virialDensityRepairStrength = 0.02;
params.virialDensityRepairK = 1.0;
params.virialDensityRepairSmoothPasses = 1;
params.virialDensityRepairMaxDisplacementFraction = 0.05;
params.virialDensityRepairRestoreVelocity = true;
params.makeFigures = false;
params.saveFigures = false;
params.writeSummary = false;

cmp = run_compare_density_projection_repair_poiseuille(params);

result = struct();
result.initialStdClassic = cmp.summary.initialStdClassic;
result.initialStdProjected = cmp.summary.initialStdProjected;
result.initialStdRepair = cmp.summary.initialStdRepair;
result.meanStdClassic = cmp.summary.meanStdClassic;
result.meanStdProjected = cmp.summary.meanStdProjected;
result.meanStdRepair = cmp.summary.meanStdRepair;
result.meanLowKClassic = cmp.summary.meanLowKClassic;
result.meanLowKProjected = cmp.summary.meanLowKProjected;
result.meanLowKRepair = cmp.summary.meanLowKRepair;
result.meanRhoTransportClassic = cmp.summary.meanRhoTransportClassic;
result.meanRhoTransportProjected = cmp.summary.meanRhoTransportProjected;
result.meanRhoTransportRepair = cmp.summary.meanRhoTransportRepair;
result.meanRepairDisplacementRms = cmp.summary.meanRepairDisplacementRms;
result.meanRmsDivAfterRepair = cmp.summary.meanRmsDivAfterRepair;

result.passFinite = all(isfinite([result.meanStdClassic, result.meanStdProjected, result.meanStdRepair, ...
    result.meanLowKClassic, result.meanLowKProjected, result.meanLowKRepair, ...
    result.meanRhoTransportClassic, result.meanRhoTransportProjected, result.meanRhoTransportRepair]));
result.passMaps = ~isempty(cmp.classic.NMaps) && ~isempty(cmp.projected.NMaps) && ~isempty(cmp.repair.NMaps);
result.passSameInitial = result.initialStdClassic == result.initialStdProjected && result.initialStdProjected == result.initialStdRepair;
result.passExactInitial = result.initialStdClassic == 0;
result.passRepairRan = result.meanRepairDisplacementRms > 0;
result.passTransportFinite = isfinite(result.meanRhoTransportRepair);
result.passed = result.passFinite && result.passMaps && result.passSameInitial && ...
    result.passExactInitial && result.passRepairRan && result.passTransportFinite;

fprintf('\n=== test_density_homogeneity_repair_compare_short ===\n');
fprintf('initial std C/P/R      : %.6g / %.6g / %.6g\n', result.initialStdClassic, result.initialStdProjected, result.initialStdRepair);
fprintf('mean std C/P/R         : %.6g / %.6g / %.6g\n', result.meanStdClassic, result.meanStdProjected, result.meanStdRepair);
fprintf('mean lowK C/P/R        : %.6g / %.6g / %.6g\n', result.meanLowKClassic, result.meanLowKProjected, result.meanLowKRepair);
fprintf('rho transport C/P/R    : %.6g / %.6g / %.6g\n', result.meanRhoTransportClassic, result.meanRhoTransportProjected, result.meanRhoTransportRepair);
fprintf('repair dx/div          : %.6g / %.6g\n', result.meanRepairDisplacementRms, result.meanRmsDivAfterRepair);
fprintf('passFinite             : %d\n', result.passFinite);
fprintf('passMaps               : %d\n', result.passMaps);
fprintf('passSameInitial        : %d\n', result.passSameInitial);
fprintf('passExactInitial       : %d\n', result.passExactInitial);
fprintf('passRepairRan          : %d\n', result.passRepairRan);
fprintf('passTransportFinite    : %d\n', result.passTransportFinite);
fprintf('passed                 : %d\n', result.passed);

if ~result.passed
    error('test_density_homogeneity_repair_compare_short failed.');
end
end
