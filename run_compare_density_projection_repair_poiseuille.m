function cmp = run_compare_density_projection_repair_poiseuille(params)
%RUN_COMPARE_DENSITY_PROJECTION_REPAIR_POISEUILLE Q7c density-repair benchmark.
%
%   cmp = run_compare_density_projection_repair_poiseuille(params)
%
% Runs the same Poiseuille case three times with identical initial seed:
%   A. SRC/MPCD without applied projection             projectionStrength = 0
%   B. pressure projection only                       projectionStrength = projectedStrength
%   C. pressure projection + density position repair  projectionStrength = projectedStrength
%
% The density repair moves particles slightly along -grad(Pvir), then restores
% the pressure-projected cell-mean velocity field. It is designed to repair rho
% without directly using a virial velocity kick as the hydrodynamic flow.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);

base = params;
base.storeDensityMaps = true;
base.makeFigures = false;
base.computeFullDiagnosticsEveryStep = false;
base.useVirialDensityKick = false;
base.virialDensityKickStrength = 0.0;

classicParams = base;
classicParams.projectionEnable = true;
classicParams.projectionStrength = 0.0;
classicParams.useVirialDensityRepair = false;
classicParams.virialDensityRepairStrength = 0.0;
classicParams.seed = params.seed;

projectedParams = base;
projectedParams.projectionEnable = true;
projectedParams.projectionStrength = params.projectedStrength;
projectedParams.useVirialDensityRepair = false;
projectedParams.virialDensityRepairStrength = 0.0;
projectedParams.seed = params.seed;

repairParams = base;
repairParams.projectionEnable = true;
repairParams.projectionStrength = params.projectedStrength;
repairParams.useVirialDensityRepair = true;
repairParams.virialDensityRepairStrength = params.virialDensityRepairStrength;
repairParams.virialDensityRepairK = params.virialDensityRepairK;
repairParams.virialDensityRepairSmoothPasses = params.virialDensityRepairSmoothPasses;
repairParams.virialDensityRepairMinCellCount = params.virialDensityRepairMinCellCount;
repairParams.virialDensityRepairMaxDisplacementFraction = params.virialDensityRepairMaxDisplacementFraction;
repairParams.virialDensityRepairRestoreVelocity = params.virialDensityRepairRestoreVelocity;
repairParams.virialDensityRepairInterpolationMethod = params.virialDensityRepairInterpolationMethod;
repairParams.seed = params.seed;

fprintf('\n=== Q7c density comparison: classic projectionStrength=0 ===\n');
outClassic = run_projection_poiseuille_long_demo(classicParams);

fprintf('\n=== Q7c density comparison: pressure projection only strength=%.6g ===\n', projectedParams.projectionStrength);
outProjected = run_projection_poiseuille_long_demo(projectedParams);

fprintf('\n=== Q7c density comparison: projection + position repair strength=%.6g, K=%.6g ===\n', ...
    repairParams.virialDensityRepairStrength, repairParams.virialDensityRepairK);
outRepair = run_projection_poiseuille_long_demo(repairParams);

metricsClassic = projection_density_homogeneity_metrics(outClassic.NMaps, classicParams, ...
    'sampleTimes', outClassic.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex);
metricsProjected = projection_density_homogeneity_metrics(outProjected.NMaps, projectedParams, ...
    'sampleTimes', outProjected.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex);
metricsRepair = projection_density_homogeneity_metrics(outRepair.NMaps, repairParams, ...
    'sampleTimes', outRepair.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex);

summary = compare_three_summaries(metricsClassic, metricsProjected, metricsRepair, outClassic, outProjected, outRepair);

cmp = struct();
cmp.params = params;
cmp.classic = outClassic;
cmp.projected = outProjected;
cmp.repair = outRepair;
cmp.metricsClassic = metricsClassic;
cmp.metricsProjected = metricsProjected;
cmp.metricsRepair = metricsRepair;
cmp.summary = summary;

fprintf('\n=== Q7c density / position-repair homogeneity summary ===\n');
fprintf('seed                                      : %d\n', params.seed);
fprintf('steps/sampleEvery                        : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf('grid, gamma                              : %d x %d, %.6g\n', params.Nx, params.Ny, params.gamma);
fprintf('initialPopulationMode                    : %s\n', char(params.initialPopulationMode));
fprintf('repair strength, K, smooth, maxFrac      : %.6g / %.6g / %d / %.6g\n', ...
    params.virialDensityRepairStrength, params.virialDensityRepairK, ...
    params.virialDensityRepairSmoothPasses, params.virialDensityRepairMaxDisplacementFraction);
fprintf('initial std C/P/R                        : %.6g / %.6g / %.6g\n', ...
    summary.initialStdClassic, summary.initialStdProjected, summary.initialStdRepair);
fprintf('mean std(N) C/P/R                        : %.6g / %.6g / %.6g\n', ...
    summary.meanStdClassic, summary.meanStdProjected, summary.meanStdRepair);
fprintf('mean std ratios P/C, R/C, R/P            : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanStdProjectedClassic, summary.ratioMeanStdRepairClassic, summary.ratioMeanStdRepairProjected);
fprintf('mean out-band C/P/R                      : %.6g / %.6g / %.6g\n', ...
    summary.meanOutBandClassic, summary.meanOutBandProjected, summary.meanOutBandRepair);
fprintf('mean low-k C/P/R                         : %.6g / %.6g / %.6g\n', ...
    summary.meanLowKClassic, summary.meanLowKProjected, summary.meanLowKRepair);
fprintf('low-k ratios P/C, R/C, R/P               : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanLowKProjectedClassic, summary.ratioMeanLowKRepairClassic, summary.ratioMeanLowKRepairProjected);
fprintf('mean rho transport C/P/R                 : %.6g / %.6g / %.6g\n', ...
    summary.meanRhoTransportClassic, summary.meanRhoTransportProjected, summary.meanRhoTransportRepair);
fprintf('rho transport ratios P/C, R/C, R/P       : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanRhoTransportProjectedClassic, summary.ratioMeanRhoTransportRepairClassic, summary.ratioMeanRhoTransportRepairProjected);
fprintf('time-avg rel RMS C/P/R                   : %.6g / %.6g / %.6g\n', ...
    summary.timeAvgRelRmsClassic, summary.timeAvgRelRmsProjected, summary.timeAvgRelRmsRepair);
fprintf('nu_eff C/P/R                             : %.6g / %.6g / %.6g\n', ...
    summary.nuEffClassic, summary.nuEffProjected, summary.nuEffRepair);
fprintf('R2 C/P/R                                 : %.6g / %.6g / %.6g\n', ...
    summary.R2Classic, summary.R2Projected, summary.R2Repair);
fprintf('repair mean dx rms, restore residual     : %.6g / %.6g\n', ...
    summary.meanRepairDisplacementRms, summary.meanRepairVelocityRestoreResidualAfterRms);
fprintf('mean div after repair                    : %.6g\n', summary.meanRmsDivAfterRepair);
fprintf('elapsed C/P/R [s]                        : %.2f / %.2f / %.2f\n', ...
    outClassic.elapsedWallClock, outProjected.elapsedWallClock, outRepair.elapsedWallClock);

if params.makeFigures
    plot_density_homogeneity_comparison(cmp, 'saveFigures', params.saveFigures, 'outputDir', params.outputDir);
end

if params.writeSummary
    write_summary_file(cmp);
end
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 16);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 1.0);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 2.0e-2);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'nSteps', 10000);
params = set_default(params, 'sampleEvery', 100);
params = set_default(params, 'progressEvery', 1000);
params = set_default(params, 'maxWallClockSeconds', Inf);
params = set_default(params, 'projectedStrength', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'projectionTransportDiagnosticsEnable', true);
params = set_default(params, 'densityTransportDiagnosticsEnable', true);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'thermostatMinParticlesPerCell', 2);
params = set_default(params, 'thermostatMaxScale', 5.0);
params = set_default(params, 'wallModeY', 'thermalize');
params = set_default(params, 'wallSigma', sqrt(params.kBT));
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'seed', 11);
params = set_default(params, 'excludeWallCellsFit', 2);
params = set_default(params, 'fitStartFraction', 0.5);
params = set_default(params, 'densityBandFraction', 0.20);
params = set_default(params, 'densityExcludeWallCells', 1);
params = set_default(params, 'lowKMaxIndex', 2);
params = set_default(params, 'virialDensityRepairStrength', 0.02);
params = set_default(params, 'virialDensityRepairK', 1.0);
params = set_default(params, 'virialDensityRepairSmoothPasses', 1);
params = set_default(params, 'virialDensityRepairMinCellCount', 1.0);
params = set_default(params, 'virialDensityRepairMaxDisplacementFraction', 0.05);
params = set_default(params, 'virialDensityRepairRestoreVelocity', true);
params = set_default(params, 'virialDensityRepairInterpolationMethod', 'nearest');
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'saveFigures', false);
params = set_default(params, 'outputDir', 'density_projection_repair_output');
params = set_default(params, 'writeSummary', false);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function summary = compare_three_summaries(mC, mP, mR, outC, outP, outR)
summary = struct();
summary.initialStdClassic = mC.summary.initialStdN;
summary.initialStdProjected = mP.summary.initialStdN;
summary.initialStdRepair = mR.summary.initialStdN;
summary.meanStdClassic = mC.summary.meanStdN;
summary.meanStdProjected = mP.summary.meanStdN;
summary.meanStdRepair = mR.summary.meanStdN;
summary.ratioMeanStdProjectedClassic = safe_ratio(summary.meanStdProjected, summary.meanStdClassic);
summary.ratioMeanStdRepairClassic = safe_ratio(summary.meanStdRepair, summary.meanStdClassic);
summary.ratioMeanStdRepairProjected = safe_ratio(summary.meanStdRepair, summary.meanStdProjected);
summary.finalStdClassic = mC.summary.finalStdN;
summary.finalStdProjected = mP.summary.finalStdN;
summary.finalStdRepair = mR.summary.finalStdN;
summary.meanRmsRelClassic = mC.summary.meanRmsRel;
summary.meanRmsRelProjected = mP.summary.meanRmsRel;
summary.meanRmsRelRepair = mR.summary.meanRmsRel;
summary.meanOutBandClassic = mC.summary.meanOutBandFraction;
summary.meanOutBandProjected = mP.summary.meanOutBandFraction;
summary.meanOutBandRepair = mR.summary.meanOutBandFraction;
summary.meanLowKClassic = mC.summary.meanLowKEnergy;
summary.meanLowKProjected = mP.summary.meanLowKEnergy;
summary.meanLowKRepair = mR.summary.meanLowKEnergy;
summary.ratioMeanLowKProjectedClassic = safe_ratio(summary.meanLowKProjected, summary.meanLowKClassic);
summary.ratioMeanLowKRepairClassic = safe_ratio(summary.meanLowKRepair, summary.meanLowKClassic);
summary.ratioMeanLowKRepairProjected = safe_ratio(summary.meanLowKRepair, summary.meanLowKProjected);
summary.timeAvgRelRmsClassic = mC.summary.timeAvgRelRms;
summary.timeAvgRelRmsProjected = mP.summary.timeAvgRelRms;
summary.timeAvgRelRmsRepair = mR.summary.timeAvgRelRms;
summary.meanRhoTransportClassic = outC.summary.meanDensityTransportProjectedRms;
summary.meanRhoTransportProjected = outP.summary.meanDensityTransportProjectedRms;
summary.meanRhoTransportRepair = outR.summary.meanDensityTransportProjectedRms;
summary.ratioMeanRhoTransportProjectedClassic = safe_ratio(summary.meanRhoTransportProjected, summary.meanRhoTransportClassic);
summary.ratioMeanRhoTransportRepairClassic = safe_ratio(summary.meanRhoTransportRepair, summary.meanRhoTransportClassic);
summary.ratioMeanRhoTransportRepairProjected = safe_ratio(summary.meanRhoTransportRepair, summary.meanRhoTransportProjected);
summary.nuEffClassic = outC.viscosity.nuEff;
summary.nuEffProjected = outP.viscosity.nuEff;
summary.nuEffRepair = outR.viscosity.nuEff;
summary.R2Classic = outC.viscosity.R2;
summary.R2Projected = outP.viscosity.R2;
summary.R2Repair = outR.viscosity.R2;
summary.signalToNoiseClassic = outC.viscosity.signalToNoise;
summary.signalToNoiseProjected = outP.viscosity.signalToNoise;
summary.signalToNoiseRepair = outR.viscosity.signalToNoise;
summary.meanRepairDisplacementRms = outR.summary.meanDensityRepairDisplacementRms;
summary.meanRepairStdBefore = outR.summary.meanDensityRepairStdBefore;
summary.meanRepairStdAfter = outR.summary.meanDensityRepairStdAfter;
summary.meanRepairVelocityRestoreResidualAfterRms = outR.summary.meanDensityRepairVelocityRestoreResidualAfterRms;
summary.meanRmsDivAfterRepair = outR.summary.meanRmsDivAfterDensityRepair;
summary.lastRepairDisplacementRms = outR.summary.lastDensityRepairDisplacementRms;
summary.lastRmsDivAfterRepair = outR.summary.lastRmsDivAfterDensityRepair;
end

function r = safe_ratio(a, b)
r = a / max(b, eps);
end

function write_summary_file(cmp)
outDir = cmp.params.outputDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fname = fullfile(outDir, 'density_projection_repair_summary.txt');
fid = fopen(fname, 'w');
if fid < 0
    warning('Could not write summary file: %s', fname);
    return;
end
cleaner = onCleanup(@() fclose(fid));
S = cmp.summary;
fprintf(fid, 'Q7c density projection + position repair comparison summary\n');
fprintf(fid, '=========================================================\n\n');
fprintf(fid, 'seed = %d\n', cmp.params.seed);
fprintf(fid, 'Nx = %d\nNy = %d\ngamma = %.16g\nnSteps = %d\nsampleEvery = %d\n', ...
    cmp.params.Nx, cmp.params.Ny, cmp.params.gamma, cmp.params.nSteps, cmp.params.sampleEvery);
fprintf(fid, 'initialPopulationMode = %s\n', char(cmp.params.initialPopulationMode));
fprintf(fid, 'virialDensityRepairStrength = %.16g\nvirialDensityRepairK = %.16g\nvirialDensityRepairSmoothPasses = %d\n', ...
    cmp.params.virialDensityRepairStrength, cmp.params.virialDensityRepairK, cmp.params.virialDensityRepairSmoothPasses);
fprintf(fid, 'meanStdClassic = %.16g\nmeanStdProjected = %.16g\nmeanStdRepair = %.16g\n', ...
    S.meanStdClassic, S.meanStdProjected, S.meanStdRepair);
fprintf(fid, 'ratioMeanStdProjectedClassic = %.16g\nratioMeanStdRepairClassic = %.16g\nratioMeanStdRepairProjected = %.16g\n', ...
    S.ratioMeanStdProjectedClassic, S.ratioMeanStdRepairClassic, S.ratioMeanStdRepairProjected);
fprintf(fid, 'meanLowKClassic = %.16g\nmeanLowKProjected = %.16g\nmeanLowKRepair = %.16g\n', ...
    S.meanLowKClassic, S.meanLowKProjected, S.meanLowKRepair);
fprintf(fid, 'meanRhoTransportClassic = %.16g\nmeanRhoTransportProjected = %.16g\nmeanRhoTransportRepair = %.16g\n', ...
    S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.meanRhoTransportRepair);
fprintf(fid, 'nuEffClassic = %.16g\nnuEffProjected = %.16g\nnuEffRepair = %.16g\n', ...
    S.nuEffClassic, S.nuEffProjected, S.nuEffRepair);
fprintf(fid, 'R2Classic = %.16g\nR2Projected = %.16g\nR2Repair = %.16g\n', ...
    S.R2Classic, S.R2Projected, S.R2Repair);
fprintf(fid, 'meanRepairDisplacementRms = %.16g\nmeanRmsDivAfterRepair = %.16g\n', ...
    S.meanRepairDisplacementRms, S.meanRmsDivAfterRepair);
fprintf('summary file written: %s\n', fname);
end
