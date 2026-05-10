function cmp = run_compare_density_projection_poiseuille(params)
%RUN_COMPARE_DENSITY_PROJECTION_POISEUILLE Compare density homogeneity with/without projection.
%
%   cmp = run_compare_density_projection_poiseuille(params)
%
% Runs two paired Poiseuille simulations with the same initial random seed:
%   classic   : projectionStrength = 0
%   projected : projectionStrength = 1, unless params.projectedStrength is set
%
% The comparison focuses on cell-wise density/population evolution, not on
% perfect Poiseuille convergence. It stores full cell population maps at sample
% times and computes temporal metrics, density maps and low-k density energy.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);

base = params;
base.storeDensityMaps = true;
base.makeFigures = false;
base.computeFullDiagnosticsEveryStep = false;

classicParams = base;
classicParams.projectionEnable = true;
classicParams.projectionStrength = 0.0;
classicParams.seed = params.seed;

projectedParams = base;
projectedParams.projectionEnable = true;
projectedParams.projectionStrength = params.projectedStrength;
projectedParams.seed = params.seed;

fprintf('\n=== density homogeneity comparison: classic projectionStrength=0 ===\n');
outClassic = run_projection_poiseuille_long_demo(classicParams);

fprintf('\n=== density homogeneity comparison: projected projectionStrength=%.6g ===\n', projectedParams.projectionStrength);
outProjected = run_projection_poiseuille_long_demo(projectedParams);

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

summary = compare_summaries(metricsClassic, metricsProjected, outClassic, outProjected);

cmp = struct();
cmp.params = params;
cmp.classic = outClassic;
cmp.projected = outProjected;
cmp.metricsClassic = metricsClassic;
cmp.metricsProjected = metricsProjected;
cmp.summary = summary;

fprintf('\n=== density homogeneity summary ===\n');
fprintf('seed                                  : %d\n', params.seed);
fprintf('steps/sampleEvery                    : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf('grid, gamma                          : %d x %d, %.6g\n', params.Nx, params.Ny, params.gamma);
fprintf('initialPopulationMode                : %s\n', char(params.initialPopulationMode));
fprintf('initial std(N) classic/projected     : %.6g / %.6g\n', ...
    summary.initialStdClassic, summary.initialStdProjected);
fprintf('mean std(N) classic/projected        : %.6g / %.6g  ratio=%.6g\n', ...
    summary.meanStdClassic, summary.meanStdProjected, summary.ratioMeanStd);
fprintf('final std(N) classic/projected       : %.6g / %.6g  ratio=%.6g\n', ...
    summary.finalStdClassic, summary.finalStdProjected, summary.ratioFinalStd);
fprintf('mean rmsRel classic/projected        : %.6g / %.6g  ratio=%.6g\n', ...
    summary.meanRmsRelClassic, summary.meanRmsRelProjected, summary.ratioMeanRmsRel);
fprintf('mean out-band classic/projected      : %.6g / %.6g  ratio=%.6g\n', ...
    summary.meanOutBandClassic, summary.meanOutBandProjected, summary.ratioMeanOutBand);
fprintf('mean low-k energy classic/projected  : %.6g / %.6g  ratio=%.6g\n', ...
    summary.meanLowKClassic, summary.meanLowKProjected, summary.ratioMeanLowK);
fprintf('mean rho transport classic/projected : %.6g / %.6g  ratio=%.6g\n', ...
    summary.meanRhoTransportClassic, summary.meanRhoTransportProjected, summary.ratioMeanRhoTransport);
fprintf('time-avg rel RMS classic/projected   : %.6g / %.6g  ratio=%.6g\n', ...
    summary.timeAvgRelRmsClassic, summary.timeAvgRelRmsProjected, summary.ratioTimeAvgRelRms);
fprintf('elapsed classic/projected [s]        : %.2f / %.2f\n', outClassic.elapsedWallClock, outProjected.elapsedWallClock);

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
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'saveFigures', false);
params = set_default(params, 'outputDir', 'density_projection_output');
params = set_default(params, 'writeSummary', false);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function summary = compare_summaries(mC, mP, outC, outP)
summary = struct();
summary.initialStdClassic = mC.summary.initialStdN;
summary.initialStdProjected = mP.summary.initialStdN;
summary.initialRmsRelClassic = mC.summary.initialRmsRel;
summary.initialRmsRelProjected = mP.summary.initialRmsRel;
summary.meanStdClassic = mC.summary.meanStdN;
summary.meanStdProjected = mP.summary.meanStdN;
summary.ratioMeanStd = safe_ratio(summary.meanStdProjected, summary.meanStdClassic);
summary.finalStdClassic = mC.summary.finalStdN;
summary.finalStdProjected = mP.summary.finalStdN;
summary.ratioFinalStd = safe_ratio(summary.finalStdProjected, summary.finalStdClassic);
summary.meanRmsRelClassic = mC.summary.meanRmsRel;
summary.meanRmsRelProjected = mP.summary.meanRmsRel;
summary.ratioMeanRmsRel = safe_ratio(summary.meanRmsRelProjected, summary.meanRmsRelClassic);
summary.meanOutBandClassic = mC.summary.meanOutBandFraction;
summary.meanOutBandProjected = mP.summary.meanOutBandFraction;
summary.ratioMeanOutBand = safe_ratio(summary.meanOutBandProjected, summary.meanOutBandClassic);
summary.meanLowKClassic = mC.summary.meanLowKEnergy;
summary.meanLowKProjected = mP.summary.meanLowKEnergy;
summary.ratioMeanLowK = safe_ratio(summary.meanLowKProjected, summary.meanLowKClassic);
summary.timeAvgRelRmsClassic = mC.summary.timeAvgRelRms;
summary.timeAvgRelRmsProjected = mP.summary.timeAvgRelRms;
summary.ratioTimeAvgRelRms = safe_ratio(summary.timeAvgRelRmsProjected, summary.timeAvgRelRmsClassic);
summary.meanRhoTransportClassic = outC.summary.meanDensityTransportProjectedRms;
summary.meanRhoTransportProjected = outP.summary.meanDensityTransportProjectedRms;
summary.ratioMeanRhoTransport = safe_ratio(summary.meanRhoTransportProjected, summary.meanRhoTransportClassic);
summary.nuEffClassic = outC.viscosity.nuEff;
summary.nuEffProjected = outP.viscosity.nuEff;
summary.R2Classic = outC.viscosity.R2;
summary.R2Projected = outP.viscosity.R2;
summary.signalToNoiseClassic = outC.viscosity.signalToNoise;
summary.signalToNoiseProjected = outP.viscosity.signalToNoise;
end

function r = safe_ratio(a, b)
r = a / max(b, eps);
end

function write_summary_file(cmp)
outDir = cmp.params.outputDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fname = fullfile(outDir, 'density_projection_summary.txt');
fid = fopen(fname, 'w');
if fid < 0
    warning('Could not write summary file: %s', fname);
    return;
end
cleaner = onCleanup(@() fclose(fid));
S = cmp.summary;
fprintf(fid, 'Density projection comparison summary\n');
fprintf(fid, '=====================================\n\n');
fprintf(fid, 'seed = %d\n', cmp.params.seed);
fprintf(fid, 'Nx = %d\nNy = %d\ngamma = %.16g\nnSteps = %d\nsampleEvery = %d\n', ...
    cmp.params.Nx, cmp.params.Ny, cmp.params.gamma, cmp.params.nSteps, cmp.params.sampleEvery);
fprintf(fid, 'initialPopulationMode = %s\n', char(cmp.params.initialPopulationMode));
fprintf(fid, 'initialStdClassic = %.16g\ninitialStdProjected = %.16g\n', S.initialStdClassic, S.initialStdProjected);
fprintf(fid, 'meanStdClassic = %.16g\nmeanStdProjected = %.16g\nratioMeanStd = %.16g\n', S.meanStdClassic, S.meanStdProjected, S.ratioMeanStd);
fprintf(fid, 'meanRmsRelClassic = %.16g\nmeanRmsRelProjected = %.16g\nratioMeanRmsRel = %.16g\n', S.meanRmsRelClassic, S.meanRmsRelProjected, S.ratioMeanRmsRel);
fprintf(fid, 'meanLowKClassic = %.16g\nmeanLowKProjected = %.16g\nratioMeanLowK = %.16g\n', S.meanLowKClassic, S.meanLowKProjected, S.ratioMeanLowK);
fprintf(fid, 'meanRhoTransportClassic = %.16g\nmeanRhoTransportProjected = %.16g\nratioMeanRhoTransport = %.16g\n', S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.ratioMeanRhoTransport);
fprintf('summary file written: %s\n', fname);
end
