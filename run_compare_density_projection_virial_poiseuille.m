function cmp = run_compare_density_projection_virial_poiseuille(params)
%RUN_COMPARE_DENSITY_PROJECTION_VIRIAL_POISEUILLE Q7 three-way density benchmark.
%
%   cmp = run_compare_density_projection_virial_poiseuille(params)
%
% Runs the same Poiseuille case three times, with identical initial seed:
%   A. SRC/MPCD without applied projection       projectionStrength = 0
%   B. pressure projection only                 projectionStrength = projectedStrength
%   C. pressure projection + weak virial kick   projectionStrength = projectedStrength
%
% The virial kick is a local post-projection velocity correction based on
% Pvir = virialK * (N - gamma). It does not move particles directly.

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
classicParams.useVirialDensityKick = false;
classicParams.virialDensityKickStrength = 0.0;
classicParams.seed = params.seed;

projectedParams = base;
projectedParams.projectionEnable = true;
projectedParams.projectionStrength = params.projectedStrength;
projectedParams.useVirialDensityKick = false;
projectedParams.virialDensityKickStrength = 0.0;
projectedParams.seed = params.seed;

virialParams = base;
virialParams.projectionEnable = true;
virialParams.projectionStrength = params.projectedStrength;
virialParams.useVirialDensityKick = true;
virialParams.virialDensityKickStrength = params.virialDensityKickStrength;
virialParams.virialK = params.virialK;
virialParams.virialSmoothPasses = params.virialSmoothPasses;
virialParams.virialMinCellCount = params.virialMinCellCount;
virialParams.virialMaxParticleKick = params.virialMaxParticleKick;
virialParams.virialInterpolationMethod = params.virialInterpolationMethod;
virialParams.seed = params.seed;

fprintf('\n=== Q7 density comparison: classic projectionStrength=0 ===\n');
outClassic = run_projection_poiseuille_long_demo(classicParams);

fprintf('\n=== Q7 density comparison: pressure projection only strength=%.6g ===\n', projectedParams.projectionStrength);
outProjected = run_projection_poiseuille_long_demo(projectedParams);
% 
fprintf('\n=== Q7 density comparison: projection + virial kick strength=%.6g, K=%.6g ===\n', ...
    virialParams.virialDensityKickStrength, virialParams.virialK);
outVirial = run_projection_poiseuille_long_demo(virialParams);

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
metricsVirial = projection_density_homogeneity_metrics(outVirial.NMaps, virialParams, ...
    'sampleTimes', outVirial.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex);

summary = compare_three_summaries(metricsClassic, metricsProjected, metricsVirial, outClassic, outProjected, outVirial);

cmp = struct();
cmp.params = params;
cmp.classic = outClassic;
cmp.projected = outProjected;
cmp.virial = outVirial;
cmp.metricsClassic = metricsClassic;
cmp.metricsProjected = metricsProjected;
cmp.metricsVirial = metricsVirial;
cmp.summary = summary;

fprintf('\n=== Q7 density / virial homogeneity summary ===\n');
fprintf('seed                                      : %d\n', params.seed);
fprintf('steps/sampleEvery                        : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf('grid, gamma                              : %d x %d, %.6g\n', params.Nx, params.Ny, params.gamma);
fprintf('initialPopulationMode                    : %s\n', char(params.initialPopulationMode));
fprintf('virial strength, K, smooth               : %.6g / %.6g / %d\n', ...
    params.virialDensityKickStrength, params.virialK, params.virialSmoothPasses);
fprintf('initial std C/P/V                        : %.6g / %.6g / %.6g\n', ...
    summary.initialStdClassic, summary.initialStdProjected, summary.initialStdVirial);
fprintf('mean std(N) C/P/V                        : %.6g / %.6g / %.6g\n', ...
    summary.meanStdClassic, summary.meanStdProjected, summary.meanStdVirial);
fprintf('mean std ratios P/C, V/C, V/P            : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanStdProjectedClassic, summary.ratioMeanStdVirialClassic, summary.ratioMeanStdVirialProjected);
fprintf('mean out-band C/P/V                      : %.6g / %.6g / %.6g\n', ...
    summary.meanOutBandClassic, summary.meanOutBandProjected, summary.meanOutBandVirial);
fprintf('mean low-k C/P/V                         : %.6g / %.6g / %.6g\n', ...
    summary.meanLowKClassic, summary.meanLowKProjected, summary.meanLowKVirial);
fprintf('low-k ratios P/C, V/C, V/P               : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanLowKProjectedClassic, summary.ratioMeanLowKVirialClassic, summary.ratioMeanLowKVirialProjected);
fprintf('mean rho transport C/P/V                 : %.6g / %.6g / %.6g\n', ...
    summary.meanRhoTransportClassic, summary.meanRhoTransportProjected, summary.meanRhoTransportVirial);
fprintf('rho transport ratios P/C, V/C, V/P       : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanRhoTransportProjectedClassic, summary.ratioMeanRhoTransportVirialClassic, summary.ratioMeanRhoTransportVirialProjected);
fprintf('time-avg rel RMS C/P/V                   : %.6g / %.6g / %.6g\n', ...
    summary.timeAvgRelRmsClassic, summary.timeAvgRelRmsProjected, summary.timeAvgRelRmsVirial);
fprintf('nu_eff C/P/V                             : %.6g / %.6g / %.6g\n', ...
    summary.nuEffClassic, summary.nuEffProjected, summary.nuEffVirial);
fprintf('R2 C/P/V                                 : %.6g / %.6g / %.6g\n', ...
    summary.R2Classic, summary.R2Projected, summary.R2Virial);
fprintf('mean virial dv rms, div after virial     : %.6g / %.6g\n', ...
    summary.meanVirialDVRms, summary.meanRmsDivAfterVirialRaw);
fprintf('elapsed C/P/V [s]                        : %.2f / %.2f / %.2f\n', ...
    outClassic.elapsedWallClock, outProjected.elapsedWallClock, outVirial.elapsedWallClock);

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
params = set_default(params, 'virialDensityKickStrength', 0.05);
params = set_default(params, 'virialK', 1.0);
params = set_default(params, 'virialSmoothPasses', 0);
params = set_default(params, 'virialMinCellCount', 1.0);
params = set_default(params, 'virialMaxParticleKick', Inf);
params = set_default(params, 'virialInterpolationMethod', params.projectionInterpolationMethod);
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'saveFigures', false);
params = set_default(params, 'outputDir', 'density_projection_virial_output');
params = set_default(params, 'writeSummary', false);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function summary = compare_three_summaries(mC, mP, mV, outC, outP, outV)
summary = struct();
summary.initialStdClassic = mC.summary.initialStdN;
summary.initialStdProjected = mP.summary.initialStdN;
summary.initialStdVirial = mV.summary.initialStdN;
summary.meanStdClassic = mC.summary.meanStdN;
summary.meanStdProjected = mP.summary.meanStdN;
summary.meanStdVirial = mV.summary.meanStdN;
summary.ratioMeanStdProjectedClassic = safe_ratio(summary.meanStdProjected, summary.meanStdClassic);
summary.ratioMeanStdVirialClassic = safe_ratio(summary.meanStdVirial, summary.meanStdClassic);
summary.ratioMeanStdVirialProjected = safe_ratio(summary.meanStdVirial, summary.meanStdProjected);
summary.finalStdClassic = mC.summary.finalStdN;
summary.finalStdProjected = mP.summary.finalStdN;
summary.finalStdVirial = mV.summary.finalStdN;
summary.meanRmsRelClassic = mC.summary.meanRmsRel;
summary.meanRmsRelProjected = mP.summary.meanRmsRel;
summary.meanRmsRelVirial = mV.summary.meanRmsRel;
summary.ratioMeanRmsRelProjectedClassic = safe_ratio(summary.meanRmsRelProjected, summary.meanRmsRelClassic);
summary.ratioMeanRmsRelVirialClassic = safe_ratio(summary.meanRmsRelVirial, summary.meanRmsRelClassic);
summary.ratioMeanRmsRelVirialProjected = safe_ratio(summary.meanRmsRelVirial, summary.meanRmsRelProjected);
summary.meanOutBandClassic = mC.summary.meanOutBandFraction;
summary.meanOutBandProjected = mP.summary.meanOutBandFraction;
summary.meanOutBandVirial = mV.summary.meanOutBandFraction;
summary.ratioMeanOutBandProjectedClassic = safe_ratio(summary.meanOutBandProjected, summary.meanOutBandClassic);
summary.ratioMeanOutBandVirialClassic = safe_ratio(summary.meanOutBandVirial, summary.meanOutBandClassic);
summary.ratioMeanOutBandVirialProjected = safe_ratio(summary.meanOutBandVirial, summary.meanOutBandProjected);
summary.meanLowKClassic = mC.summary.meanLowKEnergy;
summary.meanLowKProjected = mP.summary.meanLowKEnergy;
summary.meanLowKVirial = mV.summary.meanLowKEnergy;
summary.ratioMeanLowKProjectedClassic = safe_ratio(summary.meanLowKProjected, summary.meanLowKClassic);
summary.ratioMeanLowKVirialClassic = safe_ratio(summary.meanLowKVirial, summary.meanLowKClassic);
summary.ratioMeanLowKVirialProjected = safe_ratio(summary.meanLowKVirial, summary.meanLowKProjected);
summary.timeAvgRelRmsClassic = mC.summary.timeAvgRelRms;
summary.timeAvgRelRmsProjected = mP.summary.timeAvgRelRms;
summary.timeAvgRelRmsVirial = mV.summary.timeAvgRelRms;
summary.ratioTimeAvgRelRmsProjectedClassic = safe_ratio(summary.timeAvgRelRmsProjected, summary.timeAvgRelRmsClassic);
summary.ratioTimeAvgRelRmsVirialClassic = safe_ratio(summary.timeAvgRelRmsVirial, summary.timeAvgRelRmsClassic);
summary.ratioTimeAvgRelRmsVirialProjected = safe_ratio(summary.timeAvgRelRmsVirial, summary.timeAvgRelRmsProjected);
summary.meanRhoTransportClassic = outC.summary.meanDensityTransportProjectedRms;
summary.meanRhoTransportProjected = outP.summary.meanDensityTransportProjectedRms;
summary.meanRhoTransportVirial = outV.summary.meanDensityTransportProjectedRms;
summary.ratioMeanRhoTransportProjectedClassic = safe_ratio(summary.meanRhoTransportProjected, summary.meanRhoTransportClassic);
summary.ratioMeanRhoTransportVirialClassic = safe_ratio(summary.meanRhoTransportVirial, summary.meanRhoTransportClassic);
summary.ratioMeanRhoTransportVirialProjected = safe_ratio(summary.meanRhoTransportVirial, summary.meanRhoTransportProjected);
summary.meanRhoTransportPressureOnlyVirialRun = outV.summary.meanDensityTransportPressureOnlyProjectedRms;
summary.ratioMeanRhoTransportVirialFinalPressureOnly = safe_ratio(summary.meanRhoTransportVirial, summary.meanRhoTransportPressureOnlyVirialRun);
summary.nuEffClassic = outC.viscosity.nuEff;
summary.nuEffProjected = outP.viscosity.nuEff;
summary.nuEffVirial = outV.viscosity.nuEff;
summary.R2Classic = outC.viscosity.R2;
summary.R2Projected = outP.viscosity.R2;
summary.R2Virial = outV.viscosity.R2;
summary.signalToNoiseClassic = outC.viscosity.signalToNoise;
summary.signalToNoiseProjected = outP.viscosity.signalToNoise;
summary.signalToNoiseVirial = outV.viscosity.signalToNoise;
summary.meanVirialDVRms = outV.summary.meanVirialDVRms;
summary.meanVirialPvirRms = outV.summary.meanVirialPvirRms;
summary.meanRmsDivAfterVirialRaw = outV.summary.meanRmsDivAfterVirialRaw;
summary.lastVirialDVRms = outV.summary.lastVirialDVRms;
summary.lastRmsDivAfterVirialRaw = outV.summary.lastRmsDivAfterVirialRaw;
end

function r = safe_ratio(a, b)
r = a / max(b, eps);
end

function write_summary_file(cmp)
outDir = cmp.params.outputDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fname = fullfile(outDir, 'density_projection_virial_summary.txt');
fid = fopen(fname, 'w');
if fid < 0
    warning('Could not write summary file: %s', fname);
    return;
end
cleaner = onCleanup(@() fclose(fid));
S = cmp.summary;
fprintf(fid, 'Q7 density projection + virial comparison summary\n');
fprintf(fid, '================================================\n\n');
fprintf(fid, 'seed = %d\n', cmp.params.seed);
fprintf(fid, 'Nx = %d\nNy = %d\ngamma = %.16g\nnSteps = %d\nsampleEvery = %d\n', ...
    cmp.params.Nx, cmp.params.Ny, cmp.params.gamma, cmp.params.nSteps, cmp.params.sampleEvery);
fprintf(fid, 'initialPopulationMode = %s\n', char(cmp.params.initialPopulationMode));
fprintf(fid, 'virialDensityKickStrength = %.16g\nvirialK = %.16g\nvirialSmoothPasses = %d\n', ...
    cmp.params.virialDensityKickStrength, cmp.params.virialK, cmp.params.virialSmoothPasses);
fprintf(fid, 'meanStdClassic = %.16g\nmeanStdProjected = %.16g\nmeanStdVirial = %.16g\n', ...
    S.meanStdClassic, S.meanStdProjected, S.meanStdVirial);
fprintf(fid, 'ratioMeanStdProjectedClassic = %.16g\nratioMeanStdVirialClassic = %.16g\nratioMeanStdVirialProjected = %.16g\n', ...
    S.ratioMeanStdProjectedClassic, S.ratioMeanStdVirialClassic, S.ratioMeanStdVirialProjected);
fprintf(fid, 'meanLowKClassic = %.16g\nmeanLowKProjected = %.16g\nmeanLowKVirial = %.16g\n', ...
    S.meanLowKClassic, S.meanLowKProjected, S.meanLowKVirial);
fprintf(fid, 'meanRhoTransportClassic = %.16g\nmeanRhoTransportProjected = %.16g\nmeanRhoTransportVirial = %.16g\n', ...
    S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.meanRhoTransportVirial);
fprintf(fid, 'nuEffClassic = %.16g\nnuEffProjected = %.16g\nnuEffVirial = %.16g\n', ...
    S.nuEffClassic, S.nuEffProjected, S.nuEffVirial);
fprintf(fid, 'R2Classic = %.16g\nR2Projected = %.16g\nR2Virial = %.16g\n', ...
    S.R2Classic, S.R2Projected, S.R2Virial);
fprintf(fid, 'meanVirialDVRms = %.16g\nmeanRmsDivAfterVirialRaw = %.16g\n', ...
    S.meanVirialDVRms, S.meanRmsDivAfterVirialRaw);
fprintf('summary file written: %s\n', fname);
end
