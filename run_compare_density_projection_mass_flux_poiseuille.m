function cmp = run_compare_density_projection_mass_flux_poiseuille(params)
%RUN_COMPARE_DENSITY_PROJECTION_MASS_FLUX_POISEUILLE Compare Q6 velocity projection and Q8 mass-flux projection.
%
%   cmp = run_compare_density_projection_mass_flux_poiseuille(params)
%
% Runs three paired Poiseuille simulations with the same seed:
%   C : SRC classic, projectionStrength = 0
%   P : Q6 velocity projection, div(u) = 0
%   M : Q8 mass-flux projection, div(N u) = target
%
% Q8 modes:
%   params.massFluxProjectionMode = 'conservative'      -> target = 0
%   params.massFluxProjectionMode = 'relax_to_uniform'  -> target = beta/dt*(N-gamma)

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
classicParams.massFluxProjectionMode = 'off';
classicParams.seed = params.seed;

projectedParams = base;
projectedParams.projectionEnable = true;
projectedParams.projectionStrength = params.projectedStrength;
projectedParams.massFluxProjectionMode = 'off';
projectedParams.seed = params.seed;

massParams = base;
massParams.projectionEnable = true;
if params.massFluxApplyAfterVelocityProjection
    massParams.projectionStrength = params.projectedStrength;
else
    massParams.projectionStrength = 0.0; % pure mass-flux branch
end
massParams.massFluxProjectionMode = params.massFluxProjectionMode;
massParams.massFluxProjectionStrength = params.massFluxProjectionStrength;
massParams.massFluxDensityRelaxationBeta = params.massFluxDensityRelaxationBeta;
massParams.massFluxApplyAfterVelocityProjection = params.massFluxApplyAfterVelocityProjection;
massParams.massFluxTargetFilter = params.massFluxTargetFilter;
massParams.massFluxLowKMaxIndex = params.massFluxLowKMaxIndex;
massParams.seed = params.seed;

fprintf('\n=== Q8 density comparison: classic projectionStrength=0 ===\n');
outClassic = run_projection_poiseuille_long_demo(classicParams);

fprintf('\n=== Q8 density comparison: velocity projection strength=%.6g ===\n', projectedParams.projectionStrength);
outProjected = run_projection_poiseuille_long_demo(projectedParams);

fprintf('\n=== Q8 density comparison: mass-flux projection mode=%s strength=%.6g beta=%.6g ===\n', ...
    char(massParams.massFluxProjectionMode), massParams.massFluxProjectionStrength, massParams.massFluxDensityRelaxationBeta);
outMassFlux = run_projection_poiseuille_long_demo(massParams);

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
metricsMassFlux = projection_density_homogeneity_metrics(outMassFlux.NMaps, massParams, ...
    'sampleTimes', outMassFlux.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex);

summary = compare_summaries(metricsClassic, metricsProjected, metricsMassFlux, outClassic, outProjected, outMassFlux);

cmp = struct();
cmp.params = params;
cmp.classic = outClassic;
cmp.projected = outProjected;
cmp.massFlux = outMassFlux;
cmp.metricsClassic = metricsClassic;
cmp.metricsProjected = metricsProjected;
cmp.metricsMassFlux = metricsMassFlux;
cmp.summary = summary;

fprintf('\n=== Q8 density / mass-flux homogeneity summary ===\n');
fprintf('seed                                      : %d\n', params.seed);
fprintf('steps/sampleEvery                        : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf('grid, gamma                              : %d x %d, %.6g\n', params.Nx, params.Ny, params.gamma);
fprintf('initialPopulationMode                    : %s\n', char(params.initialPopulationMode));
fprintf('mass-flux mode, strength, beta           : %s / %.6g / %.6g\n', ...
    char(params.massFluxProjectionMode), params.massFluxProjectionStrength, params.massFluxDensityRelaxationBeta);
fprintf('mass-flux hybrid/filter/lowK             : afterVelocity=%d / %s / %d\n', ...
    params.massFluxApplyAfterVelocityProjection, char(params.massFluxTargetFilter), params.massFluxLowKMaxIndex);
fprintf('initial std C/P/M                        : %.6g / %.6g / %.6g\n', ...
    summary.initialStdClassic, summary.initialStdProjected, summary.initialStdMassFlux);
fprintf('mean std(N) C/P/M                        : %.6g / %.6g / %.6g\n', ...
    summary.meanStdClassic, summary.meanStdProjected, summary.meanStdMassFlux);
fprintf('mean std ratios P/C, M/C, M/P            : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanStdProjectedClassic, summary.ratioMeanStdMassFluxClassic, summary.ratioMeanStdMassFluxProjected);
fprintf('mean out-band C/P/M                      : %.6g / %.6g / %.6g\n', ...
    summary.meanOutBandClassic, summary.meanOutBandProjected, summary.meanOutBandMassFlux);
fprintf('mean low-k C/P/M                         : %.6g / %.6g / %.6g\n', ...
    summary.meanLowKClassic, summary.meanLowKProjected, summary.meanLowKMassFlux);
fprintf('low-k ratios P/C, M/C, M/P               : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanLowKProjectedClassic, summary.ratioMeanLowKMassFluxClassic, summary.ratioMeanLowKMassFluxProjected);
fprintf('mean rho transport C/P/M                 : %.6g / %.6g / %.6g\n', ...
    summary.meanRhoTransportClassic, summary.meanRhoTransportProjected, summary.meanRhoTransportMassFlux);
fprintf('rho transport ratios P/C, M/C, M/P       : %.6g / %.6g / %.6g\n', ...
    summary.ratioMeanRhoTransportProjectedClassic, summary.ratioMeanRhoTransportMassFluxClassic, summary.ratioMeanRhoTransportMassFluxProjected);
fprintf('time-avg rel RMS C/P/M                   : %.6g / %.6g / %.6g\n', ...
    summary.timeAvgRelRmsClassic, summary.timeAvgRelRmsProjected, summary.timeAvgRelRmsMassFlux);
fprintf('nu_eff C/P/M                             : %.6g / %.6g / %.6g\n', ...
    summary.nuEffClassic, summary.nuEffProjected, summary.nuEffMassFlux);
fprintf('R2 C/P/M                                 : %.6g / %.6g / %.6g\n', ...
    summary.R2Classic, summary.R2Projected, summary.R2MassFlux);
fprintf('mean mass-flux div before/target/after/res M: %.6g / %.6g / %.6g / %.6g\n', ...
    summary.meanMassFluxDivBeforeMassFlux, summary.meanMassFluxDivTargetMassFlux, ...
    summary.meanMassFluxDivAfterMassFlux, summary.meanMassFluxDivResidualMassFlux);
fprintf('elapsed C/P/M [s]                        : %.2f / %.2f / %.2f\n', ...
    outClassic.elapsedWallClock, outProjected.elapsedWallClock, outMassFlux.elapsedWallClock);

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
params = set_default(params, 'massFluxProjectionMode', 'conservative');
params = set_default(params, 'massFluxProjectionStrength', 1.0);
params = set_default(params, 'massFluxDensityRelaxationBeta', 0.0);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', false);
params = set_default(params, 'massFluxTargetFilter', 'none');
params = set_default(params, 'massFluxLowKMaxIndex', params.lowKMaxIndex);
params = set_default(params, 'massFluxProjectionRegularization', 1e-12);
params = set_default(params, 'massFluxMinCellCount', 1.0);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'saveFigures', false);
params = set_default(params, 'outputDir', 'density_mass_flux_projection_output');
params = set_default(params, 'writeSummary', true);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function S = compare_summaries(mC, mP, mM, outC, outP, outM)
S = struct();
S.initialStdClassic = mC.summary.initialStdN;
S.initialStdProjected = mP.summary.initialStdN;
S.initialStdMassFlux = mM.summary.initialStdN;
S.meanStdClassic = mC.summary.meanStdN;
S.meanStdProjected = mP.summary.meanStdN;
S.meanStdMassFlux = mM.summary.meanStdN;
S.ratioMeanStdProjectedClassic = safe_ratio(S.meanStdProjected, S.meanStdClassic);
S.ratioMeanStdMassFluxClassic = safe_ratio(S.meanStdMassFlux, S.meanStdClassic);
S.ratioMeanStdMassFluxProjected = safe_ratio(S.meanStdMassFlux, S.meanStdProjected);
S.finalStdClassic = mC.summary.finalStdN;
S.finalStdProjected = mP.summary.finalStdN;
S.finalStdMassFlux = mM.summary.finalStdN;
S.meanOutBandClassic = mC.summary.meanOutBandFraction;
S.meanOutBandProjected = mP.summary.meanOutBandFraction;
S.meanOutBandMassFlux = mM.summary.meanOutBandFraction;
S.meanLowKClassic = mC.summary.meanLowKEnergy;
S.meanLowKProjected = mP.summary.meanLowKEnergy;
S.meanLowKMassFlux = mM.summary.meanLowKEnergy;
S.ratioMeanLowKProjectedClassic = safe_ratio(S.meanLowKProjected, S.meanLowKClassic);
S.ratioMeanLowKMassFluxClassic = safe_ratio(S.meanLowKMassFlux, S.meanLowKClassic);
S.ratioMeanLowKMassFluxProjected = safe_ratio(S.meanLowKMassFlux, S.meanLowKProjected);
S.timeAvgRelRmsClassic = mC.summary.timeAvgRelRms;
S.timeAvgRelRmsProjected = mP.summary.timeAvgRelRms;
S.timeAvgRelRmsMassFlux = mM.summary.timeAvgRelRms;
S.meanRhoTransportClassic = outC.summary.meanDensityTransportProjectedRms;
S.meanRhoTransportProjected = outP.summary.meanDensityTransportProjectedRms;
S.meanRhoTransportMassFlux = outM.summary.meanDensityTransportProjectedRms;
S.ratioMeanRhoTransportProjectedClassic = safe_ratio(S.meanRhoTransportProjected, S.meanRhoTransportClassic);
S.ratioMeanRhoTransportMassFluxClassic = safe_ratio(S.meanRhoTransportMassFlux, S.meanRhoTransportClassic);
S.ratioMeanRhoTransportMassFluxProjected = safe_ratio(S.meanRhoTransportMassFlux, S.meanRhoTransportProjected);
S.nuEffClassic = outC.viscosity.nuEff;
S.nuEffProjected = outP.viscosity.nuEff;
S.nuEffMassFlux = outM.viscosity.nuEff;
S.R2Classic = outC.viscosity.R2;
S.R2Projected = outP.viscosity.R2;
S.R2MassFlux = outM.viscosity.R2;
S.signalToNoiseClassic = outC.viscosity.signalToNoise;
S.signalToNoiseProjected = outP.viscosity.signalToNoise;
S.signalToNoiseMassFlux = outM.viscosity.signalToNoise;
S.meanMassFluxDivBeforeMassFlux = outM.summary.meanMassFluxDivBefore;
S.meanMassFluxDivTargetMassFlux = outM.summary.meanMassFluxDivTarget;
S.meanMassFluxDivAfterMassFlux = outM.summary.meanMassFluxDivProjectedAfter;
S.meanMassFluxDivResidualMassFlux = outM.summary.meanMassFluxDivResidual;
end

function r = safe_ratio(a, b)
r = a / max(b, eps);
end

function write_summary_file(cmp)
outDir = cmp.params.outputDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fname = fullfile(outDir, 'density_mass_flux_projection_summary.txt');
fid = fopen(fname, 'w');
if fid < 0
    warning('Could not write summary file: %s', fname);
    return;
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
S = cmp.summary;
fprintf(fid, 'Q8 mass-flux projection comparison summary\n');
fprintf(fid, '==========================================\n\n');
fprintf(fid, 'seed = %d\n', cmp.params.seed);
fprintf(fid, 'Nx = %d\nNy = %d\ngamma = %.16g\nnSteps = %d\nsampleEvery = %d\n', ...
    cmp.params.Nx, cmp.params.Ny, cmp.params.gamma, cmp.params.nSteps, cmp.params.sampleEvery);
fprintf(fid, 'initialPopulationMode = %s\n', char(cmp.params.initialPopulationMode));
fprintf(fid, 'massFluxProjectionMode = %s\n', char(cmp.params.massFluxProjectionMode));
fprintf(fid, 'massFluxProjectionStrength = %.16g\n', cmp.params.massFluxProjectionStrength);
fprintf(fid, 'massFluxDensityRelaxationBeta = %.16g\n', cmp.params.massFluxDensityRelaxationBeta);
fprintf(fid, 'massFluxApplyAfterVelocityProjection = %d\n', cmp.params.massFluxApplyAfterVelocityProjection);
fprintf(fid, 'massFluxTargetFilter = %s\n', char(cmp.params.massFluxTargetFilter));
fprintf(fid, 'massFluxLowKMaxIndex = %d\n\n', cmp.params.massFluxLowKMaxIndex);
fprintf(fid, 'meanStd C/P/M = %.16g / %.16g / %.16g\n', S.meanStdClassic, S.meanStdProjected, S.meanStdMassFlux);
fprintf(fid, 'meanStd ratios P/C, M/C, M/P = %.16g / %.16g / %.16g\n', ...
    S.ratioMeanStdProjectedClassic, S.ratioMeanStdMassFluxClassic, S.ratioMeanStdMassFluxProjected);
fprintf(fid, 'meanLowK C/P/M = %.16g / %.16g / %.16g\n', S.meanLowKClassic, S.meanLowKProjected, S.meanLowKMassFlux);
fprintf(fid, 'meanRhoTransport C/P/M = %.16g / %.16g / %.16g\n', ...
    S.meanRhoTransportClassic, S.meanRhoTransportProjected, S.meanRhoTransportMassFlux);
fprintf(fid, 'rhoTransport ratios P/C, M/C, M/P = %.16g / %.16g / %.16g\n', ...
    S.ratioMeanRhoTransportProjectedClassic, S.ratioMeanRhoTransportMassFluxClassic, S.ratioMeanRhoTransportMassFluxProjected);
fprintf(fid, 'nuEff C/P/M = %.16g / %.16g / %.16g\n', S.nuEffClassic, S.nuEffProjected, S.nuEffMassFlux);
fprintf(fid, 'R2 C/P/M = %.16g / %.16g / %.16g\n', S.R2Classic, S.R2Projected, S.R2MassFlux);
fprintf(fid, 'massFluxDiv before/target/after/res M = %.16g / %.16g / %.16g / %.16g\n', ...
    S.meanMassFluxDivBeforeMassFlux, S.meanMassFluxDivTargetMassFlux, ...
    S.meanMassFluxDivAfterMassFlux, S.meanMassFluxDivResidualMassFlux);
fprintf('summary file written: %s\n', fname);
end
