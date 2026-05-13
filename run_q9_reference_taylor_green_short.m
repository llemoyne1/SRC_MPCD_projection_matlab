%RUN_Q9_REFERENCE_TAYLOR_GREEN_SHORT Reference Taylor-Green Q6/Q9 comparison.
%
% This is the first structure-preservation test after Poiseuille and piston:
% a fully periodic Taylor-Green vortex with no walls and no obstacle.

clearvars
clc

outputDir = sprintf('q9_reference_taylor_green_short_%s', datestr(now, 'yyyymmdd_HHMMSS'));
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
logFile = fullfile(outputDir, 'console_log.txt');
try
    diary(logFile);
catch ME
    warning('Could not start diary at %s: %s', logFile, ME.message);
end
cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('=== Q9 reference Taylor-Green short: Q6 vs Q9 ===\n');

params = struct();
params.Lx = 1.0;
params.Ly = 1.0;
params.Nx = 32;
params.Ny = 32;
params.gamma = 20;
params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;
params.taylorGreenAmplitude = 0.30;
params.taylorGreenModeX = 1;
params.taylorGreenModeY = 1;
params.taylorGreenThermalNoise = true;

params.dt = 2.0e-3;
params.nSteps = 500;
params.sampleEvery = 25;
params.progressEvery = 1000;
params.kBT = 0.01;
params.alphaDeg = 90;
params.bodyForceX = 0.0;
params.bodyForceY = 0.0;
params.useRandomGridShift = true;

params.projectionStrength = 1.0;
params.projectionInterpolationMethod = 'nearest';
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;

params.thermostatAfterProjection = true;
params.thermostatTargetKBT = params.kBT;
params.thermostatStrength = 1.0;
params.thermostatMinParticlesPerCell = 3;

params.densityBandFraction = 0.20;
params.lowKMaxIndex = 2;
params.storeDensityMaps = true;
params.storeVorticityMaps = true;
params.computeFullDiagnosticsEveryStep = false;
params.makeFigures = true;
params.makeComparisonFigures = true;
params.maxWallClockSeconds = Inf;

cmp = run_compare_projection_taylor_green_q6_q9(params);
summaryTable = cmp.summaryTable;

matFile = fullfile(outputDir, 'q9_reference_taylor_green_short.mat');
csvFile = fullfile(outputDir, 'q9_reference_taylor_green_short_summary.csv');
txtFile = fullfile(outputDir, 'q9_reference_taylor_green_short_summary.txt');

save(matFile, 'params', 'cmp', 'summaryTable', '-v7.3');
writetable(summaryTable, csvFile);
write_summary_txt(txtFile, params, summaryTable, matFile, csvFile);

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);

function write_summary_txt(txtFile, params, T, matFile, csvFile)
fid = fopen(txtFile, 'w');
if fid < 0
    error('Could not open summary file for writing: %s', txtFile);
end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '=== Q9 reference Taylor-Green short: Q6 vs Q9 ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'steps/sampleEvery                 : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %.12g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'initialPopulationMode             : %s\n', params.initialPopulationMode);
fprintf(fid, 'TG amplitude/mode                 : %.12g / (%d,%d)\n', params.taylorGreenAmplitude, params.taylorGreenModeX, params.taylorGreenModeY);
fprintf(fid, 'dt, kBT, alphaDeg                 : %.12g, %.12g, %.12g\n', params.dt, params.kBT, params.alphaDeg);

fprintf(fid, '\n--- Q9 parameters ---\n');
fprintf(fid, 'projectionStrength                : %.12g\n', params.projectionStrength);
fprintf(fid, 'massFluxProjectionMode            : %s\n', params.massFluxProjectionMode);
fprintf(fid, 'massFluxProjectionStrength        : %.12g\n', params.massFluxProjectionStrength);
fprintf(fid, 'massFluxDensityRelaxationBeta     : %.12g\n', params.massFluxDensityRelaxationBeta);
fprintf(fid, 'massFluxApplyAfterVelocityProjection : %d\n', params.massFluxApplyAfterVelocityProjection);
fprintf(fid, 'massFluxFinalVelocityProjectionCleanup : %d\n', params.massFluxFinalVelocityProjectionCleanup);
fprintf(fid, 'massFluxFinalVelocityProjectionStrength : %.12g\n', params.massFluxFinalVelocityProjectionStrength);
fprintf(fid, 'massFluxTargetFilter              : %s\n', params.massFluxTargetFilter);
fprintf(fid, 'massFluxLowKMaxIndex              : %d\n', params.massFluxLowKMaxIndex);

fprintf(fid, '\n--- Summary table ---\n');
for i = 1:height(T)
    fprintf(fid, '\nmethod: %s\n', string(T.method(i)));
    fprintf(fid, 'actual last step / stopped early  : %d / %d\n', T.actualLastStep(i), T.stoppedEarly(i));
    fprintf(fid, 'mean std(N)                       : %.12g\n', T.meanStdN(i));
    fprintf(fid, 'mean out-band                     : %.12g\n', T.meanOutBand(i));
    fprintf(fid, 'mean low-k density energy         : %.12g\n', T.meanLowK(i));
    fprintf(fid, 'time-avg rel RMS                  : %.12g\n', T.timeAvgRelRms(i));
    fprintf(fid, 'mean div(u) after particles       : %.12g\n', T.meanDivUAfterParticles(i));
    fprintf(fid, 'mean mass-flux residual           : %.12g\n', T.meanMassFluxResidual(i));
    fprintf(fid, 'mean mass-flux after particles    : %.12g\n', T.meanMassFluxAfterParticles(i));
    fprintf(fid, 'mean/final TG amplitude           : %.12g / %.12g\n', T.meanTGAmplitude(i), T.finalTGAmplitude(i));
    fprintf(fid, 'amplitude retention               : %.12g\n', T.amplitudeRetention(i));
    fprintf(fid, 'mean/final TG mode energy         : %.12g / %.12g\n', T.meanTGModeEnergy(i), T.finalTGModeEnergy(i));
    fprintf(fid, 'mode energy retention             : %.12g\n', T.modeEnergyRetention(i));
    fprintf(fid, 'mean/final TG coherence           : %.12g / %.12g\n', T.meanTGCoherence(i), T.finalTGCoherence(i));
    fprintf(fid, 'mean/final enstrophy              : %.12g / %.12g\n', T.meanEnstrophy(i), T.finalEnstrophy(i));
    fprintf(fid, 'mean/final high-k vel fraction    : %.12g / %.12g\n', T.meanHighKVelocityFraction(i), T.finalHighKVelocityFraction(i));
    fprintf(fid, 'mean/final kBT cell               : %.12g / %.12g\n', T.meanKBTCell(i), T.finalKBTCell(i));
end

fprintf(fid, '\n--- Ratios Q9/Q6 ---\n');
fprintf(fid, 'std(N)                            : %.12g\n', T.ratio_meanStdN_vs_Q6(2));
fprintf(fid, 'out-band                          : %.12g\n', T.ratio_meanOutBand_vs_Q6(2));
fprintf(fid, 'low-k density energy              : %.12g\n', T.ratio_meanLowK_vs_Q6(2));
fprintf(fid, 'div(u) after particles            : %.12g\n', T.ratio_meanDivUAfterParticles_vs_Q6(2));
fprintf(fid, 'final TG amplitude                : %.12g\n', T.ratio_finalTGAmplitude_vs_Q6(2));
fprintf(fid, 'amplitude retention               : %.12g\n', T.ratio_amplitudeRetention_vs_Q6(2));
fprintf(fid, 'final TG mode energy              : %.12g\n', T.ratio_finalTGModeEnergy_vs_Q6(2));
fprintf(fid, 'mode energy retention             : %.12g\n', T.ratio_modeEnergyRetention_vs_Q6(2));
fprintf(fid, 'final TG coherence                : %.12g\n', T.ratio_finalTGCoherence_vs_Q6(2));
fprintf(fid, 'final enstrophy                   : %.12g\n', T.ratio_finalEnstrophy_vs_Q6(2));
fprintf(fid, 'final high-k vel fraction         : %.12g\n', T.ratio_finalHighKVelocityFraction_vs_Q6(2));

fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', txtFile);
end
