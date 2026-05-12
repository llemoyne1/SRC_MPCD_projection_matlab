%% run_q9_reference_vonkarman_short.m
% First Q9-native cylinder/von-Karman reference check.
%
% This script compares Q6 and Q9 on a fixed circular obstacle in an
% x-periodic channel.  The goal is not yet a calibrated Strouhal/Re study;
% this is a non-regression test for density low-k control versus preservation
% of vorticity and wake fluctuations.

clear functions
close all
% clc

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q9_reference_vonkarman_short_' tag];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on
cleanupObj = onCleanup(@() diary('off'));

fprintf('\n=== Q9 reference von Karman short: Q6 vs Q9 ===\n');

%% === Parameters ===
params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Nx = 48;
params.Ny = 24;
params.gamma = 10;

params.nSteps = 5000;
params.sampleEvery = 50;
params.progressEvery = 500;
params.maxWallClockSeconds = Inf;

params.dt = 1.0e-3;
params.kBT = 0.2;
params.alphaDeg = 90;
params.bodyForceX = 0.02;
params.bodyForceY = 0.0;

params.wallModeY = 'thermalize';
params.thermostatAfterProjection = true;

params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;
params.initialMeanVelocityX = 0.08;
params.initialMeanVelocityY = 0.0;

%% === Cylinder geometry ===
params.cylinderCenterX = 0.35 * params.Lx;
params.cylinderCenterY = 0.50 * params.Ly;
params.cylinderRadius = 0.10 * params.Ly;
params.cylinderWallMode = 'bounceback';
params.wakeProbeDistanceDiameters = 2.0;
params.wakeProbeYOffsetDiameters = 0.35;
params.sheddingTransientFraction = 0.40;
params.targetStrouhal = 0.18;
params.targetRe = NaN;
params.nuEffForRe = NaN;  % set to calibrated nu_eff to print Re = U D / nu

%% === Q9 method ===
params.projectionStrength = 1.0;
params.projectedStrength = 1.0;
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;

%% === Diagnostics ===
params.storeDensityMaps = true;
params.storeVorticityMaps = true;
params.densityBandFraction = 0.20;
params.densityExcludeWallCells = 0;
params.lowKMaxIndex = 2;
params.makeFigures = true;

%% === Run comparison ===
tic
cmp = run_compare_projection_cylinder_q6_q9(params);
elapsedScript = toc;
summary = cmp.summary;
summary.elapsedScript = repmat(elapsedScript, height(summary), 1);

%% === Save data ===
matFile = fullfile(outputDir, 'q9_reference_vonkarman_short.mat');
csvFile = fullfile(outputDir, 'q9_reference_vonkarman_short_summary.csv');
txtFile = fullfile(outputDir, 'q9_reference_vonkarman_short_summary.txt');

save(matFile, 'params', 'cmp', 'summary', '-v7.3');
writetable(summary, csvFile);

fid = fopen(txtFile, 'w');
fprintf(fid, '=== Q9 reference von Karman short: Q6 vs Q9 ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'steps/sampleEvery                 : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'initialPopulationMode             : %s\n', params.initialPopulationMode);
fprintf(fid, 'initialMeanVelocityX/bodyForceX   : %.12g / %.12g\n', params.initialMeanVelocityX, params.bodyForceX);
fprintf(fid, 'cylinder center/radius/mode       : (%.12g, %.12g) / %.12g / %s\n', ...
    params.cylinderCenterX, params.cylinderCenterY, params.cylinderRadius, params.cylinderWallMode);
fprintf(fid, 'wake probe x/D, y/D offset        : %.12g / %.12g\n', ...
    params.wakeProbeDistanceDiameters, params.wakeProbeYOffsetDiameters);
fprintf(fid, 'shedding transient fraction       : %.12g\n', params.sheddingTransientFraction);
fprintf(fid, 'nuEffForRe / targetRe / targetSt  : %.12g / %.12g / %.12g\n', ...
    params.nuEffForRe, params.targetRe, params.targetStrouhal);

fprintf(fid, '\n--- Q9 parameters ---\n');
fprintf(fid, 'projectionStrength                : %.6g\n', params.projectionStrength);
fprintf(fid, 'massFluxProjectionMode            : %s\n', params.massFluxProjectionMode);
fprintf(fid, 'massFluxProjectionStrength        : %.6g\n', params.massFluxProjectionStrength);
fprintf(fid, 'massFluxDensityRelaxationBeta     : %.6g\n', params.massFluxDensityRelaxationBeta);
fprintf(fid, 'massFluxApplyAfterVelocityProjection : %d\n', params.massFluxApplyAfterVelocityProjection);
fprintf(fid, 'massFluxFinalVelocityProjectionCleanup : %d\n', params.massFluxFinalVelocityProjectionCleanup);
fprintf(fid, 'massFluxFinalVelocityProjectionStrength : %.6g\n', params.massFluxFinalVelocityProjectionStrength);
fprintf(fid, 'massFluxTargetFilter              : %s\n', params.massFluxTargetFilter);
fprintf(fid, 'massFluxLowKMaxIndex              : %d\n', params.massFluxLowKMaxIndex);

fprintf(fid, '\n--- Summary table ---\n');
for i = 1:height(summary)
    fprintf(fid, '\nmethod: %s\n', char(summary.method(i)));
    fprintf(fid, 'mean std(N)                       : %.12g\n', summary.meanStdN(i));
    fprintf(fid, 'mean out-band                     : %.12g\n', summary.meanOutBand(i));
    fprintf(fid, 'mean low-k energy                 : %.12g\n', summary.meanLowK(i));
    fprintf(fid, 'time-avg rel RMS                  : %.12g\n', summary.timeAvgRelRms(i));
    fprintf(fid, 'mean rho transport                : %.12g\n', summary.meanRhoTransport(i));
    fprintf(fid, 'mean div(u) after particles       : %.12g\n', summary.meanDivUAfter(i));
    fprintf(fid, 'mean mass-flux residual           : %.12g\n', summary.meanMassFluxResidual(i));
    fprintf(fid, 'mean mass-flux after particles    : %.12g\n', summary.meanMassFluxParticleAfter(i));
    fprintf(fid, 'mean Ux                           : %.12g\n', summary.meanUx(i));
    fprintf(fid, 'final Ux                          : %.12g\n', summary.finalUx(i));
    fprintf(fid, 'mean enstrophy                    : %.12g\n', summary.meanEnstrophy(i));
    fprintf(fid, 'final enstrophy                   : %.12g\n', summary.finalEnstrophy(i));
    fprintf(fid, 'wake omega RMS                    : %.12g\n', summary.wakeOmegaRms(i));
    fprintf(fid, 'wake Uy RMS                       : %.12g\n', summary.wakeUyRms(i));
    fprintf(fid, 'wake omega upper RMS              : %.12g\n', summary.wakeOmegaUpperRms(i));
    fprintf(fid, 'wake Uy upper RMS                 : %.12g\n', summary.wakeUyUpperRms(i));
    fprintf(fid, 'wake omega antisym RMS            : %.12g\n', summary.wakeOmegaAntiSymRms(i));
    fprintf(fid, 'wake Uy antisym RMS               : %.12g\n', summary.wakeUyAntiSymRms(i));
    fprintf(fid, 'wake omega sign changes           : %.12g\n', summary.wakeOmegaSignChanges(i));
    fprintf(fid, 'wake Uy upper sign changes        : %.12g\n', summary.wakeUyUpperSignChanges(i));
    fprintf(fid, 'wake omega antisym sign changes   : %.12g\n', summary.wakeOmegaAntiSymSignChanges(i));
    fprintf(fid, 'shedding f, St from upper Uy      : %.12g / %.12g\n', ...
        summary.sheddingFrequencyWakeUyUpper(i), summary.strouhalWakeUyUpper(i));
    fprintf(fid, 'shedding SNR upper Uy             : %.12g\n', summary.sheddingSNRWakeUyUpper(i));
    fprintf(fid, 'shedding f, St omega antisym      : %.12g / %.12g\n', ...
        summary.sheddingFrequencyWakeOmegaAntiSym(i), summary.strouhalWakeOmegaAntiSym(i));
    fprintf(fid, 'shedding SNR omega antisym        : %.12g\n', summary.sheddingSNROmegaAntiSym(i));
    fprintf(fid, 'estimated Re mean/final           : %.12g / %.12g\n', summary.reEffMeanUx(i), summary.reEffFinalUx(i));
    fprintf(fid, 'mean cylinder hits                : %.12g\n', summary.meanCylinderHits(i));
    fprintf(fid, 'mean kBT cell                     : %.12g\n', summary.meanKBTCell(i));
end

fprintf(fid, '\n--- Ratios Q9/Q6 ---\n');
fprintf(fid, 'std(N)                            : %.12g\n', summary.ratioStd_vs_Q6(2));
fprintf(fid, 'out-band                          : %.12g\n', summary.ratioOutBand_vs_Q6(2));
fprintf(fid, 'low-k energy                      : %.12g\n', summary.ratioLowK_vs_Q6(2));
fprintf(fid, 'rho transport                     : %.12g\n', summary.ratioRhoTransport_vs_Q6(2));
fprintf(fid, 'div(u) after particles            : %.12g\n', summary.ratioDivUAfter_vs_Q6(2));
fprintf(fid, 'enstrophy                         : %.12g\n', summary.ratioEnstrophy_vs_Q6(2));
fprintf(fid, 'wake omega RMS                    : %.12g\n', summary.ratioWakeOmegaRms_vs_Q6(2));
fprintf(fid, 'wake Uy RMS                       : %.12g\n', summary.ratioWakeUyRms_vs_Q6(2));
fprintf(fid, 'wake omega upper RMS              : %.12g\n', summary.ratioWakeOmegaUpperRms_vs_Q6(2));
fprintf(fid, 'wake Uy upper RMS                 : %.12g\n', summary.ratioWakeUyUpperRms_vs_Q6(2));
fprintf(fid, 'wake omega antisym RMS            : %.12g\n', summary.ratioWakeOmegaAntiSymRms_vs_Q6(2));
fprintf(fid, 'wake Uy antisym RMS               : %.12g\n', summary.ratioWakeUyAntiSymRms_vs_Q6(2));

fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', txtFile);
fclose(fid);

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);

diary off
