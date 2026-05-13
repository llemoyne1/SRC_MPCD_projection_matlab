%% run_q9_reference_piston_short.m
% First Q9-native piston reference check.
%
% This script compares Q6 and Q9 on a slow 10% moving-piston compression.
% It is intentionally shorter than the 30000-step Poiseuille reference and is
% meant as a development/non-regression case before longer piston campaigns.

clear functions
close all
%clc

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q9_reference_piston_short_' tag];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on
cleanupObj = onCleanup(@() diary('off'));

fprintf('\n=== Q9 reference piston short: Q6 vs Q9 ===\n');

%% === Parameters ===
params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Nx = 32;
params.Ny = 16;
params.gamma = 20;

params.nSteps = 30000;
params.sampleEvery = 50;
params.progressEvery = 500;
params.maxWallClockSeconds = Inf;

params.dt = 1.0e-3;
params.kBT = 1.0;
params.alphaDeg = 90;
params.bodyForceX = 0.0;
params.bodyForceY = 0.0;

params.wallModeY = 'thermalize';
params.pistonTangentialMode = 'specular';
params.thermostatAfterProjection = true;

params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';

%% === Piston kinematics ===
params.pistonY0 = params.Ly;
params.pistonCompressionTarget = 0.10;
params.pistonYMin = params.Ly * (1.0 - params.pistonCompressionTarget);
params.pistonVy = -(params.pistonY0 - params.pistonYMin) / (params.nSteps * params.dt);
params.pistonStopOnMin = true;
params.pistonAffineReposition = true;

%% === Q9 method ===
params.projectionStrength = 1.0;
params.projectedStrength = 1.0;
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxFinalVelocityProjectionCleanup = false;
params.massFluxFinalVelocityProjectionStrength = 1.0;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;

%% === Density diagnostics ===
params.storeDensityMaps = true;
params.densityBandFraction = 0.20;
params.densityExcludeWallCells = 0;
params.lowKMaxIndex = 2;
params.makeFigures = true;

%% === Run comparison ===
tic
cmp = run_compare_projection_piston_q6_q9(params);
elapsedScript = toc;
summary = cmp.summary;
summary.elapsedScript = repmat(elapsedScript, height(summary), 1);

%% === Save data ===
matFile = fullfile(outputDir, 'q9_reference_piston_short.mat');
csvFile = fullfile(outputDir, 'q9_reference_piston_short_summary.csv');
txtFile = fullfile(outputDir, 'q9_reference_piston_short_summary.txt');

save(matFile, 'params', 'cmp', 'summary', '-v7.3');
writetable(summary, csvFile);

fid = fopen(txtFile, 'w');
fprintf(fid, '=== Q9 reference piston short: Q6 vs Q9 ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'steps/sampleEvery                 : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'initialPopulationMode             : %s\n', params.initialPopulationMode);
fprintf(fid, 'piston y0/ymin/vy                 : %.12g / %.12g / %.12g\n', params.pistonY0, params.pistonYMin, params.pistonVy);
fprintf(fid, 'target compression                : %.12g\n', params.pistonCompressionTarget);

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
    fprintf(fid, 'final compression                 : %.12g\n', summary.finalCompression(i));
    fprintf(fid, 'mean std(N)                       : %.12g\n', summary.meanStdN(i));
    fprintf(fid, 'mean out-band                     : %.12g\n', summary.meanOutBand(i));
    fprintf(fid, 'mean low-k energy                 : %.12g\n', summary.meanLowK(i));
    fprintf(fid, 'time-avg rel RMS                  : %.12g\n', summary.timeAvgRelRms(i));
    fprintf(fid, 'mean rho transport                : %.12g\n', summary.meanRhoTransport(i));
    fprintf(fid, 'mean div(u) after particles       : %.12g\n', summary.meanDivUAfter(i));
    fprintf(fid, 'mean mass-flux residual           : %.12g\n', summary.meanMassFluxResidual(i));
    fprintf(fid, 'mean mass-flux after particles    : %.12g\n', summary.meanMassFluxParticleAfter(i));
    fprintf(fid, 'mean Pkin                         : %.12g\n', summary.meanPkin(i));
    fprintf(fid, 'final Pkin                        : %.12g\n', summary.finalPkin(i));
end

fprintf(fid, '\n--- Ratios Q9/Q6 ---\n');
fprintf(fid, 'std(N)                            : %.12g\n', summary.ratioStd_vs_Q6(2));
fprintf(fid, 'out-band                          : %.12g\n', summary.ratioOutBand_vs_Q6(2));
fprintf(fid, 'low-k energy                      : %.12g\n', summary.ratioLowK_vs_Q6(2));
fprintf(fid, 'div(u) after particles            : %.12g\n', summary.ratioDivUAfter_vs_Q6(2));
fprintf(fid, 'rho transport                     : %.12g\n', summary.ratioRhoTransport_vs_Q6(2));
fprintf(fid, 'mass-flux residual                : %.12g\n', summary.ratioMassFluxResidual_vs_Q6(2));
fprintf(fid, 'mass-flux after particles         : %.12g\n', summary.ratioMassFluxParticleAfter_vs_Q6(2));

fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', txtFile);
fclose(fid);

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);

diary off
