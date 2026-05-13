%RUN_Q9_REFERENCE_STEP_SHORT Reference grid-aligned step Q6/Q9 comparison.
%
% This is the first separated-flow structure test after Taylor-Green.  The
% step is rectangular and grid-aligned to avoid the curved-wall difficulty of
% the cylinder/von-Karman case.

clearvars
%clc

outputDir = sprintf('q9_reference_step_short_%s', datestr(now, 'yyyymmdd_HHMMSS'));
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

fprintf('=== Q9 reference step short: Q6 vs Q9 ===\n');

params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Nx = 48;
params.Ny = 24;
params.gamma = 20;
params.seed = 11;
params.initialPopulationMode = 'exact_per_fluid_cell';
params.initialVelocityZeroGlobalMean = true;
params.initialMeanVelocityX = 0.04;
params.initialMeanVelocityY = 0.0;

params.stepX0 = 0.15 * params.Lx;
params.stepX1 = 0.35 * params.Lx;
params.stepHeight = 0.25 * params.Ly;
params.stepWallMode = 'bounceback';
params.wallModeY = 'bounceback';

params.dt = 1.0e-3;
params.nSteps = 10000;
params.sampleEvery = 50;
params.progressEvery = 100;
params.kBT = 0.02;
params.alphaDeg = 90;
params.bodyForceX = 0.0;
params.bodyForceY = 0.0;
params.useRandomGridShiftX = true;
params.useRandomGridShiftY = false;

params.meanFlowControlMode = 'relax_to_target';
params.targetMeanVelocityX = params.initialMeanVelocityX;
params.meanFlowRelaxationTau = 0.5;
params.meanFlowCorrectionMax = 5.0e-5;
params.meanFlowControlRemoveMeanY = true;

params.projectionStrength = 1.0;
params.projectionInterpolationMethod = 'nearest';
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 5e-4;
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
params.abortOnUnstable = true;
params.maxStableAbsMeanUx = 0.5;
params.maxStableAbsMeanUy = 0.5;
params.maxStableMeanEnstrophy = 5.0e3;

cmp = run_compare_projection_step_q6_q9(params);
summaryTable = cmp.summaryTable;

matFile = fullfile(outputDir, 'q9_reference_step_short.mat');
csvFile = fullfile(outputDir, 'q9_reference_step_short_summary.csv');
txtFile = fullfile(outputDir, 'q9_reference_step_short_summary.txt');

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

fprintf(fid, '=== Q9 reference step short: Q6 vs Q9 ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'steps/sampleEvery                 : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %.12g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'initialPopulationMode             : %s\n', params.initialPopulationMode);
fprintf(fid, 'initialMeanVelocityX              : %.12g\n', params.initialMeanVelocityX);
fprintf(fid, 'step x0/x1/height/mode            : %.12g / %.12g / %.12g / %s\n', ...
    params.stepX0, params.stepX1, params.stepHeight, params.stepWallMode);
fprintf(fid, 'mean-flow control mode/target/tau : %s / %.12g / %.12g\n', ...
    params.meanFlowControlMode, params.targetMeanVelocityX, params.meanFlowRelaxationTau);
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
    fprintf(fid, 'mean rho transport                : %.12g\n', T.meanRhoTransport(i));
    fprintf(fid, 'mean div(u) after particles       : %.12g\n', T.meanDivUAfter(i));
    fprintf(fid, 'mean mass-flux residual           : %.12g\n', T.meanMassFluxResidual(i));
    fprintf(fid, 'mean mass-flux after particles    : %.12g\n', T.meanMassFluxAfter(i));
    fprintf(fid, 'mean/final Ux                     : %.12g / %.12g\n', T.meanUx(i), T.finalUx(i));
    fprintf(fid, 'mean/final enstrophy              : %.12g / %.12g\n', T.meanEnstrophy(i), T.finalEnstrophy(i));
    fprintf(fid, 'mean/final shear omega RMS        : %.12g / %.12g\n', T.meanShearOmegaRms(i), T.finalShearOmegaRms(i));
    fprintf(fid, 'mean/final recirc length          : %.12g / %.12g\n', T.meanRecirculationLength(i), T.finalRecirculationLength(i));
    fprintf(fid, 'mean/final recirc area            : %.12g / %.12g\n', T.meanRecirculationArea(i), T.finalRecirculationArea(i));
    fprintf(fid, 'mean recirc min Ux                : %.12g\n', T.meanRecirculationMinUx(i));
    fprintf(fid, 'probe omega/Uy RMS                : %.12g / %.12g\n', T.probeOmegaRms(i), T.probeUyRms(i));
    fprintf(fid, 'mean step hits                    : %.12g\n', T.meanStepHits(i));
    fprintf(fid, 'mean/final kBT cell               : %.12g / %.12g\n', T.meanKBTCell(i), T.finalKBTCell(i));
end

fprintf(fid, '\n--- Ratios Q9/Q6 ---\n');
fprintf(fid, 'std(N)                            : %.12g\n', T.ratioStd_vs_Q6(2));
fprintf(fid, 'out-band                          : %.12g\n', T.ratioOutBand_vs_Q6(2));
fprintf(fid, 'low-k density energy              : %.12g\n', T.ratioLowK_vs_Q6(2));
fprintf(fid, 'rho transport                     : %.12g\n', T.ratioRhoTransport_vs_Q6(2));
fprintf(fid, 'div(u) after particles            : %.12g\n', T.ratioDivUAfter_vs_Q6(2));
fprintf(fid, 'mean enstrophy                    : %.12g\n', T.ratioEnstrophy_vs_Q6(2));
fprintf(fid, 'mean shear omega RMS              : %.12g\n', T.ratioShearOmegaRms_vs_Q6(2));
fprintf(fid, 'mean recirc length                : %.12g\n', T.ratioRecircLength_vs_Q6(2));
fprintf(fid, 'mean recirc area                  : %.12g\n', T.ratioRecircArea_vs_Q6(2));
fprintf(fid, 'probe omega RMS                   : %.12g\n', T.ratioProbeOmegaRms_vs_Q6(2));
fprintf(fid, 'probe Uy RMS                      : %.12g\n', T.ratioProbeUyRms_vs_Q6(2));

fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', txtFile);
end
