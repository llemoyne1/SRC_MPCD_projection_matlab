%RUN_Q9_REFERENCE_STEP_VISUAL_DEMO Live visualization for validated Q9 step-channel case.
%
% This script is intended for qualitative inspection/illustration.  It uses
% the validated step-channel Q9 parameters and displays particles, density,
% velocity/quiver, and vorticity during the run.
%
% To compare Q6 and Q9 live, set runMode = 'compare'.  For a lighter run,
% keep the default runMode = 'q9_only'.

clearvars
clc

runMode = 'q9_only';   % 'q9_only' or 'compare'

outputDir = sprintf('q9_reference_step_visual_%s', datestr(now, 'yyyymmdd_HHMMSS'));
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

fprintf('=== Q9 reference step visual demo ===\n');

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
params.nSteps = 3000;
params.sampleEvery = 25;
params.progressEvery = 250;
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
params.massFluxDensityRelaxationBeta = 5.0e-4;
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
params.storeVelocityMaps = true;
params.computeFullDiagnosticsEveryStep = false;
params.makeFigures = true;
params.maxWallClockSeconds = Inf;
params.abortOnUnstable = true;
params.maxStableAbsMeanUx = 0.5;
params.maxStableAbsMeanUy = 0.5;
params.maxStableMeanEnstrophy = 5.0e3;

% Live visualization controls.
params.visualEnable = true;
params.visualEvery = 25;
params.visualFigureId = 200;
params.visualFigureName = 'Step-channel live Q9';
params.visualTitleSuffix = 'Q9 validated step';
params.visualMaxParticles = 3500;
params.visualParticleSize = 5;
params.visualDensityCLim = 0.5;
params.visualOmegaCLim = NaN;       % auto symmetric clipping
params.visualOmegaClipSigma = 3.0;
params.visualSpeedCLim = NaN;       % auto
params.visualQuiverStrideX = 4;
params.visualQuiverStrideY = 3;
params.visualQuiverScale = 1.5;
params.visualPause = 0.01;
params.visualSaveFrames = false;    % true to export PNG frames
params.visualFrameDir = fullfile(outputDir, 'frames');
params.visualFramePrefix = 'step_q9';

switch lower(runMode)
    case 'q9_only'
        outQ9 = run_projection_step_channel_demo(params); %#ok<NASGU>
        matFile = fullfile(outputDir, 'q9_reference_step_visual.mat');
        save(matFile, 'params', 'outQ9', '-v7.3');
        write_visual_summary(fullfile(outputDir, 'q9_reference_step_visual_summary.txt'), params, outQ9, matFile);
    case 'compare'
        params.makeComparisonFigures = true;
        cmp = run_compare_projection_step_q6_q9(params); %#ok<NASGU>
        matFile = fullfile(outputDir, 'q9_reference_step_visual_compare.mat');
        save(matFile, 'params', 'cmp', '-v7.3');
        fprintf('Comparison saved to:\n%s\n', matFile);
    otherwise
        error('Unknown runMode: %s', runMode);
end

fprintf('Output directory:\n%s\n', outputDir);

function write_visual_summary(txtFile, params, out, matFile)
fid = fopen(txtFile, 'w');
if fid < 0
    warning('Could not write summary file: %s', txtFile);
    return;
end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
s = out.summary;
fprintf(fid, '=== Q9 step visual demo ===\n');
fprintf(fid, 'steps/sampleEvery/visualEvery       : %d / %d / %d\n', params.nSteps, params.sampleEvery, params.visualEvery);
fprintf(fid, 'grid, gamma                         : %d x %d, %.12g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'step x0/x1/height/mode              : %.12g / %.12g / %.12g / %s\n', ...
    params.stepX0, params.stepX1, params.stepHeight, params.stepWallMode);
fprintf(fid, 'massFlux beta / lowKMax / cleanup   : %.12g / %d / %.12g\n', ...
    params.massFluxDensityRelaxationBeta, params.massFluxLowKMaxIndex, ...
    params.massFluxFinalVelocityProjectionStrength * double(params.massFluxFinalVelocityProjectionCleanup));
fprintf(fid, 'actual last step / stopped early    : %d / %d\n', s.actualLastStep, s.stoppedEarly);
fprintf(fid, 'mean low-k density energy           : %.12g\n', s.meanLowKEnergy);
fprintf(fid, 'mean div(u) after particles         : %.12g\n', s.meanRmsDivParticleAfter);
fprintf(fid, 'mean/final enstrophy                : %.12g / %.12g\n', s.meanEnstrophy, s.finalEnstrophy);
fprintf(fid, 'mean/final shear omega RMS          : %.12g / %.12g\n', s.meanShearOmegaRms, s.finalShearOmegaRms);
fprintf(fid, 'mean/final recirculation length     : %.12g / %.12g\n', s.meanRecirculationLength, s.finalRecirculationLength);
fprintf(fid, 'MAT file                            : %s\n', matFile);
end
