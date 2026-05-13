%% run_q9_vonkarman_re_st_sweep.m
% Q9-only exploratory sweep for cylinder wake Re/St diagnostics.
%
% This script is intentionally lighter than a full Q6/Q9 comparison. It runs
% several Q9 cylinder cases, extracts off-center wake-probe spectra, and saves
% a compact table with estimated Strouhal and optional Reynolds number.
%
% Important: Re is only printed if params.nuEffForRe is provided. For a true
% comparison with classical cylinder Re, first calibrate nu_eff with the same
% gamma, kBT, alphaDeg, projection and Q9 settings.

clear functions
close all
% clc

tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q9_vonkarman_re_st_sweep_' tag];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on
cleanupObj = onCleanup(@() diary('off'));

fprintf('\n=== Q9 von Karman Re/St exploratory sweep ===\n');

%% === Base parameters ===
base = struct();
base.Lx = 2.0;
base.Ly = 1.0;
base.Nx = 64;
base.Ny = 32;
base.gamma = 10;

base.nSteps = 20000;
base.sampleEvery = 100;
base.progressEvery = 1000;
base.maxWallClockSeconds = Inf;

base.dt = 1.0e-3;
base.kBT = 0.2;
base.alphaDeg = 90;
base.bodyForceY = 0.0;
% Long periodic-cylinder runs are run with bulk-velocity control rather than
% large constant forcing.  This avoids indefinite energy injection while
% preserving a controlled mean advection past the obstacle.
base.meanFlowControlMode = 'relax_to_target';
base.targetMeanVelocityX = 0.01;
base.meanFlowRelaxationTau = 1.0;
base.meanFlowCorrectionMax = 5.0e-6;
base.meanFlowControlRemoveMeanY = true;
base.abortOnUnstable = true;
base.maxStableAbsMeanUx = 0.30;
base.maxStableAbsMeanUy = 0.30;
base.maxStableMeanEnstrophy = 1.0e3;
base.sheddingUseStrouhalBand = true;
base.sheddingStrouhalMin = 0.05;
base.sheddingStrouhalMax = 0.50;
base.wallModeY = 'thermalize';
base.thermostatAfterProjection = true;

base.seed = 11;
base.initialPopulationMode = 'exact_per_cell';
base.initialVelocityZeroGlobalMean = true;
base.initialMeanVelocityY = 0.01;

base.cylinderCenterX = 0.35 * base.Lx;
base.cylinderCenterY = 0.50 * base.Ly;
base.cylinderWallMode = 'specular';%'bounceback';
base.wakeProbeDistanceDiameters = 2.0;
base.wakeProbeYOffsetDiameters = 0.35;
base.sheddingTransientFraction = 0.40;

% Optional physical interpretation. Leave NaN until nu_eff is calibrated for
% this exact numerical fluid. Example placeholder from a previous Poiseuille
% Q9 run would be around 0.0126, but that was not calibrated at gamma=10.
base.nuEffForRe = NaN;
base.targetRe = 100;
base.targetStrouhal = 0.18;

base.projectionStrength = 1.0;
base.projectedStrength = 1.0;
base.massFluxProjectionMode = 'relax_to_uniform_lowk';
base.massFluxProjectionStrength = 1.0;
base.massFluxDensityRelaxationBeta = 0.002;
base.massFluxApplyAfterVelocityProjection = true;
base.massFluxFinalVelocityProjectionCleanup = true;
base.massFluxFinalVelocityProjectionStrength = 0.5;
base.massFluxTargetFilter = 'lowpass_fft';
base.massFluxLowKMaxIndex = 2;

base.storeDensityMaps = true;
base.storeVorticityMaps = false;
base.densityBandFraction = 0.20;
base.densityExcludeWallCells = 0;
base.lowKMaxIndex = 2;
base.makeFigures = false;

%% === Cases ===
% D = 2R, nominal U is initialMeanVelocityX.  The true U used in St is the
% measured mean Ux from the run.  Start with one conservative controlled case;
% add extra cases only after this baseline remains stable for the full run.
caseName = strings(1,1);
U0 = zeros(1,1);
bodyForceX = zeros(1,1);
R = zeros(1,1);

caseName(1) = "controlled_U001_R010_specular_soft";
U0(1) = 0.02;
bodyForceX(1) = 0.0;
R(1) = 0.10 * base.Ly;

nCases = numel(caseName);
rows = cell(nCases, 1);
outList = cell(nCases, 1);

for ic = 1:nCases
    fprintf('\n--- Case %d/%d: %s ---\n', ic, nCases, caseName(ic));
    params = base;
    params.initialMeanVelocityX = U0(ic);
    params.targetMeanVelocityX = U0(ic);
    params.bodyForceX = bodyForceX(ic);
    params.cylinderRadius = R(ic);
    params.caseName = char(caseName(ic));
    params.makeFigures=1;

    try
        out = run_projection_cylinder_demo(params);
        outList{ic} = out;
        s = out.summary;
        rows{ic} = table(caseName(ic), params.nSteps, s.actualLastStep, s.stoppedEarly, params.seed, params.Nx, params.Ny, params.gamma, ...
            params.initialMeanVelocityX, params.bodyForceX, params.cylinderRadius, ...
            s.meanMeanUx, s.finalMeanUx, s.reEffMeanUx, s.reEffFinalUx, ...
            s.meanStdN, s.meanOutBand, s.meanLowKEnergy, s.timeAvgRelRms, ...
            s.meanDensityTransportProjectedRms, s.meanRmsDivParticleAfter, ...
            s.meanEnstrophy, s.finalEnstrophy, s.wakeUyUpperRms, s.wakeOmegaAntiSymRms, ...
            s.sheddingFrequencyWakeUyUpper, s.strouhalWakeUyUpper, s.sheddingSNRWakeUyUpper, ...
            s.sheddingFrequencyWakeOmegaAntiSym, s.strouhalWakeOmegaAntiSym, s.sheddingSNROmegaAntiSym, ...
            'VariableNames', {'caseName','nSteps','actualLastStep','stoppedEarly','seed','Nx','Ny','gamma', ...
            'initialMeanVelocityX','bodyForceX','cylinderRadius', ...
            'meanUx','finalUx','reEffMeanUx','reEffFinalUx', ...
            'meanStdN','meanOutBand','meanLowK','timeAvgRelRms', ...
            'meanRhoTransport','meanDivUAfter','meanEnstrophy','finalEnstrophy', ...
            'wakeUyUpperRms','wakeOmegaAntiSymRms', ...
            'shedFreqUyUpper','strouhalUyUpper','shedSNRUyUpper', ...
            'shedFreqOmegaAntiSym','strouhalOmegaAntiSym','shedSNROmegaAntiSym'});
    catch ME
        warning('Case %s failed: %s', caseName(ic), ME.message);
        rows{ic} = failure_row(caseName(ic), base, U0(ic), bodyForceX(ic), R(ic));
    end

    partial = vertcat(rows{1:ic});
    writetable(partial, fullfile(outputDir, 'q9_vonkarman_re_st_sweep_summary_partial.csv'));
    save(fullfile(outputDir, 'q9_vonkarman_re_st_sweep_partial.mat'), 'base', 'caseName', 'U0', 'bodyForceX', 'R', 'rows', 'outList', '-v7.3');
end

summary = vertcat(rows{:});
csvFile = fullfile(outputDir, 'q9_vonkarman_re_st_sweep_summary.csv');
matFile = fullfile(outputDir, 'q9_vonkarman_re_st_sweep.mat');
txtFile = fullfile(outputDir, 'q9_vonkarman_re_st_sweep_summary.txt');
writetable(summary, csvFile);
save(matFile, 'base', 'caseName', 'U0', 'bodyForceX', 'R', 'summary', 'outList', '-v7.3');

fid = fopen(txtFile, 'w');
fprintf(fid, '=== Q9 von Karman Re/St exploratory sweep ===\n');
fprintf(fid, 'nCases                            : %d\n', nCases);
fprintf(fid, 'nSteps/sampleEvery                : %d / %d\n', base.nSteps, base.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %g\n', base.Nx, base.Ny, base.gamma);
fprintf(fid, 'targetRe / targetSt / nuEffForRe  : %.12g / %.12g / %.12g\n', base.targetRe, base.targetStrouhal, base.nuEffForRe);
fprintf(fid, 'wake probe x/D, y/D offset        : %.12g / %.12g\n', base.wakeProbeDistanceDiameters, base.wakeProbeYOffsetDiameters);
fprintf(fid, 'mean-flow control mode/target/tau  : %s / %.12g / %.12g\n', char(base.meanFlowControlMode), base.targetMeanVelocityX, base.meanFlowRelaxationTau);
fprintf(fid, 'shedding St search band            : %.12g / %.12g\n', base.sheddingStrouhalMin, base.sheddingStrouhalMax);
for i = 1:height(summary)
    fprintf(fid, '\ncase: %s\n', char(summary.caseName(i)));
    fprintf(fid, 'U0, bodyForceX, R                 : %.12g / %.12g / %.12g\n', summary.initialMeanVelocityX(i), summary.bodyForceX(i), summary.cylinderRadius(i));
    fprintf(fid, 'actualLastStep / stoppedEarly     : %.12g / %d\n', summary.actualLastStep(i), logical(summary.stoppedEarly(i)));
    fprintf(fid, 'mean/final Ux                     : %.12g / %.12g\n', summary.meanUx(i), summary.finalUx(i));
    fprintf(fid, 'estimated Re mean/final           : %.12g / %.12g\n', summary.reEffMeanUx(i), summary.reEffFinalUx(i));
    fprintf(fid, 'mean low-k / divU                 : %.12g / %.12g\n', summary.meanLowK(i), summary.meanDivUAfter(i));
    fprintf(fid, 'mean/final enstrophy              : %.12g / %.12g\n', summary.meanEnstrophy(i), summary.finalEnstrophy(i));
    fprintf(fid, 'wake Uy upper RMS                 : %.12g\n', summary.wakeUyUpperRms(i));
    fprintf(fid, 'omega antisym RMS                 : %.12g\n', summary.wakeOmegaAntiSymRms(i));
    fprintf(fid, 'shedding f, St upper Uy, SNR      : %.12g / %.12g / %.12g\n', summary.shedFreqUyUpper(i), summary.strouhalUyUpper(i), summary.shedSNRUyUpper(i));
    fprintf(fid, 'shedding f, St omega antisym, SNR : %.12g / %.12g / %.12g\n', summary.shedFreqOmegaAntiSym(i), summary.strouhalOmegaAntiSym(i), summary.shedSNROmegaAntiSym(i));
end
fprintf(fid, '\n--- Files ---\n%s\n%s\n%s\n', matFile, csvFile, txtFile);
fclose(fid);

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);

diary off

function T = failure_row(caseName, base, U0, bodyForceX, R)
T = table(caseName, base.nSteps, nan, true, base.seed, base.Nx, base.Ny, base.gamma, ...
    U0, bodyForceX, R, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, ...
    nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, ...
    'VariableNames', {'caseName','nSteps','actualLastStep','stoppedEarly','seed','Nx','Ny','gamma', ...
    'initialMeanVelocityX','bodyForceX','cylinderRadius', ...
    'meanUx','finalUx','reEffMeanUx','reEffFinalUx', ...
    'meanStdN','meanOutBand','meanLowK','timeAvgRelRms', ...
    'meanRhoTransport','meanDivUAfter','meanEnstrophy','finalEnstrophy', ...
    'wakeUyUpperRms','wakeOmegaAntiSymRms', ...
    'shedFreqUyUpper','strouhalUyUpper','shedSNRUyUpper', ...
    'shedFreqOmegaAntiSym','strouhalOmegaAntiSym','shedSNROmegaAntiSym'});
end
