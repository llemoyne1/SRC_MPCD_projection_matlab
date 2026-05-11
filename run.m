%% run_q8_mass_flux_relax_beta0002_30000.m
clear functions
close all
clc

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q8_mass_flux_relax_beta0002_30000_' tag];

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on

fprintf('\n=== Q8 mass-flux relax-to-uniform beta=0.002, 30000 steps ===\n');

%% === Q6 projection-only reference, 30000 steps, seed=11 ===
refQ6 = struct();

refQ6.meanStdN = 3.54127;
refQ6.finalStdN = 3.36379;
refQ6.meanOutBandFraction = 0.202547;
refQ6.meanLowKEnergy = 0.000174986;
refQ6.meanRhoTransport = 0.00944305;
refQ6.timeAvgRelRms = 0.0232799;
refQ6.nuEff = 2.011848517677e-02;
refQ6.R2 = 0.964183;

%% === Parameters ===
params = struct();

params.Nx = 32;
params.Ny = 16;
params.gamma = 20;

params.nSteps = 30000;
params.sampleEvery = 100;
params.progressEvery = 1000;
params.maxWallClockSeconds = Inf;

params.bodyForceX = 0.02;
params.wallModeY = 'thermalize';
params.thermostatAfterProjection = true;

params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';

%% === Pure Q8 mass-flux projection ===
% Disable velocity projection div(u)=0.
params.projectionStrength = 0.0;
params.projectedStrength = 0.0;

% Enable Q8 relaxed mass-flux projection.
params.massFluxProjectionMode = 'relax_to_uniform';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;

%% === Density diagnostics ===
params.storeDensityMaps = true;
params.densityBandFraction = 0.20;
params.densityExcludeWallCells = 0;
params.lowKMaxIndex = 2;

%% === Figures / files ===
params.makeFigures = true;
params.saveFigures = true;
params.writeSummary = true;
params.outputDir = outputDir;

%% === Run ===
tic

outM = run_projection_poiseuille_long_demo(params);

elapsedScript = toc;

%% === Density metrics ===
metricsM = projection_density_homogeneity_metrics( ...
    outM.NMaps, ...
    params, ...
    'sampleTimes', outM.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex ...
);

%% === Summary table ===
summary = table();

summary.seed = params.seed;
summary.nSteps = params.nSteps;
summary.Nx = params.Nx;
summary.Ny = params.Ny;
summary.gamma = params.gamma;

summary.massFluxMode = string(params.massFluxProjectionMode);
summary.massFluxStrength = params.massFluxProjectionStrength;
summary.beta = params.massFluxDensityRelaxationBeta;

%% === Density summary ===
summary.meanStdN = metricsM.summary.meanStdN;
summary.ratioStd_vs_Q6Projection = ...
    metricsM.summary.meanStdN / refQ6.meanStdN;

summary.finalStdN = metricsM.summary.finalStdN;
summary.ratioFinalStd_vs_Q6Projection = ...
    metricsM.summary.finalStdN / refQ6.finalStdN;

summary.meanOutBand = metricsM.summary.meanOutBandFraction;
summary.ratioOutBand_vs_Q6Projection = ...
    metricsM.summary.meanOutBandFraction / refQ6.meanOutBandFraction;

summary.meanLowK = metricsM.summary.meanLowKEnergy;
summary.ratioLowK_vs_Q6Projection = ...
    metricsM.summary.meanLowKEnergy / refQ6.meanLowKEnergy;

summary.timeAvgRelRms = metricsM.summary.timeAvgRelRms;
summary.ratioTimeAvgRelRms_vs_Q6Projection = ...
    metricsM.summary.timeAvgRelRms / refQ6.timeAvgRelRms;

summary.meanRhoTransport = outM.summary.meanDensityTransportProjectedRms;
summary.ratioRhoTransport_vs_Q6Projection = ...
    outM.summary.meanDensityTransportProjectedRms / refQ6.meanRhoTransport;

%% === Hydrodynamics ===
summary.nuEff = outM.viscosity.nuEff;
summary.ratioNuEff_vs_Q6Projection = ...
    outM.viscosity.nuEff / refQ6.nuEff;

summary.R2 = outM.viscosity.R2;
summary.signalToNoise = outM.viscosity.signalToNoise;
summary.centerMinusWall = outM.viscosity.centerMinusWall;

%% === Mass-flux diagnostics ===
if isfield(outM.summary, 'meanMassFluxDivBefore')
    summary.meanMassFluxDivBefore = outM.summary.meanMassFluxDivBefore;
else
    summary.meanMassFluxDivBefore = NaN;
end

if isfield(outM.summary, 'meanMassFluxDivTarget')
    summary.meanMassFluxDivTarget = outM.summary.meanMassFluxDivTarget;
else
    summary.meanMassFluxDivTarget = NaN;
end

if isfield(outM.summary, 'meanMassFluxDivProjectedAfter')
    summary.meanMassFluxDivAfter = outM.summary.meanMassFluxDivProjectedAfter;
else
    summary.meanMassFluxDivAfter = NaN;
end

if isfield(outM.summary, 'meanMassFluxDivResidual')
    summary.meanMassFluxResidual = outM.summary.meanMassFluxDivResidual;
else
    summary.meanMassFluxResidual = NaN;
end

summary.massFluxResidualRelative = summary.meanMassFluxResidual / ...
    max([summary.meanMassFluxDivBefore, summary.meanMassFluxDivTarget, 1]);

summary.elapsedRun = outM.elapsedWallClock;
summary.elapsedScript = elapsedScript;

disp(summary)

%% === Save data ===
matFile = fullfile(outputDir, 'q8_mass_flux_relax_beta0002_30000.mat');
csvFile = fullfile(outputDir, 'q8_mass_flux_relax_beta0002_30000_summary.csv');
txtFile = fullfile(outputDir, 'q8_mass_flux_relax_beta0002_30000_summary.txt');

save(matFile, 'params', 'outM', 'metricsM', 'summary', 'refQ6', '-v7.3');
writetable(summary, csvFile);

%% === Write compact text summary ===
fid = fopen(txtFile, 'w');

fprintf(fid, '=== Q8 mass-flux relax-to-uniform beta=0.002, 30000 steps ===\n');
fprintf(fid, 'seed                       : %d\n', params.seed);
fprintf(fid, 'steps/sampleEvery          : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                : %d x %d, %g\n', params.Nx, params.Ny, params.gamma);

fprintf(fid, 'massFlux mode/strength/beta: %s / %.6g / %.6g\n', ...
    params.massFluxProjectionMode, ...
    params.massFluxProjectionStrength, ...
    params.massFluxDensityRelaxationBeta);

fprintf(fid, '\n--- Density metrics ---\n');

fprintf(fid, 'mean std(N)                : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanStdN, ...
    summary.ratioStd_vs_Q6Projection);

fprintf(fid, 'final std(N)               : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.finalStdN, ...
    summary.ratioFinalStd_vs_Q6Projection);

fprintf(fid, 'mean out-band              : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanOutBand, ...
    summary.ratioOutBand_vs_Q6Projection);

fprintf(fid, 'mean low-k energy          : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanLowK, ...
    summary.ratioLowK_vs_Q6Projection);

fprintf(fid, 'time-avg rel RMS           : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.timeAvgRelRms, ...
    summary.ratioTimeAvgRelRms_vs_Q6Projection);

fprintf(fid, 'mean rho transport         : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanRhoTransport, ...
    summary.ratioRhoTransport_vs_Q6Projection);

fprintf(fid, '\n--- Hydrodynamics ---\n');

fprintf(fid, 'nu_eff                     : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.nuEff, ...
    summary.ratioNuEff_vs_Q6Projection);

fprintf(fid, 'R2                         : %.12g\n', summary.R2);
fprintf(fid, 'signal/noise               : %.12g\n', summary.signalToNoise);
fprintf(fid, 'center minus wall U        : %.12g\n', summary.centerMinusWall);

fprintf(fid, '\n--- Mass-flux diagnostics ---\n');

fprintf(fid, 'mean div(Nu) before        : %.12g\n', ...
    summary.meanMassFluxDivBefore);

fprintf(fid, 'mean div(Nu) target        : %.12g\n', ...
    summary.meanMassFluxDivTarget);

fprintf(fid, 'mean div(Nu) after         : %.12g\n', ...
    summary.meanMassFluxDivAfter);

fprintf(fid, 'mean residual              : %.12g\n', ...
    summary.meanMassFluxResidual);

fprintf(fid, 'relative residual          : %.12g\n', ...
    summary.massFluxResidualRelative);

fprintf(fid, '\n--- Timing ---\n');

fprintf(fid, 'elapsed run [s]            : %.3f\n', summary.elapsedRun);
fprintf(fid, 'elapsed script [s]         : %.3f\n', summary.elapsedScript);

fprintf(fid, '\nFiles:\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', txtFile);

fclose(fid);

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);

diary off