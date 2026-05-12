%% run_q9_reference_lowk_mass_flux_30000.m
clear functions
close all
clc

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q9_reference_lowk_mass_flux_30000_' tag];

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on

fprintf('\n=== Q9 reference: velocity projection + low-k mass-flux relaxation ===\n');

%% === Q6 projection-only reference, 30000 steps, seed=11 ===
% Reference from validated Q6 run:
% projectionStrength = 1
% massFluxProjectionMode = off
% initialPopulationMode = exact_per_cell
% seed = 11
% nSteps = 30000

refQ6 = struct();

refQ6.meanStdN = 3.541271774591202;
refQ6.finalStdN = 3.36379;
refQ6.meanOutBandFraction = 0.202547;
refQ6.meanLowKEnergy = 0.000174985711550722;
refQ6.meanRhoTransport = 0.009443045495482883;
refQ6.timeAvgRelRms = 0.0232799;
refQ6.nuEff = 0.02011848517676771;
refQ6.R2 = 0.9641825015280132;

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

%% === Q9 method ===
% Q9 = Q6 velocity projection + low-k mass-flux density relaxation.

params.projectionStrength = 1.0;
params.projectedStrength = 1.0;

params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;

params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;

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

%% === Run Q9 only ===
tic

outQ9 = run_projection_poiseuille_long_demo(params);

elapsedScript = toc;

%% === Density metrics ===
metricsQ9 = projection_density_homogeneity_metrics( ...
    outQ9.NMaps, ...
    params, ...
    'sampleTimes', outQ9.sampleTimes, ...
    'targetGamma', params.gamma, ...
    'bandFraction', params.densityBandFraction, ...
    'excludeWallCells', params.densityExcludeWallCells, ...
    'lowKMaxIndex', params.lowKMaxIndex ...
);

%% === Summary table ===
summary = table();

summary.method = string('Q9_lowk_mass_flux');
summary.seed = params.seed;
summary.nSteps = params.nSteps;
summary.Nx = params.Nx;
summary.Ny = params.Ny;
summary.gamma = params.gamma;

summary.projectionStrength = params.projectionStrength;
summary.massFluxMode = string(params.massFluxProjectionMode);
summary.massFluxStrength = params.massFluxProjectionStrength;
summary.beta = params.massFluxDensityRelaxationBeta;
summary.massFluxApplyAfterVelocityProjection = params.massFluxApplyAfterVelocityProjection;
summary.massFluxTargetFilter = string(params.massFluxTargetFilter);
summary.massFluxLowKMaxIndex = params.massFluxLowKMaxIndex;

%% === Density summary ===
summary.meanStdN = metricsQ9.summary.meanStdN;
summary.ratioStd_vs_Q6Projection = ...
    metricsQ9.summary.meanStdN / refQ6.meanStdN;

summary.finalStdN = metricsQ9.summary.finalStdN;
summary.ratioFinalStd_vs_Q6Projection = ...
    metricsQ9.summary.finalStdN / refQ6.finalStdN;

summary.meanOutBand = metricsQ9.summary.meanOutBandFraction;
summary.ratioOutBand_vs_Q6Projection = ...
    metricsQ9.summary.meanOutBandFraction / refQ6.meanOutBandFraction;

summary.meanLowK = metricsQ9.summary.meanLowKEnergy;
summary.ratioLowK_vs_Q6Projection = ...
    metricsQ9.summary.meanLowKEnergy / refQ6.meanLowKEnergy;

summary.timeAvgRelRms = metricsQ9.summary.timeAvgRelRms;
summary.ratioTimeAvgRelRms_vs_Q6Projection = ...
    metricsQ9.summary.timeAvgRelRms / refQ6.timeAvgRelRms;

summary.meanRhoTransport = outQ9.summary.meanDensityTransportProjectedRms;
summary.ratioRhoTransport_vs_Q6Projection = ...
    outQ9.summary.meanDensityTransportProjectedRms / refQ6.meanRhoTransport;

%% === Hydrodynamics ===
summary.nuEff = outQ9.viscosity.nuEff;
summary.ratioNuEff_vs_Q6Projection = ...
    outQ9.viscosity.nuEff / refQ6.nuEff;

summary.R2 = outQ9.viscosity.R2;
summary.signalToNoise = outQ9.viscosity.signalToNoise;
summary.centerMinusWall = outQ9.viscosity.centerMinusWall;

%% === Mass-flux diagnostics ===
if isfield(outQ9.summary, 'meanMassFluxDivBefore')
    summary.meanMassFluxDivBefore = outQ9.summary.meanMassFluxDivBefore;
else
    summary.meanMassFluxDivBefore = NaN;
end

if isfield(outQ9.summary, 'meanMassFluxDivTarget')
    summary.meanMassFluxDivTarget = outQ9.summary.meanMassFluxDivTarget;
else
    summary.meanMassFluxDivTarget = NaN;
end

if isfield(outQ9.summary, 'meanMassFluxDivProjectedAfter')
    summary.meanMassFluxDivAfter = outQ9.summary.meanMassFluxDivProjectedAfter;
else
    summary.meanMassFluxDivAfter = NaN;
end

if isfield(outQ9.summary, 'meanMassFluxDivResidual')
    summary.meanMassFluxResidual = outQ9.summary.meanMassFluxDivResidual;
else
    summary.meanMassFluxResidual = NaN;
end

summary.massFluxResidualRelative = summary.meanMassFluxResidual / ...
    max([summary.meanMassFluxDivBefore, summary.meanMassFluxDivTarget, 1]);

summary.elapsedRun = outQ9.elapsedWallClock;
summary.elapsedScript = elapsedScript;

disp(summary)

%% === Save data ===
matFile = fullfile(outputDir, 'q9_reference_lowk_mass_flux_30000.mat');
csvFile = fullfile(outputDir, 'q9_reference_lowk_mass_flux_30000_summary.csv');
txtFile = fullfile(outputDir, 'q9_reference_lowk_mass_flux_30000_summary.txt');

save(matFile, 'params', 'outQ9', 'metricsQ9', 'summary', 'refQ6', '-v7.3');
writetable(summary, csvFile);

%% === Write compact text summary ===
fid = fopen(txtFile, 'w');

fprintf(fid, '=== Q9 reference: velocity projection + low-k mass-flux relaxation ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'steps/sampleEvery                 : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'initialPopulationMode             : %s\n', params.initialPopulationMode);

fprintf(fid, '\n--- Method parameters ---\n');
fprintf(fid, 'projectionStrength                : %.6g\n', params.projectionStrength);
fprintf(fid, 'massFluxProjectionMode            : %s\n', params.massFluxProjectionMode);
fprintf(fid, 'massFluxProjectionStrength        : %.6g\n', params.massFluxProjectionStrength);
fprintf(fid, 'massFluxDensityRelaxationBeta     : %.6g\n', params.massFluxDensityRelaxationBeta);
fprintf(fid, 'massFluxApplyAfterVelocityProjection : %d\n', params.massFluxApplyAfterVelocityProjection);
fprintf(fid, 'massFluxTargetFilter              : %s\n', params.massFluxTargetFilter);
fprintf(fid, 'massFluxLowKMaxIndex              : %d\n', params.massFluxLowKMaxIndex);

fprintf(fid, '\n--- Density metrics ---\n');
fprintf(fid, 'mean std(N)                       : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanStdN, summary.ratioStd_vs_Q6Projection);
fprintf(fid, 'final std(N)                      : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.finalStdN, summary.ratioFinalStd_vs_Q6Projection);
fprintf(fid, 'mean out-band                     : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanOutBand, summary.ratioOutBand_vs_Q6Projection);
fprintf(fid, 'mean low-k energy                 : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanLowK, summary.ratioLowK_vs_Q6Projection);
fprintf(fid, 'time-avg rel RMS                  : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.timeAvgRelRms, summary.ratioTimeAvgRelRms_vs_Q6Projection);
fprintf(fid, 'mean rho transport                : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.meanRhoTransport, summary.ratioRhoTransport_vs_Q6Projection);

fprintf(fid, '\n--- Hydrodynamics ---\n');
fprintf(fid, 'nu_eff                            : %.12g   ratio vs Q6 P = %.6g\n', ...
    summary.nuEff, summary.ratioNuEff_vs_Q6Projection);
fprintf(fid, 'R2                                : %.12g\n', summary.R2);
fprintf(fid, 'signal/noise                      : %.12g\n', summary.signalToNoise);
fprintf(fid, 'center minus wall U               : %.12g\n', summary.centerMinusWall);

fprintf(fid, '\n--- Mass-flux diagnostics ---\n');
fprintf(fid, 'mean div(Nu) before               : %.12g\n', summary.meanMassFluxDivBefore);
fprintf(fid, 'mean div(Nu) target               : %.12g\n', summary.meanMassFluxDivTarget);
fprintf(fid, 'mean div(Nu) after                : %.12g\n', summary.meanMassFluxDivAfter);
fprintf(fid, 'mean residual                     : %.12g\n', summary.meanMassFluxResidual);
fprintf(fid, 'relative residual                 : %.12g\n', summary.massFluxResidualRelative);

fprintf(fid, '\n--- Timing ---\n');
fprintf(fid, 'elapsed run [s]                   : %.3f\n', summary.elapsedRun);
fprintf(fid, 'elapsed script [s]                : %.3f\n', summary.elapsedScript);

fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', txtFile);

fclose(fid);

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);

diary off