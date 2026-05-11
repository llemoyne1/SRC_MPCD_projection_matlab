%% run_q9_lowk_mass_flux_3000_check.m
clear functions
close all
clc

tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q9_lowk_mass_flux_3000_check_' tag];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
params = struct();

params.Nx = 32;
params.Ny = 16;
params.gamma = 20;

params.nSteps = 3000;
params.sampleEvery = 100;
params.progressEvery = 500;
params.maxWallClockSeconds = Inf;

params.bodyForceX = 0.02;
params.wallModeY = 'thermalize';
params.thermostatAfterProjection = true;

params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';

%% Q6 velocity projection active
params.projectionStrength = 1.0;
params.projectedStrength = 1.0;

%% Q9 low-k mass-flux correction
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;

params.massFluxApplyAfterVelocityProjection = true;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;

%% Diagnostics
params.storeDensityMaps = true;
params.densityBandFraction = 0.20;
params.densityExcludeWallCells = 0;
params.lowKMaxIndex = 2;

params.makeFigures = false;
params.saveFigures = false;
params.writeSummary = true;
params.outputDir = outputDir;

diary(fullfile(outputDir, 'console_log.txt'));
diary on

fprintf('\n=== Q9 low-k mass-flux 3000-step check ===\n');

cmp = run_compare_density_projection_lowk_mass_flux_poiseuille(params);

diary off