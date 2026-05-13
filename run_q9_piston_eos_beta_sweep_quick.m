%% run_q9_piston_eos_beta_sweep_quick.m
% Cheaper first-pass beta sweep. Use this to check monotonicity before the
% resolved run.

pistonEosBetaList = [0, 5e-4, 1e-3, 2e-3, 4e-3];
pistonEosSweepUserParams = struct();
pistonEosSweepUserParams.rhoRatioList = [1.00 1.05 1.10];
pistonEosSweepUserParams.stepsPerPlateau = 3000;
pistonEosSweepUserParams.sampleEvery = 50;
pistonEosSweepUserParams.plateauDiscardFraction = 0.50;

run_q9_piston_eos_beta_sweep
