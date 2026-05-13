%% run_q9_piston_eos_staircase_resolved3.m
% Three-plateau, better-resolved EOS staircase run.
%
% This wrapper uses fewer density levels but longer plateaux, so that the
% top-wall mechanical pressure has a more meaningful standard error.

pistonEosOutputPrefix = 'q9_piston_eos_staircase_resolved3';
pistonEosIncludeDecompression = false;
pistonEosRhoRatioList = [1.00 1.05 1.10];

pistonEosUserParams = struct();
pistonEosUserParams.stepsPerPlateau = 6000;
pistonEosUserParams.sampleEvery = 50;
pistonEosUserParams.plateauDiscardFraction = 0.50;
pistonEosUserParams.useCumulativeWallPressureForFit = true;
pistonEosUserParams.makeFigures = true;

run_q9_piston_eos_staircase
