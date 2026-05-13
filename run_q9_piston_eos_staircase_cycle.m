%% run_q9_piston_eos_staircase_cycle.m
% Up/down staircase wrapper for Q9 effective-EOS diagnostics.
%
% This uses the same quasi-static plateau protocol as
% run_q9_piston_eos_staircase.m, but adds decompression plateaux to diagnose
% hysteresis at much lower cost than a slow continuously moving piston.

pistonEosIncludeDecompression = true;
pistonEosOutputPrefix = 'q9_piston_eos_staircase_cycle';

% Optional default can be overridden before running this wrapper by defining
% pistonEosRhoRatioList or pistonEosUserParams in the caller workspace.
run_q9_piston_eos_staircase
