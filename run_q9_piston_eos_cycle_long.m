%% run_q9_piston_eos_cycle_long.m
% Long quasi-static piston cycle for Q9 effective-EOS diagnostics.
%
% This wrapper uses the improved diagnostics from run_q9_piston_eos_cycle.m
% with a longer compression/hold/decompression/hold protocol.
%
% Main outputs added by the v2 diagnostics:
%   - baseline-subtracted pressures relative to the final relaxed hold,
%   - plateau EOS estimate between compressed and relaxed holds,
%   - compression/decompression hysteresis diagnostics,
%   - relative-pressure EOS figures.

pistonEosCyclePreset = 'long';
pistonEosOutputPrefix = 'q9_piston_eos_cycle_long';

% Optional local overrides can be uncommented here, or set in the workspace
% before running this wrapper.
% pistonEosUserParams = struct();
% pistonEosUserParams.pistonCompressionTarget = 0.05;
% pistonEosUserParams.massFluxDensityRelaxationBeta = 0.002;
% pistonEosUserParams.massFluxLowKMaxIndex = 2;
% pistonEosUserParams.massFluxFinalVelocityProjectionStrength = 0.5;

run_q9_piston_eos_cycle
