function [stateOut, diag] = resamp_step_projection_periodic_weighted(state, params)
%RESAMP_STEP_PROJECTION_PERIODIC_WEIGHTED Weighted periodic MPCD step plus Q6 projection.
%
% This is deliberately Q6-only.  Q9 is kept out of the first weighted-mass
% infrastructure patch because the new direction is to first establish
% variable particle masses and weighted hydrodynamic deposits.

[stateClassic, classicDiag] = resamp_step_classic_periodic_weighted(state, params);
[stateOut, projDiag] = resamp_apply_q6_periodic_weighted(stateClassic, params);
diag = projDiag;
diag.classic = classicDiag;
end
