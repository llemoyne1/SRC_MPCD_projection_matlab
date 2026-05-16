function [stateOut, diag] = mpcd_step_projection_periodic_q9_forced(state, params)
%MPCD_STEP_PROJECTION_PERIODIC_Q9_FORCED Periodic MPCD step plus Q6/Q9 projection.
[stateClassic, classicDiag] = mpcd_step_classic_periodic_forced(state, params);
[stateOut, projDiag] = mpcd_apply_q9_projection_periodic(stateClassic, params);
diag = projDiag;
diag.classic = classicDiag;
end
