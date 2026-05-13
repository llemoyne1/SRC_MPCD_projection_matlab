function [stateOut, diag] = mpcd_step_projection_periodic_q9(state, params)
%MPCD_STEP_PROJECTION_PERIODIC_Q9 Classic periodic MPCD step plus Q6/Q9.
%
% Q6 is obtained with massFluxProjectionMode='off'.  Q9 is obtained with
% massFluxProjectionMode='relax_to_uniform_lowk'.

[stateClassic, classicDiag] = mpcd_step_classic_periodic(state, params);
[stateOut, projDiag] = mpcd_apply_q9_projection_periodic(stateClassic, params);

diag = projDiag;
diag.classic = classicDiag;
end
