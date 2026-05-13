function [stateOut, diag] = mpcd_step_projection_step_channel(state, params)
%MPCD_STEP_PROJECTION_STEP_CHANNEL Rectangular-step SRC/MPCD step followed by Q6/Q9.

[stateClassic, classicDiag] = mpcd_step_classic_step_channel(state, params);
[stateOut, projDiag] = mpcd_apply_q9_projection_channel(stateClassic, params);

diag = projDiag;
diag.classicDiag = classicDiag;
diag.stepInfo = classicDiag.stepInfo;
diag.wallInfoY = classicDiag.wallInfoY;
diag.wallInfo = classicDiag.wallInfo;
diag.stepHits = classicDiag.stepInfo.nStepHits;
diag.stepDPx = classicDiag.stepInfo.dPxStep;
diag.stepDPy = classicDiag.stepInfo.dPyStep;
diag.classicMeanVx = classicDiag.meanVx;
diag.classicMeanVy = classicDiag.meanVy;
end
