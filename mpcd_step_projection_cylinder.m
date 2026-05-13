function [stateOut, diag] = mpcd_step_projection_cylinder(state, params)
%MPCD_STEP_PROJECTION_CYLINDER Cylinder SRC/MPCD step followed by Q6/Q9.
%
%   [stateOut, diag] = mpcd_step_projection_cylinder(state, params)

[stateClassic, classicDiag] = mpcd_step_classic_cylinder(state, params);
[stateOut, projDiag] = mpcd_apply_q9_projection_channel(stateClassic, params);

diag = projDiag;
diag.classicDiag = classicDiag;
diag.cylinderInfo = classicDiag.cylinderInfo;
diag.wallInfoY = classicDiag.wallInfoY;
diag.wallInfo = classicDiag.wallInfo;
diag.cylinderHits = classicDiag.cylinderInfo.nCylinderHits;
diag.cylinderDPx = classicDiag.cylinderInfo.dPxCylinder;
diag.cylinderDPy = classicDiag.cylinderInfo.dPyCylinder;
diag.classicMeanVx = classicDiag.meanVx;
diag.classicMeanVy = classicDiag.meanVy;
end
