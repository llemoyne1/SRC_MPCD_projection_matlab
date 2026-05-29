function [stateOut, diag] = resamp_step_projection_poiseuille_weighted(state, params)
%RESAMP_STEP_PROJECTION_POISEUILLE_WEIGHTED Weighted channel step plus Q6.

pClassic=params; pClassic.thermostatAfterStep=false;
[stateClassic,classicDiag]=resamp_step_classic_poiseuille_weighted(state,pClassic);
[stateOut,projDiag]=resamp_apply_q6_channel_weighted(stateClassic,params);
diag=projDiag; diag.classic=classicDiag; diag.wallInfo=classicDiag.wallInfo; diag.classicMeanVx=classicDiag.meanVxWeighted; diag.classicMeanVy=classicDiag.meanVyWeighted;
end
