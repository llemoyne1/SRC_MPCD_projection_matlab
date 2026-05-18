function [stateOut, diag] = mpcd_step_projection_poiseuille_q9(state, params)
%MPCD_STEP_PROJECTION_POISEUILLE_Q9 Poiseuille SRC/MPCD step followed by Q6/Q9.
%
% This wrapper reuses the corrected channel projection block from the recent
% semi-parallel snapshot.  To make comparisons with the classic method fair,
% the classic substep is run without a thermostat and the common thermostat is
% applied only after the Q6/Q9 projection through mpcd_apply_q9_projection_channel.

pClassic = params;
pClassic.thermostatAfterStep = false;
[stateClassic, classicDiag] = mpcd_step_classic_poiseuille(state, pClassic);

pProj = params;
pProj.thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', ...
    get_param(params, 'thermostatAfterStep', false)));
[stateOut, projDiag] = mpcd_apply_q9_projection_channel(stateClassic, pProj);

diag = projDiag;
diag.classic = classicDiag;
diag.classicMeanVx = classicDiag.meanVx;
diag.classicMeanVy = classicDiag.meanVy;
diag.wallInfo = classicDiag.wallInfo;
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
