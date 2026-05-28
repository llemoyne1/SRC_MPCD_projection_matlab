function [vOut, info] = projection_apply_taylor_green_forcing(x, v, params)
%PROJECTION_APPLY_TAYLOR_GREEN_FORCING Apply divergence-free TG body-force kick.
%
% The acceleration has the same mode shape as the Taylor-Green velocity:
%   f_x = F sin(kx x) cos(ky y)
%   f_y = -F kx/ky cos(kx x) sin(ky y)
%
% The kick can be zero-meaned over particles to prevent sampling noise from
% introducing a net force in a nominally zero-mean periodic forcing.

enabled = logical(get_param(params, 'taylorGreenForceEnable', false));
F = get_param(params, 'taylorGreenForceAmplitude', 0.0);
dt = get_param(params, 'dt', 1.0);
info = struct('enabled', enabled, 'amplitude', F, 'meanKickRaw', [0 0], ...
    'meanKickApplied', [0 0], 'rmsKickApplied', 0, 'zeroMeanKick', false);

if ~enabled || F == 0
    vOut = v;
    return;
end

[fx, fy] = projection_taylor_green_mode_at_points(x, params, F);
dv = dt * [fx, fy];
info.meanKickRaw = mean(dv, 1, 'omitnan');
zeroMeanKick = logical(get_param(params, 'taylorGreenForceZeroMeanKick', true));
info.zeroMeanKick = zeroMeanKick;
if zeroMeanKick
    dv = dv - mean(dv, 1, 'omitnan');
end
info.meanKickApplied = mean(dv, 1, 'omitnan');
info.rmsKickApplied = sqrt(mean(sum(dv.^2, 2), 'omitnan'));
vOut = v + dv;
end

function v = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end
