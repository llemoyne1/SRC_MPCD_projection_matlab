function [x, v, info] = mpcd_apply_wall_bc_y(x, v, params)
%MPCD_APPLY_WALL_BC_Y Apply simple y-wall conditions for Poiseuille prototype.
%
%   [x, v, info] = mpcd_apply_wall_bc_y(x, v, params)
%
% x is periodic elsewhere; this function only handles y=0 and y=Ly walls.
%
% Optional params fields:
%   wallModeY     'bounceback', 'specular', or 'thermalize', default 'bounceback'
%   Ubottom       default 0
%   Utop          default 0
%   kBT           default 1
%   wallSigma     default sqrt(kBT)
%
% bounceback reverses the tangential component relative to the wall velocity:
%   vx <- 2*Uwall - vx, vy reflected toward the domain.

Ly = params.Ly;
mode = lower(char(string(get_param(params, 'wallModeY', 'bounceback'))));
Ubottom = get_param(params, 'Ubottom', 0.0);
Utop = get_param(params, 'Utop', 0.0);
wallSigma = get_param(params, 'wallSigma', sqrt(max(get_param(params, 'kBT', 1.0), 0)));

info = struct('nBot', 0, 'nTop', 0, 'dEwall', 0, ...
              'dPxBot', 0, 'dPxTop', 0, 'dPyBot', 0, 'dPyTop', 0);

% A few passes are enough for the small displacements expected here.
for pass = 1:6
    hitB = x(:, 2) < 0;
    if any(hitB)
        ids = find(hitB);
        vOld = v(ids, :);
        x(ids, 2) = -x(ids, 2);
        switch mode
            case 'thermalize'
                v(ids, 1) = Ubottom + wallSigma * randn(numel(ids), 1);
                v(ids, 2) = abs(wallSigma * randn(numel(ids), 1));
            case 'bounceback'
                v(ids, 1) = 2*Ubottom - v(ids, 1);
                v(ids, 2) = abs(v(ids, 2));
            case 'specular'
                v(ids, 2) = abs(v(ids, 2));
            otherwise
                error('Unknown wallModeY: %s', mode);
        end
        dv = v(ids, :) - vOld;
        info.nBot = info.nBot + numel(ids);
        info.dEwall = info.dEwall + 0.5 * sum(sum(v(ids, :).^2 - vOld.^2, 2));
        info.dPxBot = info.dPxBot + sum(dv(:, 1));
        info.dPyBot = info.dPyBot + sum(abs(dv(:, 2)));
    end

    hitT = x(:, 2) > Ly;
    if any(hitT)
        ids = find(hitT);
        vOld = v(ids, :);
        x(ids, 2) = 2*Ly - x(ids, 2);
        switch mode
            case 'thermalize'
                v(ids, 1) = Utop + wallSigma * randn(numel(ids), 1);
                v(ids, 2) = -abs(wallSigma * randn(numel(ids), 1));
            case 'bounceback'
                v(ids, 1) = 2*Utop - v(ids, 1);
                v(ids, 2) = -abs(v(ids, 2));
            case 'specular'
                v(ids, 2) = -abs(v(ids, 2));
            otherwise
                error('Unknown wallModeY: %s', mode);
        end
        dv = v(ids, :) - vOld;
        info.nTop = info.nTop + numel(ids);
        info.dEwall = info.dEwall + 0.5 * sum(sum(v(ids, :).^2 - vOld.^2, 2));
        info.dPxTop = info.dPxTop + sum(dv(:, 1));
        info.dPyTop = info.dPyTop + sum(abs(dv(:, 2)));
    end

    if ~any(hitB) && ~any(hitT)
        break;
    end
end

x(:, 2) = min(max(x(:, 2), 0), Ly);
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
