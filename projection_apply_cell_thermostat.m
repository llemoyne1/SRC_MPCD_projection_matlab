function [vOut, info] = projection_apply_cell_thermostat(x, v, params, varargin)
%PROJECTION_APPLY_CELL_THERMOSTAT Momentum-preserving local thermostat.
%
%   [vOut, info] = projection_apply_cell_thermostat(x, v, params)
%
% For each populated grid cell, rescales the particle velocity fluctuations
% around the cell mean velocity:
%      v_i = U_cell + s_cell * (v_i - U_cell)
% This preserves the cell momentum exactly and therefore does not change the
% deposited hydrodynamic velocity field when nearest-cell assignment is used.

if nargin < 3
    error('Usage: [vOut, info] = projection_apply_cell_thermostat(x, v, params, ...)');
end

periodicX = true;
periodicY = false;
targetKBT = get_param(params, 'thermostatTargetKBT', get_param(params, 'kBT', 1.0));
strength = get_param(params, 'thermostatStrength', 1.0);
minParticles = get_param(params, 'thermostatMinParticlesPerCell', 2);
maxScale = get_param(params, 'thermostatMaxScale', 10.0);
minKBT = get_param(params, 'thermostatMinKBT', 1e-14);
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "targetkbt"
            targetKBT = val;
        case "strength"
            strength = val;
        case "minparticles"
            minParticles = val;
        case "maxscale"
            maxScale = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

if targetKBT <= 0 || strength <= 0
    vOut = v;
    info = empty_info();
    return;
end

ids = local_cell_ids(x, params, periodicX, periodicY);
Nc = params.Nx * params.Ny;
vOut = v;
scales = nan(Nc, 1);
kBTBefore = nan(Nc, 1);
kBTAfter = nan(Nc, 1);
nThermo = 0;

for c = 1:Nc
    pids = find(ids == c);
    if numel(pids) < minParticles
        continue;
    end
    u = mean(vOut(pids, :), 1);
    rel = vOut(pids, :) - u;
    current = 0.5 * mean(sum(rel.^2, 2));
    if ~isfinite(current) || current < minKBT
        continue;
    end
    targetMixed = (1 - strength) * current + strength * targetKBT;
    scale = sqrt(max(targetMixed, minKBT) / current);
    scale = min(max(scale, 1/maxScale), maxScale);
    vOut(pids, :) = u + scale * rel;
    scales(c) = scale;
    kBTBefore(c) = current;
    relAfter = vOut(pids, :) - mean(vOut(pids, :), 1);
    kBTAfter(c) = 0.5 * mean(sum(relAfter.^2, 2));
    nThermo = nThermo + 1;
end

info = struct();
info.enabled = true;
info.targetKBT = targetKBT;
info.strength = strength;
info.minParticlesPerCell = minParticles;
info.maxScale = maxScale;
info.nThermostattedCells = nThermo;
info.meanScale = mean(scales, 'omitnan');
info.minScale = min(scales, [], 'omitnan');
info.maxScaleApplied = max(scales, [], 'omitnan');
info.meanKBTBefore = mean(kBTBefore, 'omitnan');
info.meanKBTAfter = mean(kBTAfter, 'omitnan');
info.rmsVelocityChange = sqrt(mean(sum((vOut - v).^2, 2)));
end

function info = empty_info()
info = struct('enabled', false, 'targetKBT', NaN, 'strength', NaN, ...
    'minParticlesPerCell', NaN, 'maxScale', NaN, 'nThermostattedCells', 0, ...
    'meanScale', NaN, 'minScale', NaN, 'maxScaleApplied', NaN, ...
    'meanKBTBefore', NaN, 'meanKBTAfter', NaN, 'rmsVelocityChange', 0.0);
end

function ids = local_cell_ids(x, params, periodicX, periodicY)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;

xp = x(:,1);
yp = x(:,2);
if periodicX
    xp = mod(xp, Lx);
else
    xp = min(max(xp, 0), Lx - eps(Lx));
end
if periodicY
    yp = mod(yp, Ly);
else
    yp = min(max(yp, 0), Ly - eps(Ly));
end
ix = floor(xp / dx) + 1;
iy = floor(yp / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
ids = iy + Ny * (ix - 1);
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
