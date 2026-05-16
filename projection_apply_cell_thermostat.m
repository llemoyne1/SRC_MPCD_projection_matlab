function [vOut, info] = projection_apply_cell_thermostat(x, v, params, varargin)
%PROJECTION_APPLY_CELL_THERMOSTAT Momentum-preserving local thermostat.
%
%   [vOut, info] = projection_apply_cell_thermostat(x, v, params, ...)
%
% For each populated grid cell, rescales particle velocity fluctuations
% around the cell mean velocity:
%      v_i <- U_cell + s_cell * (v_i - U_cell)
%
% This preserves the cell momentum exactly up to roundoff.  This version is
% vectorized using accumarray and is intended to replace the legacy loop
% implementation that performed a costly find(ids==c) for every cell.
%
% Supported options:
%   'periodicX'     true/false
%   'periodicY'     true/false
%   'targetKBT'     scalar
%   'strength'      scalar in [0,1] usually
%   'minParticles'  minimum population per cell
%   'maxScale'      bound on thermostat rescale factor
%
% Notes:
% - The operation is local to each cell.
% - The hydrodynamic cell velocity is preserved.
% - The global momentum is therefore also preserved, up to roundoff.

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

if targetKBT <= 0 || strength <= 0 || isempty(v)
    vOut = v;
    info = empty_info();
    return;
end

ids = local_cell_ids(x, params, periodicX, periodicY);
Nc = params.Nx * params.Ny;
Np = size(v, 1);

% Population and cell means. accumarray is much faster than find(ids==c)
% inside a cell loop, especially for moderately large particle counts.
nCell = accumarray(ids, 1, [Nc 1], @sum, 0);
sumVx = accumarray(ids, v(:,1), [Nc 1], @sum, 0);
sumVy = accumarray(ids, v(:,2), [Nc 1], @sum, 0);

Ux = zeros(Nc,1);
Uy = zeros(Nc,1);
populated = nCell > 0;
Ux(populated) = sumVx(populated) ./ nCell(populated);
Uy(populated) = sumVy(populated) ./ nCell(populated);

UxP = Ux(ids);
UyP = Uy(ids);
relx = v(:,1) - UxP;
rely = v(:,2) - UyP;
rel2 = relx.^2 + rely.^2;

sumRel2 = accumarray(ids, rel2, [Nc 1], @sum, 0);
kBTBefore = nan(Nc,1);
kBTBefore(populated) = 0.5 * sumRel2(populated) ./ nCell(populated);

valid = nCell >= minParticles & isfinite(kBTBefore) & kBTBefore > minKBT;

scales = nan(Nc,1);
scaleApply = ones(Nc,1);
if any(valid)
    targetMixed = (1 - strength) .* kBTBefore(valid) + strength .* targetKBT;
    s = sqrt(max(targetMixed, minKBT) ./ kBTBefore(valid));
    s = min(max(s, 1/maxScale), maxScale);
    scales(valid) = s;
    scaleApply(valid) = s;
end

sP = scaleApply(ids);
validP = valid(ids);
vOut = v;
if any(validP)
    vOut(validP,1) = UxP(validP) + sP(validP) .* relx(validP);
    vOut(validP,2) = UyP(validP) + sP(validP) .* rely(validP);
end

% Since the cell mean is preserved, kBT_after = scale^2 kBT_before for
% thermostatted cells.  Avoid a second deposition pass.
kBTAfter = nan(Nc,1);
kBTAfter(valid) = (scales(valid).^2) .* kBTBefore(valid);

info = struct();
info.enabled = true;
info.targetKBT = targetKBT;
info.strength = strength;
info.minParticlesPerCell = minParticles;
info.maxScale = maxScale;
info.nThermostattedCells = nnz(valid);
info.meanScale = mean(scales, 'omitnan');
info.minScale = min(scales, [], 'omitnan');
info.maxScaleApplied = max(scales, [], 'omitnan');
info.meanKBTBefore = mean(kBTBefore(valid), 'omitnan');
info.meanKBTAfter = mean(kBTAfter(valid), 'omitnan');
if Np > 0
    info.rmsVelocityChange = sqrt(mean(sum((vOut - v).^2, 2), 'omitnan'));
else
    info.rmsVelocityChange = 0.0;
end
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
