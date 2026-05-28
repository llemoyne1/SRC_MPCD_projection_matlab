function [vOut, info] = resamp_apply_cell_thermostat_weighted(x, v, m, params, varargin)
%RESAMP_APPLY_CELL_THERMOSTAT_WEIGHTED Momentum-preserving weighted thermostat.
%
% Per cell, rescales fluctuations around the weighted mean velocity:
%      v_i <- U_cell + s_cell (v_i - U_cell)
% with U_cell = sum(m_i v_i)/sum(m_i).  This preserves weighted cell
% momentum exactly up to roundoff.

if nargin < 4
    error('Usage: [vOut, info] = resamp_apply_cell_thermostat_weighted(x, v, m, params, ...)');
end
periodicX = true;
periodicY = true;
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

m = m(:);
if targetKBT <= 0 || strength <= 0 || isempty(v)
    vOut = v;
    info = empty_info();
    return;
end
if numel(m) ~= size(v,1)
    error('m must have one entry per particle.');
end

ids = resamp_cell_ids_periodic(x, params, 'periodicX', periodicX, 'periodicY', periodicY);
Nc = params.Nx * params.Ny;
Np = size(v, 1);

nCell = accumarray(ids, 1, [Nc 1], @sum, 0);
M = accumarray(ids, m, [Nc 1], @sum, 0);
Px = accumarray(ids, m .* v(:,1), [Nc 1], @sum, 0);
Py = accumarray(ids, m .* v(:,2), [Nc 1], @sum, 0);

Ux = zeros(Nc,1);
Uy = zeros(Nc,1);
populated = M > eps;
Ux(populated) = Px(populated) ./ M(populated);
Uy(populated) = Py(populated) ./ M(populated);

UxP = Ux(ids);
UyP = Uy(ids);
relx = v(:,1) - UxP;
rely = v(:,2) - UyP;
rel2 = relx.^2 + rely.^2;

sumMRel2 = accumarray(ids, m .* rel2, [Nc 1], @sum, 0);
kBTBefore = nan(Nc,1);
kBTBefore(populated) = 0.5 * sumMRel2(populated) ./ M(populated);

valid = nCell >= minParticles & populated & isfinite(kBTBefore) & kBTBefore > minKBT;
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

function v = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end
