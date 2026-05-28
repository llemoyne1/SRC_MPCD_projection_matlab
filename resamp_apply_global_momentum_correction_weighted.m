function [vOut, info] = resamp_apply_global_momentum_correction_weighted(vBefore, dvRaw, m, params, varargin)
%RESAMP_APPLY_GLOBAL_MOMENTUM_CORRECTION_WEIGHTED Remove weighted global momentum drift.
%
% This is the mass-weighted analogue of projection_apply_global_momentum_correction.
% It subtracts a uniform velocity offset from active particles so that
%       sum_i m_i dv_i = 0
% over the active set.

if nargin < 4 || isempty(params)
    params = struct();
end
stage = 'projection';
activeMask = [];
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "stage"
            stage = char(string(val));
        case {"activemask", "mask"}
            activeMask = logical(val(:));
        otherwise
            error('Unknown option: %s', string(key));
    end
end

if ~isequal(size(vBefore), size(dvRaw)) || size(vBefore,2) ~= 2
    error('vBefore and dvRaw must be matching Np-by-2 arrays.');
end
m = m(:);
Np = size(vBefore, 1);
if numel(m) ~= Np
    error('m must have one entry per particle.');
end
if isempty(activeMask)
    activeMask = true(Np, 1);
end
if numel(activeMask) ~= Np
    error('activeMask must have Np entries.');
end
activeMask = activeMask & isfinite(m) & m > 0;

info = empty_info(stage, params, Np, nnz(activeMask));
if Np == 0 || nnz(activeMask) == 0
    vOut = vBefore + dvRaw;
    return;
end

enabled = logical(get_param(params, 'projectionMomentumCorrectionEnable', true));
mode = char(string(get_param(params, 'projectionMomentumCorrectionMode', 'particle_global_exact')));
info.enabled = enabled;
info.mode = mode;

ma = m(activeMask);
Mactive = sum(ma, 'omitnan');
rawDeltaP = [sum(ma .* dvRaw(activeMask,1), 'omitnan'), ...
             sum(ma .* dvRaw(activeMask,2), 'omitnan')];
rawMeanDV = rawDeltaP ./ max(Mactive, eps);

if enabled && strcmpi(mode, 'particle_global_exact')
    dvCorr = dvRaw;
    dvCorr(activeMask,:) = dvCorr(activeMask,:) - rawMeanDV;
    correctionVector = -rawMeanDV;
else
    dvCorr = dvRaw;
    correctionVector = [0 0];
end

vOut = vBefore + dvCorr;
residualDeltaP = [sum(ma .* (vOut(activeMask,1) - vBefore(activeMask,1)), 'omitnan'), ...
                  sum(ma .* (vOut(activeMask,2) - vBefore(activeMask,2)), 'omitnan')];
correctionDeltaP = [sum(ma .* (dvCorr(activeMask,1) - dvRaw(activeMask,1)), 'omitnan'), ...
                    sum(ma .* (dvCorr(activeMask,2) - dvRaw(activeMask,2)), 'omitnan')];

info.applied = enabled && strcmpi(mode, 'particle_global_exact');
info.totalActiveMass = Mactive;
info.rawDeltaP = rawDeltaP;
info.rawDeltaPNorm = norm(rawDeltaP);
info.rawMeanDV = rawMeanDV;
info.rawMeanDVNorm = norm(rawMeanDV);
info.correctionVector = correctionVector;
info.correctionVectorNorm = norm(correctionVector);
info.correctionDeltaP = correctionDeltaP;
info.correctionDeltaPNorm = norm(correctionDeltaP);
info.residualDeltaP = residualDeltaP;
info.residualDeltaPNorm = norm(residualDeltaP);
info.residualMeanDV = residualDeltaP ./ max(Mactive, eps);
info.residualMeanDVNorm = norm(info.residualMeanDV);
info.dvRawRms = sqrt(mean(sum(dvRaw(activeMask,:).^2, 2), 'omitnan'));
info.dvCorrectedRms = sqrt(mean(sum(dvCorr(activeMask,:).^2, 2), 'omitnan'));
end

function info = empty_info(stage, params, Np, Na)
info = struct();
info.stage = stage;
info.enabled = logical(get_param(params, 'projectionMomentumCorrectionEnable', true));
info.applied = false;
info.mode = char(string(get_param(params, 'projectionMomentumCorrectionMode', 'particle_global_exact')));
info.nParticles = Np;
info.nActiveParticles = Na;
info.totalActiveMass = 0;
info.rawDeltaP = [0 0];
info.rawDeltaPNorm = 0;
info.rawMeanDV = [0 0];
info.rawMeanDVNorm = 0;
info.correctionVector = [0 0];
info.correctionVectorNorm = 0;
info.correctionDeltaP = [0 0];
info.correctionDeltaPNorm = 0;
info.residualDeltaP = [0 0];
info.residualDeltaPNorm = 0;
info.residualMeanDV = [0 0];
info.residualMeanDVNorm = 0;
info.dvRawRms = 0;
info.dvCorrectedRms = 0;
end

function v = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end
