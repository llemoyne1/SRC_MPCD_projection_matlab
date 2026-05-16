function [vOut, info] = projection_apply_global_momentum_correction(vBefore, dvRaw, params, varargin)
%PROJECTION_APPLY_GLOBAL_MOMENTUM_CORRECTION Remove global momentum drift from a numerical kick.
%
%   [vOut, info] = projection_apply_global_momentum_correction(vBefore, dvRaw, params)
%
% The correction is exact for unit-mass particles:
%       dvCorrected = dvRaw - mean(dvRaw(active,:),1)
% so that sum_i dvCorrected_i = 0 over active particles.  It is intended
% for numerical projection stages only (Q6, Q9, cleanup), not for body
% forces, walls, or imposed mean-flow control.

if nargin < 3 || isempty(params)
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
Np = size(vBefore, 1);
if isempty(activeMask)
    activeMask = true(Np, 1);
end
if numel(activeMask) ~= Np
    error('activeMask must have Np entries.');
end

info = empty_info(stage, params, Np, nnz(activeMask));
if Np == 0 || nnz(activeMask) == 0
    vOut = vBefore + dvRaw;
    return;
end

enabled = logical(get_param(params, 'projectionMomentumCorrectionEnable', true));
mode = char(string(get_param(params, 'projectionMomentumCorrectionMode', 'particle_global_exact')));
info.enabled = enabled;
info.mode = mode;

rawDeltaP = sum(dvRaw(activeMask,:), 1);
rawMean = rawDeltaP ./ nnz(activeMask);

if enabled && strcmpi(mode, 'particle_global_exact')
    dvCorr = dvRaw;
    dvCorr(activeMask,:) = dvCorr(activeMask,:) - rawMean;
    correctionVector = -rawMean;
else
    dvCorr = dvRaw;
    correctionVector = [0 0];
end

vOut = vBefore + dvCorr;
residualDeltaP = sum((vOut(activeMask,:) - vBefore(activeMask,:)), 1);
correctionDeltaP = sum((dvCorr(activeMask,:) - dvRaw(activeMask,:)), 1);

info.applied = enabled && strcmpi(mode, 'particle_global_exact');
info.rawDeltaP = rawDeltaP;
info.rawDeltaPNorm = norm(rawDeltaP);
info.rawDeltaMeanV = rawMean;
info.rawDeltaMeanVNorm = norm(rawMean);
info.correctionVector = correctionVector;
info.correctionVectorNorm = norm(correctionVector);
info.correctionDeltaP = correctionDeltaP;
info.correctionDeltaPNorm = norm(correctionDeltaP);
info.residualDeltaP = residualDeltaP;
info.residualDeltaPNorm = norm(residualDeltaP);
info.residualMeanV = residualDeltaP ./ nnz(activeMask);
info.residualMeanVNorm = norm(info.residualMeanV);
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
info.rawDeltaP = [0 0];
info.rawDeltaPNorm = 0;
info.rawDeltaMeanV = [0 0];
info.rawDeltaMeanVNorm = 0;
info.correctionVector = [0 0];
info.correctionVectorNorm = 0;
info.correctionDeltaP = [0 0];
info.correctionDeltaPNorm = 0;
info.residualDeltaP = [0 0];
info.residualDeltaPNorm = 0;
info.residualMeanV = [0 0];
info.residualMeanVNorm = 0;
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
