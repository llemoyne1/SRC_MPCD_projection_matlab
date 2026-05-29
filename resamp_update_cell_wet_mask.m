function [stateOut, info] = resamp_update_cell_wet_mask(state, params, varargin)
%RESAMP_UPDATE_CELL_WET_MASK Evolve or set state.cellWetMask.
%
%   [stateOut, info] = resamp_update_cell_wet_mask(state, params, ...)
%
% The default mode preserves an existing state.cellWetMask, or initializes an
% all-wet mask when none exists.  Threshold modes are intentionally simple and
% are meant as the first movable-fluid-domain layer for later injection,
% free-surface or mobile-solid cases.
%
% Modes:
%   'none'/'fixed'          keep previous mask, or all wet if absent;
%   'all'                   set all cells wet;
%   'manual'                use supplied 'cellWetMask';
%   'mass_threshold'        wet if M >= wetMassOnThreshold;
%   'mass_hysteresis'       dry->wet if M>=on, wet->dry if M<=off;
%   'particle_threshold'    wet if N >= wetParticleOnThreshold;
%   'particle_hysteresis'   same idea using population N.

if nargin < 2
    error('Usage: [stateOut, info] = resamp_update_cell_wet_mask(state, params, ...)');
end
Nx = params.Nx;
Ny = params.Ny;
mode = 'none';
cellWetMask = [];
wetMassOnThreshold = get_param(params, 'resampWetMassOnThreshold', 0.5 * get_param(params, 'resampTargetCellMass', get_param(params,'gamma',20)));
wetMassOffThreshold = get_param(params, 'resampWetMassOffThreshold', 0.1 * get_param(params, 'resampTargetCellMass', get_param(params,'gamma',20)));
wetParticleOnThreshold = get_param(params, 'resampWetParticleOnThreshold', max(1, ceil(0.5 * get_param(params,'gamma',20))));
wetParticleOffThreshold = get_param(params, 'resampWetParticleOffThreshold', 0);
periodicX = true;
periodicY = true;
activeMask = [];

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'mode','wetmaskupdatemode'}
            mode = lower(char(string(val)));
        case {'cellwetmask','wetmask'}
            cellWetMask = logical(val);
        case {'wetmassonthreshold','massonthreshold'}
            wetMassOnThreshold = val;
        case {'wetmassoffthreshold','massoffthreshold'}
            wetMassOffThreshold = val;
        case {'wetparticleonthreshold','particleonthreshold'}
            wetParticleOnThreshold = val;
        case {'wetparticleoffthreshold','particleoffthreshold'}
            wetParticleOffThreshold = val;
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case {'activemask','particleactivemask'}
            activeMask = logical(val(:));
        otherwise
            error('Unknown option: %s', string(key));
    end
end

stateOut = state;
[prevMask, prevInfo] = resamp_cell_wet_mask(state, params, 'mode', 'auto');
newMask = prevMask;
G = [];

switch mode
    case {'none','fixed','keep','preserve'}
        newMask = prevMask;
    case {'all','all_wet'}
        newMask = true(Nx, Ny);
    case {'dry','none_wet','all_dry'}
        newMask = false(Nx, Ny);
    case {'manual','explicit'}
        if isempty(cellWetMask)
            error('manual wet mask update requires cellWetMask.');
        end
        if ~isequal(size(cellWetMask), [Nx, Ny])
            error('cellWetMask must have size Nx-by-Ny.');
        end
        newMask = logical(cellWetMask);
    case {'mass_threshold','mass'}
        G = deposit(state, params, periodicX, periodicY, activeMask);
        newMask = G.M >= wetMassOnThreshold;
    case {'mass_hysteresis','mass_hyst'}
        G = deposit(state, params, periodicX, periodicY, activeMask);
        newMask = prevMask;
        newMask(~prevMask & G.M >= wetMassOnThreshold) = true;
        newMask(prevMask & G.M <= wetMassOffThreshold) = false;
    case {'particle_threshold','population_threshold','n_threshold'}
        G = deposit(state, params, periodicX, periodicY, activeMask);
        newMask = G.N >= wetParticleOnThreshold;
    case {'particle_hysteresis','population_hysteresis','n_hysteresis'}
        G = deposit(state, params, periodicX, periodicY, activeMask);
        newMask = prevMask;
        newMask(~prevMask & G.N >= wetParticleOnThreshold) = true;
        newMask(prevMask & G.N <= wetParticleOffThreshold) = false;
    otherwise
        error('Unknown wetMaskUpdateMode: %s', mode);
end

stateOut.cellWetMask = logical(newMask);
info = struct();
info.kind = 'update_cell_wet_mask';
info.mode = mode;
info.nWetCellsBefore = prevInfo.nWetCells;
info.nWetCellsAfter = nnz(newMask);
info.nDryCellsAfter = numel(newMask) - nnz(newMask);
info.wetFractionAfter = nnz(newMask) / max(numel(newMask), 1);
info.nCellsBecameWet = nnz(~prevMask & newMask);
info.nCellsBecameDry = nnz(prevMask & ~newMask);
info.wetMassOnThreshold = wetMassOnThreshold;
info.wetMassOffThreshold = wetMassOffThreshold;
info.wetParticleOnThreshold = wetParticleOnThreshold;
info.wetParticleOffThreshold = wetParticleOffThreshold;
if ~isempty(G) && any(newMask(:))
    info.NMinWetAfter = min(double(G.N(newMask)), [], 'omitnan');
    info.NMaxWetAfter = max(double(G.N(newMask)), [], 'omitnan');
    info.MMinWetAfter = min(double(G.M(newMask)), [], 'omitnan');
    info.MMaxWetAfter = max(double(G.M(newMask)), [], 'omitnan');
else
    info.NMinWetAfter = NaN;
    info.NMaxWetAfter = NaN;
    info.MMinWetAfter = NaN;
    info.MMaxWetAfter = NaN;
end
end

function G = deposit(state, params, periodicX, periodicY, activeMask)
if isempty(activeMask)
    activeMask = resamp_active_mask(state);
end
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'activeMask', activeMask);
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
