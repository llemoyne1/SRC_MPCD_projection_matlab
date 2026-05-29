function [wetMask, info] = resamp_cell_wet_mask(state, params, varargin)
%RESAMP_CELL_WET_MASK Return the current wet/fluid-active cell mask.
%
%   wetMask = resamp_cell_wet_mask(state, params)
%
% The weighted-resampling formulation must not impose full cell mass on
% cells that are not part of the current fluid domain.  This helper provides
% one common wet/non-wet tag for remap, insertion, extraction and diagnostics.
%
% Resolution order for the default mode is:
%   1. explicit 'cellWetMask'/'wetMask' option;
%   2. state.cellWetMask, if present;
%   3. params.cellWetMask, if present;
%   4. all cells wet.
%
% Optional modes:
%   'all'                  every cell is wet;
%   'state'                require state.cellWetMask;
%   'params'               require params.cellWetMask;
%   'mass_threshold'       wet if deposited mass >= wetMassThreshold;
%   'particle_threshold'   wet if deposited population >= wetParticleThreshold.
%
% The mask has size Nx-by-Ny and is logical.  It can evolve between steps by
% updating state.cellWetMask, or by using resamp_update_cell_wet_mask.

if nargin < 2
    error('Usage: wetMask = resamp_cell_wet_mask(state, params, ...)');
end
if ~isfield(params, 'Nx') || ~isfield(params, 'Ny')
    error('params must contain Nx and Ny.');
end
Nx = params.Nx;
Ny = params.Ny;
mode = 'auto';
cellWetMask = [];
wetMassThreshold = get_param(params, 'resampWetMassThreshold', eps);
wetParticleThreshold = get_param(params, 'resampWetParticleThreshold', 1);
periodicX = true;
periodicY = true;
activeMask = [];
G = [];

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'cellwetmask','wetmask','fluidmask','activecellmask'}
            % Empty masks are common when callers simply forward an optional
            % argument.  They must not force explicit/manual mode; otherwise
            % all-wet default initialization fails before state.cellWetMask
            % exists.  A non-empty mask remains an explicit user request.
            if isempty(val)
                cellWetMask = [];
            else
                cellWetMask = logical(val);
                mode = 'explicit';
            end
        case {'mode','wetmaskmode'}
            mode = lower(char(string(val)));
        case {'wetmassthreshold','massthreshold'}
            wetMassThreshold = val;
        case {'wetparticlethreshold','particlethreshold'}
            wetParticleThreshold = val;
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case {'activemask','particleactivemask'}
            activeMask = logical(val(:));
        case {'deposit','grid','g'}
            G = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

switch mode
    case {'explicit','manual'}
        wetMask = validate_mask(cellWetMask, Nx, Ny, 'cellWetMask');
    case {'auto','state_or_params_or_all','default'}
        if ~isempty(cellWetMask)
            wetMask = validate_mask(cellWetMask, Nx, Ny, 'cellWetMask');
        elseif isstruct(state) && isfield(state, 'cellWetMask') && ~isempty(state.cellWetMask)
            wetMask = validate_mask(state.cellWetMask, Nx, Ny, 'state.cellWetMask');
        elseif isfield(params, 'cellWetMask') && ~isempty(params.cellWetMask)
            wetMask = validate_mask(params.cellWetMask, Nx, Ny, 'params.cellWetMask');
        else
            wetMask = true(Nx, Ny);
        end
    case {'all','full','all_wet'}
        wetMask = true(Nx, Ny);
    case {'none','dry','all_dry'}
        wetMask = false(Nx, Ny);
    case 'state'
        if ~isstruct(state) || ~isfield(state, 'cellWetMask') || isempty(state.cellWetMask)
            error('mode=''state'' requires state.cellWetMask.');
        end
        wetMask = validate_mask(state.cellWetMask, Nx, Ny, 'state.cellWetMask');
    case 'params'
        if ~isfield(params, 'cellWetMask') || isempty(params.cellWetMask)
            error('mode=''params'' requires params.cellWetMask.');
        end
        wetMask = validate_mask(params.cellWetMask, Nx, Ny, 'params.cellWetMask');
    case {'mass_threshold','mass','by_mass'}
        G = ensure_deposit(G, state, params, periodicX, periodicY, activeMask);
        wetMask = G.M >= wetMassThreshold;
    case {'particle_threshold','population_threshold','n_threshold','by_population'}
        G = ensure_deposit(G, state, params, periodicX, periodicY, activeMask);
        wetMask = G.N >= wetParticleThreshold;
    otherwise
        error('Unknown wetMaskMode: %s', mode);
end
wetMask = logical(wetMask);

info = struct();
info.kind = 'cell_wet_mask';
info.mode = mode;
info.nWetCells = nnz(wetMask);
info.nDryCells = numel(wetMask) - nnz(wetMask);
info.wetFraction = nnz(wetMask) / max(numel(wetMask), 1);
info.wetMassThreshold = wetMassThreshold;
info.wetParticleThreshold = wetParticleThreshold;
end

function wetMask = validate_mask(wetMask, Nx, Ny, label)
if isempty(wetMask)
    error('%s is empty.', label);
end
if ~isequal(size(wetMask), [Nx, Ny])
    error('%s must have size Nx-by-Ny = [%d %d].', label, Nx, Ny);
end
wetMask = logical(wetMask);
end

function G = ensure_deposit(G, state, params, periodicX, periodicY, activeMask)
if isstruct(G) && isfield(G, 'N') && isfield(G, 'M')
    return;
end
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
