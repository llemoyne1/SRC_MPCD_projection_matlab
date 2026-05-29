function [stateOut, diag] = resamp_extract_overpopulated_particles(state, params, varargin)
%RESAMP_EXTRACT_OVERPOPULATED_PARTICLES Recycle particles from overpopulated cells.
%
%   [stateOut, diag] = resamp_extract_overpopulated_particles(state, params)
%
% The routine deactivates storage slots in cells with N > NMax, bringing
% those cells back toward NTarget.  It does not attempt to preserve local mass
% or momentum by itself.  The intended use is:
%
%   1. save the pre-extraction hydrodynamic field U_target;
%   2. extract overpopulated particles;
%   3. insert support particles in underpopulated cells;
%   4. run resamp_local_mass_moment_remap(...,'targetVelocityMode','grid',...)
%
% so that exact cell mass and target cell momentum are restored after the
% population edit.
%
% Supported selection modes:
%   closest_to_cell_mean  remove particles with velocity closest to cell mean
%   random                remove random particles from the cell
%   lightest              remove lightest particles first
%   heaviest              remove heaviest particles first

validate_state(state);
Nx = params.Nx;
Ny = params.Ny;
Nc = Nx * Ny;

gamma = get_param(params, 'gamma', 20);
NTarget = get_param(params, 'resampNTarget', gamma);
NMin = get_param(params, 'resampNMin', ceil(0.5 * gamma));
NMax = get_param(params, 'resampNMax', ceil(1.5 * gamma));
selectionMode = 'closest_to_cell_mean';
periodicX = true;
periodicY = true;
computeDiagnostics = true;
minMass = eps;

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case 'ntarget'
            NTarget = val;
        case 'nmin'
            NMin = val;
        case 'nmax'
            NMax = val;
        case {'selectionmode','extractselectionmode'}
            selectionMode = lower(char(string(val)));
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case {'minmass','minm'}
            minMass = val;
        case 'computediagnostics'
            computeDiagnostics = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

if NTarget < 0 || NMax < 0
    error('NTarget and NMax must be non-negative.');
end
if NTarget > NMax
    error('NTarget must be <= NMax for extraction.');
end

stateOut = state;
if ~isfield(stateOut, 'active') || isempty(stateOut.active)
    error('state must be converted to a particle pool with resamp_enable_particle_pool before extraction.');
end
activeBefore = resamp_active_mask(stateOut);
Gbefore = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', minMass, 'activeMask', activeBefore);
cellId = Gbefore.cellId;
Nvec = reshape(Gbefore.N.', [Nc, 1]);
UxVec = reshape(Gbefore.Ux.', [Nc, 1]);
UyVec = reshape(Gbefore.Uy.', [Nc, 1]);

overCells = find(Nvec > NMax);
poorCells = find(Nvec < NMin);
extractedPerCell = zeros(Nc, 1);
extractedIds = [];

for jj = 1:numel(overCells)
    c = overCells(jj);
    nExtract = max(0, round(Nvec(c) - NTarget));
    if nExtract <= 0
        continue;
    end
    ids = find(activeBefore & cellId == c);
    if isempty(ids)
        continue;
    end
    nExtract = min(nExtract, max(0, numel(ids) - max(0, round(NTarget))));
    if nExtract <= 0
        continue;
    end
    chosen = choose_particles(ids, stateOut, UxVec(c), UyVec(c), nExtract, selectionMode);
    if isempty(chosen)
        continue;
    end
    stateOut.active(chosen) = false;
    stateOut.m(chosen) = 0;
    stateOut.v(chosen,:) = 0;
    stateOut.x(chosen,:) = 0;
    extractedPerCell(c) = numel(chosen);
    extractedIds = [extractedIds; chosen(:)]; %#ok<AGROW>
end

activeAfter = resamp_active_mask(stateOut);
stateOut.Nactive = nnz(activeAfter);
stateOut.Ncapacity = size(stateOut.x, 1);

Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', minMass, 'activeMask', activeAfter);
Nafter = reshape(Gafter.N.', [Nc, 1]);

diag = struct();
diag.kind = 'extract_overpopulated_particles';
diag.NTarget = NTarget;
diag.NMin = NMin;
diag.NMax = NMax;
diag.selectionMode = selectionMode;
diag.nCells = Nc;
diag.nOverCellsBefore = numel(overCells);
diag.nPoorCellsBefore = numel(poorCells);
diag.nEmptyCellsBefore = nnz(Nvec == 0);
diag.nExtractedParticles = sum(extractedPerCell);
diag.nCellsExtracted = nnz(extractedPerCell > 0);
diag.NactiveBefore = nnz(activeBefore);
diag.NactiveAfter = nnz(activeAfter);
diag.Ncapacity = size(stateOut.x, 1);
diag.NfreeBefore = nnz(~activeBefore);
diag.NfreeAfter = nnz(~activeAfter);
diag.nPoorCellsAfter = nnz(Nafter < NMin);
diag.nEmptyCellsAfter = nnz(Nafter == 0);
diag.nOverCellsAfter = nnz(Nafter > NMax);
diag.NMinAfter = min(Nafter);
diag.NMaxAfter = max(Nafter);
diag.NStdAfter = std(double(Nafter), 0, 'omitnan');
diag.extractedPerCellGrid = [];
diag.extractedIds = [];
diag.Gbefore = [];
diag.Gafter = [];
if computeDiagnostics
    diag.extractedPerCellGrid = reshape(extractedPerCell, [Ny, Nx]).';
    diag.extractedIds = extractedIds;
    diag.Gbefore = Gbefore;
    diag.Gafter = Gafter;
end
end

function idsChosen = choose_particles(ids, state, ux, uy, nExtract, selectionMode)
ids = ids(:);
if nExtract <= 0 || isempty(ids)
    idsChosen = [];
    return;
end
nExtract = min(nExtract, numel(ids));
switch selectionMode
    case {'closest_to_cell_mean','closest','closest_to_mean'}
        dvx = state.v(ids,1) - ux;
        dvy = state.v(ids,2) - uy;
        score = dvx.^2 + dvy.^2;
        [~, order] = sort(score, 'ascend');
        idsChosen = ids(order(1:nExtract));
    case 'random'
        order = randperm(numel(ids), nExtract);
        idsChosen = ids(order);
    case 'lightest'
        [~, order] = sort(state.m(ids), 'ascend');
        idsChosen = ids(order(1:nExtract));
    case 'heaviest'
        [~, order] = sort(state.m(ids), 'descend');
        idsChosen = ids(order(1:nExtract));
    otherwise
        error('Unknown extractSelectionMode: %s', selectionMode);
end
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must contain x, v and m.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
if numel(state.m) ~= size(state.x,1)
    error('state.m must have one entry per particle slot.');
end
if any(~isfinite(state.m(:))) || any(state.m(:) < 0)
    error('state.m must contain finite non-negative masses.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
