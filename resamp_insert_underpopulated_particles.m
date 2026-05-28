function [stateOut, diag] = resamp_insert_underpopulated_particles(state, params, varargin)
%RESAMP_INSERT_UNDERPOPULATED_PARTICLES Fill poor cells from a preallocated pool.
%
% This routine does not remove/fuse overpopulated particles.  It only
% activates free slots when a cell has N < NMin and fills it up to NTarget.
% New particles receive nominal mass; the existing local remap stage is then
% expected to impose exact target cell mass and preserve target cell velocity.
%
% Velocity memory:
%   state.uMemUx/uMemUy/uMemValid store the last reliable cell velocity.
%   Reliable cells are those with N >= memoryMinParticles before insertion.
%   Empty cells are initialized from memory when available, otherwise zero.

validate_state(state);
Nx = params.Nx;
Ny = params.Ny;
Nc = Nx * Ny;
Lx = params.Lx;
Ly = params.Ly;
dx = Lx / Nx;
dy = Ly / Ny;

gamma = get_param(params, 'gamma', 20);
m0 = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', 1.0));
kBT = get_param(params, 'kBT', 0.0);
NTarget = get_param(params, 'resampNTarget', gamma);
NMin = get_param(params, 'resampNMin', ceil(0.5 * gamma));
NMax = get_param(params, 'resampNMax', ceil(1.5 * gamma));
memoryMinParticles = get_param(params, 'resampMemoryMinParticles', NMin);
insertVelocityMode = 'current_or_memory_pairwise';
periodicX = true;
periodicY = true;
computeDiagnostics = true;
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
        case 'memoryminparticles'
            memoryMinParticles = val;
        case {'particlemass','m0'}
            m0 = val;
        case 'kbt'
            kBT = val;
        case 'insertvelocitymode'
            insertVelocityMode = lower(char(string(val)));
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case 'computediagnostics'
            computeDiagnostics = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

stateOut = ensure_memory_fields(state, Nx, Ny);
activeBefore = resamp_active_mask(stateOut);
if ~isfield(stateOut, 'active') || isempty(stateOut.active)
    error('state must be converted to a particle pool with resamp_enable_particle_pool before insertion.');
end

Gbefore = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', eps, 'activeMask', activeBefore);
Nvec = reshape(Gbefore.N.', [Nc, 1]);
UxVec = reshape(Gbefore.Ux.', [Nc, 1]);
UyVec = reshape(Gbefore.Uy.', [Nc, 1]);
validVec = reshape(Gbefore.valid.', [Nc, 1]);

% Update memory from reliable cells before inserting new support.
reliable = Nvec >= memoryMinParticles & validVec;
if any(reliable)
    uMemUxVec = reshape(stateOut.uMemUx.', [Nc, 1]);
    uMemUyVec = reshape(stateOut.uMemUy.', [Nc, 1]);
    uMemValidVec = reshape(stateOut.uMemValid.', [Nc, 1]);
    uMemUxVec(reliable) = UxVec(reliable);
    uMemUyVec(reliable) = UyVec(reliable);
    uMemValidVec(reliable) = true;
    stateOut.uMemUx = reshape(uMemUxVec, [Ny, Nx]).';
    stateOut.uMemUy = reshape(uMemUyVec, [Ny, Nx]).';
    stateOut.uMemValid = reshape(uMemValidVec, [Ny, Nx]).';
end

uMemUxVec = reshape(stateOut.uMemUx.', [Nc, 1]);
uMemUyVec = reshape(stateOut.uMemUy.', [Nc, 1]);
uMemValidVec = reshape(stateOut.uMemValid.', [Nc, 1]);

poorCells = find(Nvec < NMin);
overCells = find(Nvec > NMax);
freeSlots = find(~activeBefore);
freePtr = 1;
insertedPerCell = zeros(Nc, 1);
capacityHit = false;

for jj = 1:numel(poorCells)
    c = poorCells(jj);
    nAddWanted = max(0, round(NTarget - Nvec(c)));
    if nAddWanted <= 0
        continue;
    end
    nFree = numel(freeSlots) - freePtr + 1;
    if nFree <= 0
        capacityHit = true;
        break;
    end
    nAdd = min(nAddWanted, nFree);
    if nAdd < nAddWanted
        capacityHit = true;
    end
    ids = freeSlots(freePtr:(freePtr + nAdd - 1));
    freePtr = freePtr + nAdd;

    [ix, iy] = cell_index_from_id(c, Ny);
    stateOut.x(ids,1) = (ix - 1 + rand(nAdd, 1)) * dx;
    stateOut.x(ids,2) = (iy - 1 + rand(nAdd, 1)) * dy;

    switch insertVelocityMode
        case {'current_or_memory_pairwise','current_or_memory'}
            if Nvec(c) > 0 && validVec(c)
                U = [UxVec(c), UyVec(c)];
            elseif uMemValidVec(c)
                U = [uMemUxVec(c), uMemUyVec(c)];
            else
                U = [0, 0];
            end
        case {'memory_pairwise','memory'}
            if uMemValidVec(c)
                U = [uMemUxVec(c), uMemUyVec(c)];
            else
                U = [0, 0];
            end
        case {'zero','zero_pairwise'}
            U = [0, 0];
        otherwise
            error('Unknown insertVelocityMode: %s', insertVelocityMode);
    end

    fluct = zeros(nAdd, 2);
    if kBT > 0
        fluct = sqrt(kBT) * randn(nAdd, 2);
        % Enforce zero mean over inserted particles so that insertion does
        % not add a systematic momentum beyond the selected memory velocity.
        fluct(:,1) = fluct(:,1) - mean(fluct(:,1), 'omitnan');
        fluct(:,2) = fluct(:,2) - mean(fluct(:,2), 'omitnan');
    end
    stateOut.v(ids,1) = U(1) + fluct(:,1);
    stateOut.v(ids,2) = U(2) + fluct(:,2);
    stateOut.m(ids) = m0;
    stateOut.active(ids) = true;
    insertedPerCell(c) = nAdd;

    if capacityHit
        break;
    end
end

activeAfter = resamp_active_mask(stateOut);
stateOut.Nactive = nnz(activeAfter);
stateOut.Ncapacity = size(stateOut.x, 1);
Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', eps, 'activeMask', activeAfter);
Nafter = reshape(Gafter.N.', [Nc, 1]);

diag = struct();
diag.kind = 'insert_underpopulated_particles';
diag.NTarget = NTarget;
diag.NMin = NMin;
diag.NMax = NMax;
diag.memoryMinParticles = memoryMinParticles;
diag.insertVelocityMode = insertVelocityMode;
diag.nCells = Nc;
diag.nPoorCellsBefore = numel(poorCells);
diag.nOverCellsBefore = numel(overCells);
diag.nEmptyCellsBefore = nnz(Nvec == 0);
diag.nInsertedParticles = sum(insertedPerCell);
diag.nCellsInserted = nnz(insertedPerCell > 0);
diag.capacityHit = capacityHit;
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
diag.nMemoryReliableUpdated = nnz(reliable);
diag.nMemoryValid = nnz(stateOut.uMemValid);
diag.insertedPerCellGrid = [];
diag.Gbefore = [];
diag.Gafter = [];
if computeDiagnostics
    diag.insertedPerCellGrid = reshape(insertedPerCell, [Ny, Nx]).';
    diag.Gbefore = Gbefore;
    diag.Gafter = Gafter;
end
end

function state = ensure_memory_fields(state, Nx, Ny)
if ~isfield(state, 'uMemUx') || isempty(state.uMemUx) || ~isequal(size(state.uMemUx), [Nx, Ny])
    state.uMemUx = zeros(Nx, Ny);
end
if ~isfield(state, 'uMemUy') || isempty(state.uMemUy) || ~isequal(size(state.uMemUy), [Nx, Ny])
    state.uMemUy = zeros(Nx, Ny);
end
if ~isfield(state, 'uMemValid') || isempty(state.uMemValid) || ~isequal(size(state.uMemValid), [Nx, Ny])
    state.uMemValid = false(Nx, Ny);
end
end

function [ix, iy] = cell_index_from_id(c, Ny)
ix = floor((c - 1) / Ny) + 1;
iy = c - Ny * (ix - 1);
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
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
