function [stateOut, diag] = resamp_activate_particles_in_cell(state, params, cellIds, varargin)
%RESAMP_ACTIVATE_PARTICLES_IN_CELL Activate latent/free slots as fluid in cells.
%
%   Activates nPerCell particles in each requested cell.  Slots may be latent
%   (particleRole==2) or free pool (particleRole==0 / active=false).

if nargin < 3
    error('Usage: [stateOut, diag] = resamp_activate_particles_in_cell(state, params, cellIds, ...)');
end
cellIds = cellIds(:);
cellIds = cellIds(isfinite(cellIds) & cellIds >= 1 & cellIds <= params.Nx*params.Ny);
cellIds = round(cellIds(:));

nPerCell = get_param(params, 'resampNTarget', get_param(params, 'gamma', 20));
m0 = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', 1.0));
U = [0 0];
kBT = get_param(params, 'kBT', 0.0);
preferLatent = true;
periodicX = true;
periodicY = true;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'npercell','nparticles','ntarget'}
            nPerCell = val;
        case {'particlemass','m0'}
            m0 = val;
        case {'velocity','u','uin'}
            U = double(val(:)).';
            if numel(U) < 2, error('velocity must have two components.'); end
            U = U(1:2);
        case 'kbt'
            kBT = val;
        case 'preferlatent'
            preferLatent = logical(val);
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

stateOut = ensure_particle_role(state);
roleInfo = resamp_particle_role_mask(stateOut);
latentSlots = find(roleInfo.latentMask);
freeSlots = find(~stateOut.active(:) & roleInfo.role ~= 1 & ~roleInfo.latentMask);
if preferLatent
    availableSlots = [latentSlots(:); freeSlots(:)];
else
    availableSlots = [freeSlots(:); latentSlots(:)];
end

Nx = params.Nx; Ny = params.Ny;
dx = params.Lx / Nx; dy = params.Ly / Ny;
activatedPerCell = zeros(Nx*Ny, 1);
activatedIds = [];
slotPtr = 1;
capacityHit = false;

for jj = 1:numel(cellIds)
    c = cellIds(jj);
    nWant = max(0, round(nPerCell));
    if nWant == 0, continue; end
    nFree = numel(availableSlots) - slotPtr + 1;
    if nFree <= 0
        capacityHit = true;
        break;
    end
    nAdd = min(nWant, nFree);
    if nAdd < nWant
        capacityHit = true;
    end
    ids = availableSlots(slotPtr:(slotPtr+nAdd-1));
    slotPtr = slotPtr + nAdd;
    [ix, iy] = cell_index_from_id(c, Ny);
    stateOut.x(ids,1) = (ix - 1 + rand(nAdd,1)) * dx;
    stateOut.x(ids,2) = (iy - 1 + rand(nAdd,1)) * dy;
    fluct = zeros(nAdd,2);
    if kBT > 0
        fluct = sqrt(kBT) * randn(nAdd,2);
        fluct(:,1) = fluct(:,1) - mean(fluct(:,1), 'omitnan');
        fluct(:,2) = fluct(:,2) - mean(fluct(:,2), 'omitnan');
    end
    stateOut.v(ids,1) = U(1) + fluct(:,1);
    stateOut.v(ids,2) = U(2) + fluct(:,2);
    stateOut.m(ids) = m0;
    stateOut.active(ids) = true;
    stateOut.particleRole(ids) = 1;
    activatedPerCell(c) = activatedPerCell(c) + nAdd;
    activatedIds = [activatedIds; ids(:)]; %#ok<AGROW>
    if capacityHit, break; end
end
roleAfter = resamp_particle_role_mask(stateOut);
stateOut.Nactive = roleAfter.nFluid;
stateOut.Ncapacity = size(stateOut.x, 1);

% Wrap injected particle positions defensively.
if any(stateOut.active)
    fluid = roleAfter.fluidMask;
    if periodicX
        stateOut.x(fluid,1) = mod(stateOut.x(fluid,1), params.Lx);
    end
    if periodicY
        stateOut.x(fluid,2) = mod(stateOut.x(fluid,2), params.Ly);
    end
end

diag = struct();
diag.kind = 'activate_particles_in_cell';
diag.nRequestedCells = numel(cellIds);
diag.nPerCell = nPerCell;
diag.nActivatedParticles = numel(activatedIds);
diag.nCellsActivated = nnz(activatedPerCell > 0);
diag.capacityHit = capacityHit;
diag.activatedPerCellGrid = reshape(activatedPerCell, [Ny, Nx]).';
diag.activatedIds = activatedIds;
diag.nFluidAfter = roleAfter.nFluid;
diag.nLatentAfter = roleAfter.nLatent;
diag.nInactiveAfter = roleAfter.nInactive;
end

function state = ensure_particle_role(state)
Np = size(state.x,1);
if ~isfield(state,'active') || isempty(state.active)
    state.active = true(Np,1);
else
    state.active = logical(state.active(:));
end
if ~isfield(state,'particleRole') || isempty(state.particleRole)
    state.particleRole = zeros(Np,1);
    state.particleRole(state.active) = 1;
else
    state.particleRole = round(double(state.particleRole(:)));
end
end

function [ix, iy] = cell_index_from_id(c, Ny)
ix = floor((c - 1) / Ny) + 1;
iy = c - Ny * (ix - 1);
end

function v = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end
