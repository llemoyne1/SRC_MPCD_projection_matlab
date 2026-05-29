function [stateOut, info] = resamp_set_dry_cell_particles_latent(state, params, varargin)
%RESAMP_SET_DRY_CELL_PARTICLES_LATENT Convert particles in dry cells to latent slots.
%
% This is the particle-level complement of the cell wet/dry mask.  Particles
% located in non-wet cells are removed from the fluid dynamics by setting
% active=false, particleRole=2, m=0 and v=0.  Positions are kept, so the slots
% can still be interpreted as latent markers of the dry region.

if nargin < 2
    error('Usage: [stateOut, info] = resamp_set_dry_cell_particles_latent(state, params, ...)');
end
cellWetMask = [];
periodicX = true;
periodicY = true;
zeroVelocity = true;
zeroMass = true;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'cellwetmask','wetmask','fluidmask'}
            cellWetMask = logical(val);
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case 'zerovelocity'
            zeroVelocity = logical(val);
        case 'zeromass'
            zeroMass = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

stateOut = ensure_particle_role(state);
roleBefore = resamp_particle_role_mask(stateOut);
fluidBefore = roleBefore.fluidMask;
if isempty(cellWetMask)
    [cellWetMask, wetInfo] = resamp_cell_wet_mask(stateOut, params, 'mode', 'auto');
else
    [cellWetMask, wetInfo] = resamp_cell_wet_mask(stateOut, params, 'mode', 'explicit', 'cellWetMask', cellWetMask);
end

latentIds = [];
if any(fluidBefore)
    cellId = resamp_cell_ids_periodic(stateOut.x(fluidBefore,:), params, 'periodicX', periodicX, 'periodicY', periodicY);
    activeIds = find(fluidBefore);
    wetVec = reshape(cellWetMask.', [params.Nx*params.Ny, 1]);
    toLatentLocal = ~wetVec(cellId);
    latentIds = activeIds(toLatentLocal);
    if ~isempty(latentIds)
        stateOut.active(latentIds) = false;
        stateOut.particleRole(latentIds) = 2;
        if zeroMass
            stateOut.m(latentIds) = 0;
        end
        if zeroVelocity
            stateOut.v(latentIds,:) = 0;
        end
    end
end

roleAfter = resamp_particle_role_mask(stateOut);
stateOut.Nactive = roleAfter.nFluid;
stateOut.Ncapacity = size(stateOut.x, 1);

info = struct();
info.kind = 'set_dry_cell_particles_latent';
info.nWetCells = wetInfo.nWetCells;
info.nDryCells = wetInfo.nDryCells;
info.nConvertedToLatent = numel(latentIds);
info.nFluidBefore = roleBefore.nFluid;
info.nFluidAfter = roleAfter.nFluid;
info.nLatentBefore = roleBefore.nLatent;
info.nLatentAfter = roleAfter.nLatent;
info.convertedIds = latentIds(:);
end

function state = ensure_particle_role(state)
Np = size(state.x, 1);
if ~isfield(state, 'active') || isempty(state.active)
    state.active = true(Np, 1);
else
    state.active = logical(state.active(:));
end
if ~isfield(state, 'particleRole') || isempty(state.particleRole)
    state.particleRole = zeros(Np, 1);
    state.particleRole(state.active) = 1;
else
    state.particleRole = round(double(state.particleRole(:)));
end
end
