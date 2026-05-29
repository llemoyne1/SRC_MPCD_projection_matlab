function roleInfo = resamp_particle_role_mask(state)
%RESAMP_PARTICLE_ROLE_MASK Classify particle storage slots by role.
%
% Roles are intentionally minimal and compatible with existing pool states:
%   0 inactive storage/free pool
%   1 fluid particle: participates in deposit/collision/Q6/remap/thermostat
%   2 latent particle: stored marker, non-fluid, m=0/v=0 by convention
%
% Old states without state.particleRole are interpreted from state.active:
% active -> fluid, inactive -> pool.

if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must contain x, v and m.');
end
Np = size(state.x, 1);
if isfield(state, 'active') && ~isempty(state.active)
    storageActive = logical(state.active(:));
    if numel(storageActive) ~= Np
        error('state.active must have one entry per particle slot.');
    end
else
    storageActive = true(Np, 1);
end

if isfield(state, 'particleRole') && ~isempty(state.particleRole)
    role = double(state.particleRole(:));
    if numel(role) ~= Np
        error('state.particleRole must have one entry per particle slot.');
    end
    role(~isfinite(role)) = 0;
    role = round(role);
else
    role = zeros(Np, 1);
    role(storageActive) = 1;
end

finiteState = all(isfinite(state.x), 2) & all(isfinite(state.v), 2) & isfinite(state.m(:));
nonNegativeMass = state.m(:) >= 0;
fluidMask = storageActive & role == 1 & finiteState & nonNegativeMass;
latentMask = role == 2 & finiteState;
inactiveMask = ~fluidMask & ~latentMask;

roleInfo = struct();
roleInfo.role = role;
roleInfo.storageActive = storageActive;
roleInfo.fluidMask = fluidMask;
roleInfo.latentMask = latentMask;
roleInfo.inactiveMask = inactiveMask;
roleInfo.nFluid = nnz(fluidMask);
roleInfo.nLatent = nnz(latentMask);
roleInfo.nInactive = nnz(inactiveMask);
roleInfo.nSlots = Np;
end
