function activeMask = resamp_active_mask(state)
%RESAMP_ACTIVE_MASK Return the fluid-active particle mask for weighted resampling.
%
% Historical states only contain x/v/m and are interpreted as fully fluid
% active.  Pool states add state.active.  Latent-aware states may also add
% state.particleRole with roles:
%   0 inactive/free pool, 1 fluid, 2 latent/non-fluid marker.
%
% Only role==1 particles enter deposits, collisions, projection, thermostat
% and mass remapping.  Latent particles are stored but hydrodynamically inert.

if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must contain x, v and m.');
end
Np = size(state.x, 1);
if size(state.v,1) ~= Np || numel(state.m) ~= Np
    error('state.x, state.v and state.m must have matching particle count.');
end

if isfield(state, 'active') && ~isempty(state.active)
    activeMask = logical(state.active(:));
    if numel(activeMask) ~= Np
        error('state.active must have one entry per particle slot.');
    end
else
    activeMask = true(Np, 1);
end

if isfield(state, 'particleRole') && ~isempty(state.particleRole)
    role = double(state.particleRole(:));
    if numel(role) ~= Np
        error('state.particleRole must have one entry per particle slot.');
    end
    activeMask = activeMask & round(role) == 1;
end

finiteState = all(isfinite(state.x), 2) & all(isfinite(state.v), 2) & isfinite(state.m(:));
nonNegativeMass = state.m(:) >= 0;
activeMask = activeMask & finiteState & nonNegativeMass;
end
