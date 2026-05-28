function activeMask = resamp_active_mask(state)
%RESAMP_ACTIVE_MASK Return the active-particle mask for weighted resampling states.
%
% Historical weighted states only contain x/v/m and are interpreted as fully
% active.  Pool states add state.active, state.Nactive and state.Ncapacity;
% inactive slots are storage only and must not enter deposits, collisions,
% projection, thermostat or diagnostics.

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

finiteState = all(isfinite(state.x), 2) & all(isfinite(state.v), 2) & isfinite(state.m(:));
nonNegativeMass = state.m(:) >= 0;
activeMask = activeMask & finiteState & nonNegativeMass;
end
