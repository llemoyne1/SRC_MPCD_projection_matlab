function [stateOut, info] = resamp_enable_particle_pool(state, varargin)
%RESAMP_ENABLE_PARTICLE_POOL Preallocate inactive particle slots.
%
%   [stateOut, info] = resamp_enable_particle_pool(state, 'capacityFactor', 2)
%
% The active particles are copied to the beginning of the arrays.  Additional
% rows are inactive slots with m=0.  All later resamp_* routines in this
% patch ignore inactive slots through state.active.

if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must contain x, v and m.');
end
capacityFactor = 2.0;
minExtraSlots = 0;
explicitCapacity = [];
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'capacityfactor','factor'}
            capacityFactor = val;
        case {'minextraslots','minfree'}
            minExtraSlots = val;
        case {'ncapacity','capacity'}
            explicitCapacity = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

activeMask = resamp_active_mask(state);
xActive = state.x(activeMask, :);
vActive = state.v(activeMask, :);
mActive = state.m(activeMask);
Nactive = size(xActive, 1);

if isempty(explicitCapacity)
    Ncapacity = max(Nactive, ceil(capacityFactor * max(Nactive, 1)));
    Ncapacity = max(Ncapacity, Nactive + max(0, round(minExtraSlots)));
else
    Ncapacity = round(explicitCapacity);
    if Ncapacity < Nactive
        error('Requested capacity %d is smaller than active count %d.', Ncapacity, Nactive);
    end
end

stateOut = state;
stateOut.x = zeros(Ncapacity, 2);
stateOut.v = zeros(Ncapacity, 2);
stateOut.m = zeros(Ncapacity, 1);
stateOut.active = false(Ncapacity, 1);
stateOut.particleRole = zeros(Ncapacity, 1);
if Nactive > 0
    stateOut.x(1:Nactive,:) = xActive;
    stateOut.v(1:Nactive,:) = vActive;
    stateOut.m(1:Nactive) = mActive(:);
    stateOut.active(1:Nactive) = true;
    stateOut.particleRole(1:Nactive) = 1;
end
stateOut.Nactive = Nactive;
stateOut.Ncapacity = Ncapacity;
stateOut.Ninitial = Nactive;

info = resamp_particle_pool_info(stateOut);
info.capacityFactorRequested = capacityFactor;
info.minExtraSlotsRequested = minExtraSlots;
end
