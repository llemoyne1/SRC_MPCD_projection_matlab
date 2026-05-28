function info = resamp_particle_pool_info(state)
%RESAMP_PARTICLE_POOL_INFO Summarize active/free capacity of a particle pool state.

activeMask = resamp_active_mask(state);
Ncapacity = size(state.x, 1);
Nactive = nnz(activeMask);
info = struct();
info.hasPool = isfield(state, 'active');
info.Ncapacity = Ncapacity;
info.Nactive = Nactive;
info.Nfree = Ncapacity - Nactive;
info.activeFraction = Nactive / max(Ncapacity, 1);
info.freeFraction = info.Nfree / max(Ncapacity, 1);
if isfield(state, 'Ninitial') && ~isempty(state.Ninitial)
    info.Ninitial = state.Ninitial;
else
    info.Ninitial = Nactive;
end
if isfield(state, 'Ncapacity') && ~isempty(state.Ncapacity)
    info.NcapacityDeclared = state.Ncapacity;
else
    info.NcapacityDeclared = Ncapacity;
end
end
