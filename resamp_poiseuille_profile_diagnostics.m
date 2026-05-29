function prof = resamp_poiseuille_profile_diagnostics(state, params)
%RESAMP_POISEUILLE_PROFILE_DIAGNOSTICS Weighted y-profiles for channel runs.
active=resamp_active_mask(state);
G=resamp_deposit_weighted_to_grid(state.x,state.v,state.m,params,'periodicX',true,'periodicY',false,'activeMask',active,'cellWetMask',get_cell_wet_mask(state,params));
Ny=params.Ny; dy=params.Ly/Ny;
prof=struct();
prof.yCenters=((0:Ny-1)+0.5)'*dy;
prof.Ux=mean(G.Ux,1,'omitnan').';
prof.Uy=mean(G.Uy,1,'omitnan').';
prof.N=mean(G.N,1,'omitnan').';
prof.M=mean(G.M,1,'omitnan').';
prof.G=G;
prof.centerVelocity=interp1(prof.yCenters,prof.Ux,params.Ly/2,'linear','extrap');
prof.wallMeanVelocity=mean([prof.Ux(1),prof.Ux(end)],'omitnan');
prof.centerMinusWall=prof.centerVelocity-prof.wallMeanVelocity;
end
function mask=get_cell_wet_mask(state,params)
[mask,~]=resamp_cell_wet_mask(state,params,'mode','auto');
end
