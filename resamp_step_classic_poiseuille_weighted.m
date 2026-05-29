function [stateOut, diag] = resamp_step_classic_poiseuille_weighted(state, params)
%RESAMP_STEP_CLASSIC_POISEUILLE_WEIGHTED Weighted SRC/MPCD step in a channel.
%
% Periodic x, bounded y.  When wallVirtualParticlesEnable=true, the SRD
% collision uses aggregate virtual wall particles, in weighted form, to
% restore near-wall collision statistics and reduce slip.  Virtual particles
% are not transported, stored in the pool, remapped, or counted as fluid.

validate_state(state);
Lx=params.Lx; Ly=params.Ly; Nx=params.Nx; Ny=params.Ny; dt=params.dt;
bodyForceX=get_param(params,'bodyForceX',0.0); bodyForceY=get_param(params,'bodyForceY',0.0);
alphaDeg=get_param(params,'alphaDeg',90);

activeMask=resamp_active_mask(state); idx=find(activeMask);
stateOut=state;
if isempty(idx)
    diag=empty_diag(state,params); return;
end
x=state.x(idx,:); v=state.v(idx,:); m=state.m(idx);
dx=Lx/Nx; dy=Ly/Ny;

% Force kick and streaming.
v(:,1)=v(:,1)+dt*bodyForceX;
v(:,2)=v(:,2)+dt*bodyForceY;
x(:,1)=mod(x(:,1)+dt*v(:,1), Lx);
x(:,2)=x(:,2)+dt*v(:,2);
[x,v,wallInfo]=mpcd_apply_wall_bc_y(x,v,params);

Pbefore=[sum(m.*v(:,1),'omitnan'), sum(m.*v(:,2),'omitnan')];

% Weighted SRD collision.  This routine also handles the no-wall-VP case,
% so all channel collision paths use exactly the same binning semantics.
[v, collisionInfo] = resamp_srd_collision_channel_virtual_walls_weighted(x, v, m, params);
Pafter=[sum(m.*v(:,1),'omitnan'), sum(m.*v(:,2),'omitnan')];

stateOut.x(idx,:)=x; stateOut.v(idx,:)=v; stateOut.m(idx)=m;
stateOut.Nactive=nnz(resamp_active_mask(stateOut)); stateOut.Ncapacity=size(stateOut.x,1);

if logical(get_param(params,'thermostatAfterStep',false))
    [stateOut.v, thermostatInfo]=resamp_apply_cell_thermostat_weighted(stateOut.x,stateOut.v,stateOut.m,params,'periodicX',true,'periodicY',false,'activeMask',resamp_active_mask(stateOut));
else
    thermostatInfo=empty_thermostat_info();
end
md=resamp_population_mass_diagnostics(stateOut,params,'periodicX',true,'periodicY',false,'cellWetMask',get_cell_wet_mask(stateOut,params));

wvp = collisionInfo.wallVirtualParticles;

diag=struct();
diag.kind='weighted_classic_poiseuille'; diag.Np=numel(idx); diag.NpActive=md.NpActive; diag.Ncapacity=md.Ncapacity; diag.Nfree=md.Nfree;
diag.Nx=Nx; diag.Ny=Ny; diag.dx=dx; diag.dy=dy; diag.shiftX=collisionInfo.shiftX; diag.shiftY=collisionInfo.shiftY; diag.alphaDeg=alphaDeg;
diag.bodyForceX=bodyForceX; diag.bodyForceY=bodyForceY; diag.wallInfo=wallInfo; diag.thermostat=thermostatInfo;
diag.NMean=md.NMean; diag.NStd=md.NStd; diag.NMin=md.NMin; diag.NMax=md.NMax; diag.MRelRms=md.MRelRms;
diag.meanVxWeighted=md.meanVxWeighted; diag.meanVyWeighted=md.meanVyWeighted; diag.kBTWeighted=md.kBTWeighted;
diag.totalMass=md.totalMass; diag.totalMomentum=md.totalMomentum; diag.momentumBeforeCollision=Pbefore; diag.momentumAfterCollision=Pafter; diag.collisionDeltaPNorm=norm(Pafter-Pbefore);
diag.massDiagnostics=md;
diag.NMeanCollision=collisionInfo.NMeanReal; diag.NStdCollision=collisionInfo.NStdReal;
diag.MMeanCollision=collisionInfo.MMeanReal; diag.MStdCollision=collisionInfo.MStdReal;
diag.wallVirtualParticles=wvp;
diag.wallVirtualParticlesEnable=logical(get_param(params,'wallVirtualParticlesEnable',false));
diag.wallVirtualParticlesGeometryMode=collisionInfo.geometryMode;
diag.nVirtualWallCells=collisionInfo.nVirtualCells;
diag.nVirtualWallParticlesEquivalent=collisionInfo.nVirtualParticlesEquivalentTotal;
diag.virtualWallMassTotal=collisionInfo.virtualMassTotal;
diag.virtualWallMomentumXTotal=collisionInfo.virtualMomentumXTotal;
diag.virtualWallMomentumYTotal=collisionInfo.virtualMomentumYTotal;
diag.wallVPBottomParticles=get_field(wvp,'nBottomVirtualParticles',0);
diag.wallVPTopParticles=get_field(wvp,'nTopVirtualParticles',0);
diag.wallVPBottomCells=get_field(wvp,'nBottomCells',0);
diag.wallVPTopCells=get_field(wvp,'nTopCells',0);
end

function mask=get_cell_wet_mask(state,params)
[mask,~]=resamp_cell_wet_mask(state,params,'mode','auto');
end
function diag=empty_diag(state,params)
md=resamp_population_mass_diagnostics(state,params,'periodicX',true,'periodicY',false);
diag=struct('kind','weighted_classic_poiseuille_empty','Np',0,'NpActive',md.NpActive,'Ncapacity',md.Ncapacity,'Nfree',md.Nfree,'NMean',md.NMean,'NStd',md.NStd,'NMin',md.NMin,'NMax',md.NMax,'MRelRms',md.MRelRms,'kBTWeighted',md.kBTWeighted);
end
function validate_state(state)
if ~isstruct(state)||~isfield(state,'x')||~isfield(state,'v')||~isfield(state,'m'), error('state must contain x, v, m.'); end
end
function v=get_param(params,name,defaultValue)
if isfield(params,name)&&~isempty(params.(name)), v=params.(name); else, v=defaultValue; end
end
function v=get_field(s,name,defaultValue)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=defaultValue; end
end
function info=empty_thermostat_info()
info=struct('enabled',false,'targetKBT',NaN,'strength',NaN,'nThermostattedCells',0,'meanKBTBefore',NaN,'meanKBTAfter',NaN);
end
