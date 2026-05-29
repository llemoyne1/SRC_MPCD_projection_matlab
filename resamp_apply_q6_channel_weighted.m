function [stateOut, diag] = resamp_apply_q6_channel_weighted(stateClassic, params)
%RESAMP_APPLY_Q6_CHANNEL_WEIGHTED Apply periodic-x / bounded-y Q6 projection.

activeMask=resamp_active_mask(stateClassic);
cellWetMask=[];
if isfield(stateClassic,'cellWetMask') && ~isempty(stateClassic.cellWetMask)
    cellWetMask=stateClassic.cellWetMask;
end
Gbefore=resamp_deposit_weighted_to_grid(stateClassic.x,stateClassic.v,stateClassic.m,params,'periodicX',true,'periodicY',false,'minMass',eps,'activeMask',activeMask,'cellWetMask',cellWetMask);
proj=projection_project_grid_periodic_x_neumann_y(Gbefore.Ux,Gbefore.Uy,params);
stateOut=stateClassic;
dvApplied=zeros(size(stateClassic.v));
projectionStrength=get_param(params,'projectionStrength',1.0);
if logical(get_param(params,'projectionEnable',true)) && projectionStrength~=0 && any(activeMask)
    dv=projection_interpolate_grid_delta_to_particles(stateClassic.x(activeMask,:),proj.dUx,proj.dUy,params,'periodicX',true,'periodicY',false,'method',get_param(params,'projectionInterpolationMethod','nearest'));
    dvApplied(activeMask,:)=projectionStrength*dv;
    [stateOut.v,momCorr]=resamp_apply_global_momentum_correction_weighted(stateClassic.v,dvApplied,stateClassic.m,params,'stage','weighted_q6_channel','activeMask',activeMask);
else
    [stateOut.v,momCorr]=resamp_apply_global_momentum_correction_weighted(stateClassic.v,dvApplied,stateClassic.m,params,'stage','weighted_q6_channel_disabled','activeMask',activeMask);
end
stateOut.Nactive=nnz(resamp_active_mask(stateOut)); stateOut.Ncapacity=size(stateOut.x,1);
if logical(get_param(params,'thermostatAfterProjection',false))
    [stateOut.v,thermostatInfo]=resamp_apply_cell_thermostat_weighted(stateOut.x,stateOut.v,stateOut.m,params,'periodicX',true,'periodicY',false,'activeMask',resamp_active_mask(stateOut));
else
    thermostatInfo=struct('enabled',false,'meanKBTBefore',NaN,'meanKBTAfter',NaN);
end
Gafter=resamp_deposit_weighted_to_grid(stateOut.x,stateOut.v,stateOut.m,params,'periodicX',true,'periodicY',false,'minMass',eps,'activeMask',resamp_active_mask(stateOut),'cellWetMask',cellWetMask);
projAfter=projection_project_grid_periodic_x_neumann_y(Gafter.Ux,Gafter.Uy,params);
md=resamp_population_mass_diagnostics(stateOut,params,'periodicX',true,'periodicY',false,'cellWetMask',Gafter.cellWetMask);

diag=struct(); diag.kind='weighted_q6_channel'; diag.proj=proj; diag.projAfterParticles=projAfter; diag.thermostat=thermostatInfo;
diag.rmsDivBefore=proj.rmsDivBefore; diag.rmsDivAfterGrid=proj.rmsDivAfter; diag.rmsDivParticleAfter=projAfter.rmsDivBefore;
diag.divReductionGrid=proj.rmsDivAfter/max(proj.rmsDivBefore,eps); diag.divReductionParticle=diag.rmsDivParticleAfter/max(proj.rmsDivBefore,eps);
diag.dvAppliedRms=sqrt(mean(sum(dvApplied(activeMask,:).^2,2),'omitnan')); diag.momentumCorrection=momCorr;
diag.NpActive=md.NpActive; diag.Ncapacity=md.Ncapacity; diag.Nfree=md.Nfree; diag.NStdAfterProjection=md.NStd; diag.MRelRmsAfterProjection=md.MRelRms; diag.kBTWeighted=md.kBTWeighted; diag.massDiagnosticsAfterProjection=md;
end
function v=get_param(params,name,defaultValue)
if isfield(params,name)&&~isempty(params.(name)), v=params.(name); else, v=defaultValue; end
end
