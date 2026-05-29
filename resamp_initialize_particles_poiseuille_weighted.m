function [state, info] = resamp_initialize_particles_poiseuille_weighted(params, varargin)
%RESAMP_INITIALIZE_PARTICLES_POISEUILLE_WEIGHTED Initialize weighted channel particles.
%
% Periodic in x, bounded in y.  The initial state is fully wet and contains
% state.x, state.v, state.m with m(:)=resampParticleMass by default.

mode = get_param(params, 'initialPopulationMode', 'exact_per_cell');
zeroGlobalMean = logical(get_param(params, 'initialVelocityZeroGlobalMean', true));
initialMeanVelocityX = get_param(params, 'initialMeanVelocityX', 0.0);
initialMeanVelocityY = get_param(params, 'initialMeanVelocityY', 0.0);
profileEnable = logical(get_param(params, 'initialPoiseuilleProfileEnable', false));
nuGuess = get_param(params, 'initialPoiseuilleNuGuess', 0.05);
profileScale = get_param(params, 'initialPoiseuilleScale', 1.0);
profileSubtractMean = logical(get_param(params, 'initialPoiseuilleSubtractMean', false));
for k=1:2:numel(varargin)
    key=lower(char(string(varargin{k}))); val=varargin{k+1};
    switch key
        case 'initialpopulationmode', mode=char(string(val));
        case 'initialvelocityzeroglobalmean', zeroGlobalMean=logical(val);
        case 'initialmeanvelocityx', initialMeanVelocityX=val;
        case 'initialmeanvelocityy', initialMeanVelocityY=val;
        case 'initialpoiseuilleprofileenable', profileEnable=logical(val);
        case 'initialpoiseuillenuguess', nuGuess=val;
        case 'initialpoiseuillescale', profileScale=val;
        case 'initialpoiseuillesubtractmean', profileSubtractMean=logical(val);
        otherwise, error('Unknown option: %s', key);
    end
end

Nx=params.Nx; Ny=params.Ny; Lx=params.Lx; Ly=params.Ly; gamma=params.gamma; kBT=params.kBT;
dx=Lx/Nx; dy=Ly/Ny;
mode = lower(strrep(char(mode),'-','_'));
switch mode
    case {'exact_per_cell','cell_exact','uniform_by_cell'}
        gammaInt=round(gamma);
        if abs(gamma-gammaInt)>1e-12, error('exact_per_cell requires integer gamma.'); end
        Np=Nx*Ny*gammaInt; x=zeros(Np,2); p=0;
        for ix=1:Nx
            for iy=1:Ny
                ids=p+(1:gammaInt);
                x(ids,1)=(ix-1+rand(gammaInt,1))*dx;
                x(ids,2)=(iy-1+rand(gammaInt,1))*dy;
                p=p+gammaInt;
            end
        end
    case {'random','random_uniform','uniform'}
        Np=round(Nx*Ny*gamma); x=[Lx*rand(Np,1), Ly*rand(Np,1)];
    otherwise
        error('Unknown initialPopulationMode: %s', mode);
end

v=sqrt(max(kBT,0))*randn(Np,2);
if zeroGlobalMean && Np>0
    v(:,1)=v(:,1)-mean(v(:,1),'omitnan');
    v(:,2)=v(:,2)-mean(v(:,2),'omitnan');
end
v(:,1)=v(:,1)+initialMeanVelocityX;
v(:,2)=v(:,2)+initialMeanVelocityY;

profileInfo=struct('enabled',profileEnable,'nuGuess',nuGuess,'scale',profileScale,'subtractMean',profileSubtractMean,'UmaxNoSlipFormula',NaN,'meanProfileAdded',0,'centerMinusWallProfile',0);
if profileEnable && Np>0
    bodyForceX=get_param(params,'bodyForceX',0.0);
    Ubottom=get_param(params,'Ubottom',0.0); Utop=get_param(params,'Utop',0.0);
    y=x(:,2);
    prof = Ubottom + (Utop-Ubottom)*y/Ly + bodyForceX/(2*max(nuGuess,eps))*y.*(Ly-y);
    profileInfo.UmaxNoSlipFormula = bodyForceX*Ly^2/(8*max(nuGuess,eps));
    if profileSubtractMean
        prof = prof - mean(prof,'omitnan');
    end
    prof = profileScale * prof;
    v(:,1)=v(:,1)+prof;
    profileInfo.meanProfileAdded=mean(prof,'omitnan');
    profileInfo.centerMinusWallProfile=max(prof)-mean([prof(y<dy); prof(y>Ly-dy)],'omitnan');
end

m0=get_param(params,'resampParticleMass',get_param(params,'particleMass',1.0));
m=m0*ones(Np,1);
state=struct('x',x,'v',v,'m',m);
state.cellWetMask=true(Nx,Ny);
G=resamp_deposit_weighted_to_grid(state.x,state.v,state.m,params,'periodicX',true,'periodicY',false,'activeMask',true(Np,1),'cellWetMask',state.cellWetMask);
md=resamp_population_mass_diagnostics(state,params,'periodicX',true,'periodicY',false,'cellWetMask',state.cellWetMask);
info=struct('Np',Np,'mode',mode,'particleMass',m0,'nominalCellMass',gamma*m0,'profile',profileInfo,'initialDiagnostics',md,'G',G);
end

function v=get_param(params,name,defaultValue)
if isstruct(params)&&isfield(params,name)&&~isempty(params.(name)), v=params.(name); else, v=defaultValue; end
end
