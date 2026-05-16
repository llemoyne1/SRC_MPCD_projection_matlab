function out = run_projection_forced_taylor_green_demo(params)
%RUN_PROJECTION_FORCED_TAYLOR_GREEN_DEMO Forced periodic Taylor-Green benchmark.
%
% Methods:
%   classic : SRC/MPCD + TG forcing, no projection
%   q6      : classic + div(u)=0 projection
%   q9      : q6 + low-k div(Nu) relaxation

if nargin < 1 || isempty(params), params = struct(); end
params = set_default_params(params);
rng(params.seed);
[state, initInfo] = projection_initialize_particles_taylor_green_forced(params);
Np = initInfo.Np;

nSamples = floor(params.nSteps/params.sampleEvery) + 1;
series = allocate_series(nSamples);
meanAccum = init_mean_accumulator(params);
isamp = 0; actualLastStep = 0; stopReason = '';
lastStepDiag = [];
ticRun = tic; lastProgress = -inf;

for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        [sampleDiag, G] = tg_sample_diagnostics(state, params, it);
        series = store_sample(series, isamp, it, sampleDiag, lastStepDiag, params);
        if logical(params.meanFieldEnable) && it >= params.meanFieldStartStep
            meanAccum = update_mean_accumulator(meanAccum, G);
        end
    end

    if logical(params.visualEnable) && (mod(it, params.visualEvery)==0 || it==params.nSteps)
        projection_taylor_green_visualize_frame(state, params, it, ...
            'saveFrame', logical(params.visualSaveFrames), ...
            'frameDir', params.visualFrameDir, ...
            'framePrefix', params.visualFramePrefix);
        if logical(params.meanFieldEnable) && logical(params.meanFieldVisualEnable) && meanAccum.n > 0
            Gmean = mean_accumulator_to_grid(meanAccum, params);
            pMean = params;
            pMean.visualFigureId = params.meanFieldFigureId;
            pMean.visualFigureName = params.meanFieldFigureName;
            pMean.visualTitleSuffix = sprintf('n=%d start=%g', meanAccum.n, params.meanFieldStartStep*params.dt);
            projection_taylor_green_visualize_frame(Gmean, pMean, it, 'isMean', true, ...
                'saveFrame', logical(params.meanFieldSaveFrames), ...
                'frameDir', params.meanFieldFrameDir, ...
                'framePrefix', params.meanFieldFramePrefix);
        end
        pause(params.visualPause);
    end

    if logical(params.abortOnUnstable) && it > 0 && mod(it, params.sampleEvery)==0
        [unstable, reason] = check_unstable(series, isamp, params);
        if unstable
            stopReason = reason;
            actualLastStep = it;
            warning('run_projection_forced_taylor_green_demo:Unstable', ...
                'Stopping early at step %d/%d: %s', it, params.nSteps, reason);
            break;
        end
    end

    if it == params.nSteps
        actualLastStep = it;
        break;
    end

    stepParams = params;
    stepParams.computeDiagnostics = logical(params.computeFullDiagnosticsEveryStep) || mod(it+1, params.sampleEvery)==0 || mod(it+1, params.progressEvery)==0;
    switch lower(params.method)
        case 'classic'
            [state, stepDiag] = mpcd_step_classic_periodic_forced(state, stepParams);
        case {'q6','q9'}
            [state, stepDiag] = mpcd_step_projection_periodic_q9_forced(state, stepParams);
        otherwise
            error('Unknown method: %s', params.method);
    end
    lastStepDiag = stepDiag;
    actualLastStep = it + 1;

    if params.progressEvery > 0 && mod(actualLastStep, params.progressEvery)==0 && actualLastStep ~= lastProgress
        lastProgress = actualLastStep;
        fprintf('  forced TG %s step %d/%d, t=%.4g, A=%.5g, coh=%.3g, enst=%.4g, highK=%.3g, lowKdens=%.3e, kBT=%.4g, elapsed=%.1fs\n', ...
            upper(params.method), actualLastStep, params.nSteps, actualLastStep*params.dt, ...
            series.tgAmplitude(isamp), series.tgCoherence(isamp), series.enstrophy(isamp), ...
            series.highKVelocityFraction(isamp), series.lowKDensity(isamp), series.kBTCell(isamp), toc(ticRun));
    end

    if isfinite(params.maxWallClockSeconds) && toc(ticRun) > params.maxWallClockSeconds
        stopReason = 'maxWallClockSeconds';
        warning('run_projection_forced_taylor_green_demo:TimeLimit', ...
            'Stopping early at step %d/%d after %.1fs.', actualLastStep, params.nSteps, toc(ticRun));
        break;
    end
end

elapsed = toc(ticRun);
series = trim_series(series, isamp);
series = add_forced_balance_estimates(series, params);
summary = summarize_run(series, params, actualLastStep, elapsed, stopReason);

out = struct();
out.params = params;
out.initialization = initInfo;
out.state = state;
out.series = series;
out.meanFields = meanAccum;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.stopReason = stopReason;
out.elapsedWallClock = elapsed;
out.summary = summary;

fprintf('\n=== run_projection_forced_taylor_green_demo (%s) ===\n', upper(params.method));
fprintf('Np                         : %d\n', Np);
fprintf('grid/gamma                 : %d x %d / %g\n', params.Nx, params.Ny, params.gamma);
fprintf('steps completed            : %d / %d in %.2f s\n', actualLastStep, params.nSteps, elapsed);
fprintf('TG force/initial amplitude  : %.12g / %.12g\n', params.taylorGreenForceAmplitude, params.taylorGreenInitialAmplitude);
fprintf('mean/final amplitude        : %.12g / %.12g\n', summary.meanAmplitude, summary.finalAmplitude);
fprintf('mean/final coherence        : %.12g / %.12g\n', summary.meanCoherence, summary.finalCoherence);
fprintf('mean/final high-k frac      : %.12g / %.12g\n', summary.meanHighKFraction, summary.finalHighKFraction);
fprintf('mean/final low-k density    : %.12g / %.12g\n', summary.meanLowKDensity, summary.finalLowKDensity);
fprintf('mean/final kBT cell         : %.12g / %.12g\n', summary.meanKBTCell, summary.finalKBTCell);
fprintf('nu_eff fit preferred        : %.12g  (%s, window %.4g--%.4g)\n', ...
    summary.nuEffFitPreferred, summary.viscosityFitMethod, summary.viscosityFitT0, summary.viscosityFitT1);
fprintf('nu_eff derivative/plateau   : %.12g / %.12g\n', summary.nuEffFitDerivativeKnownF, summary.nuEffFitPlateau);
end

function params = set_default_params(params)
params = set_default(params,'method','q9');
params = set_default(params,'Lx',1.0); params = set_default(params,'Ly',1.0);
params = set_default(params,'Nx',64); params = set_default(params,'Ny',64);
params = set_default(params,'gamma',20); params = set_default(params,'seed',11);
params = set_default(params,'dt',1.0e-3); params = set_default(params,'kBT',0.01);
params = set_default(params,'alphaDeg',90); params = set_default(params,'bodyForceX',0); params = set_default(params,'bodyForceY',0);
params = set_default(params,'useRandomGridShift',true);
params = set_default(params,'initialPopulationMode','exact_per_cell');
params = set_default(params,'initialVelocityZeroGlobalMean',true);
params = set_default(params,'taylorGreenInitialAmplitude',0.10);
params = set_default(params,'taylorGreenAmplitude',params.taylorGreenInitialAmplitude);
params = set_default(params,'taylorGreenModeX',1); params = set_default(params,'taylorGreenModeY',1);
params = set_default(params,'taylorGreenThermalNoise',true);
params = set_default(params,'taylorGreenForceEnable',true);
params = set_default(params,'taylorGreenForceAmplitude',0.08);
params = set_default(params,'taylorGreenForceZeroMeanKick',true);
params = set_default(params,'projectionEnable',true);
params = set_default(params,'projectionStrength',1.0);
params = set_default(params,'projectionInterpolationMethod','nearest');
params = set_default(params,'massFluxProjectionMode','relax_to_uniform_lowk');
params = set_default(params,'massFluxProjectionStrength',1.0);
params = set_default(params,'massFluxDensityRelaxationBeta',5e-4);
params = set_default(params,'massFluxApplyAfterVelocityProjection',true);
params = set_default(params,'massFluxLowKMaxIndex',2);
params = set_default(params,'massFluxFinalVelocityProjectionCleanup',true);
params = set_default(params,'massFluxFinalVelocityProjectionStrength',0.5);
params = set_default(params,'projectionMomentumCorrectionEnable',true);
params = set_default(params,'projectionMomentumCorrectionMode','particle_global_exact');
params = set_default(params,'thermostatAfterStep',false);
params = set_default(params,'thermostatAfterProjection',params.thermostatAfterStep);
params = set_default(params,'nSteps',30000); params = set_default(params,'sampleEvery',100); params = set_default(params,'progressEvery',1000);
params = set_default(params,'computeFullDiagnosticsEveryStep',false);
params = set_default(params,'visualEnable',true); params = set_default(params,'visualEvery',500); params = set_default(params,'visualSaveFrames',false); params = set_default(params,'visualPause',0.001);
params = set_default(params,'visualFigureId',410); params = set_default(params,'visualFigureName','Forced Taylor-Green live'); params = set_default(params,'visualTitleSuffix','');
params = set_default(params,'visualFrameDir','tg_frames'); params = set_default(params,'visualFramePrefix','tg'); params = set_default(params,'visualMaxParticles',8000);
params = set_default(params,'visualDensityCLim',0.6); params = set_default(params,'visualOmegaCLim',NaN);
params = set_default(params,'visualQuiverStrideX',max(2,round(params.Nx/18))); params = set_default(params,'visualQuiverStrideY',max(2,round(params.Ny/12))); params = set_default(params,'visualQuiverScale',1.2);
params = set_default(params,'meanFieldEnable',true); params = set_default(params,'meanFieldVisualEnable',true); params = set_default(params,'meanFieldStartStep',round(0.2*params.nSteps));
params = set_default(params,'meanFieldFigureId',params.visualFigureId+20); params = set_default(params,'meanFieldFigureName','Forced Taylor-Green running mean'); params = set_default(params,'meanFieldSaveFrames',false); params = set_default(params,'meanFieldFrameDir','tg_mean_frames'); params = set_default(params,'meanFieldFramePrefix','tg_mean');
params = set_default(params,'abortOnUnstable',true); params = set_default(params,'maxStableKBTCell',0.2); params = set_default(params,'maxStableHighKFraction',0.85); params = set_default(params,'maxStableDensityRelRMS',1.0); params = set_default(params,'maxWallClockSeconds',Inf);
params = set_default(params,'lowKMaxIndex',2); params = set_default(params,'taylorGreenHighKCut',max(4,min(params.Nx,params.Ny)/4));
params = set_default(params,'taylorGreenViscosityFitFraction',0.60);
end

function [diag,G] = tg_sample_diagnostics(state, params, it)
G = projection_deposit_particles_to_grid(state.x, state.v, params, 'periodicX', true, 'periodicY', true, 'minCount', 1);
tg = projection_taylor_green_diagnostics(G, params);
pop = projection_population_diagnostics(state.x, params, 'periodicX', true, 'periodicY', true);
therm = projection_thermal_diagnostics(state.x, state.v, params, 'periodicX', true, 'periodicY', true, 'minCount', 1);
diag = struct(); diag.step = it; diag.t = it*params.dt; diag.G = G; diag.tg = tg; diag.pop = pop; diag.thermal = therm;
diag.meanUx = mean(state.v(:,1),'omitnan'); diag.meanUy = mean(state.v(:,2),'omitnan'); diag.kBTCell = local_get_kBT_cell(therm);
end


function val = local_get_kBT_cell(therm)
%LOCAL_GET_KBT_CELL Compatibility helper for thermal diagnostics.
% Newer semi-parallel scripts use kBTCellRelative; older snapshots used
% kBTCellMean. Keep both names valid so the Taylor--Green patch can be
% applied on either code lineage.
if isfield(therm,'kBTCellRelative')
    val = therm.kBTCellRelative;
elseif isfield(therm,'kBTCellMean')
    val = therm.kBTCellMean;
elseif isfield(therm,'thermalKineticEnergy')
    val = therm.thermalKineticEnergy;
elseif isfield(therm,'localKBTMean')
    val = therm.localKBTMean;
else
    val = NaN;
end
end

function S=allocate_series(n)
fields={'step','t','meanUx','meanUy','kBTCell','tgAmplitude','tgModeEnergy','tgCoherence','tgResidualEnergy','tgTotalHydroEnergy','enstrophy','omegaRms','highKVelocityFraction','lowKVelocityFraction','tgModeVelocityFraction','rmsDivGrid','NStd','NOutBand','lowKDensity','densityRelRMS','massFluxResidual','divParticleAfter','momentumRawDV','momentumResidualDV'};
for i=1:numel(fields), S.(fields{i})=nan(n,1); end
end
function S=store_sample(S,i,it,d,lastStepDiag,params)
S.step(i)=it; S.t(i)=d.t; S.meanUx(i)=d.meanUx; S.meanUy(i)=d.meanUy; S.kBTCell(i)=d.kBTCell;
S.tgAmplitude(i)=d.tg.modeAmplitude; S.tgModeEnergy(i)=d.tg.modeEnergy; S.tgCoherence(i)=d.tg.modeCoherence; S.tgResidualEnergy(i)=d.tg.residualEnergy; S.tgTotalHydroEnergy(i)=d.tg.totalHydroEnergy;
S.enstrophy(i)=d.tg.enstrophy; S.omegaRms(i)=d.tg.omegaRms; S.highKVelocityFraction(i)=d.tg.highKEnergyFractionVelocity; S.lowKVelocityFraction(i)=d.tg.lowKEnergyFractionVelocity; S.tgModeVelocityFraction(i)=d.tg.tgModeEnergyFractionVelocity; S.rmsDivGrid(i)=d.tg.rmsDiv;
S.NStd(i)=d.pop.stdN; S.NOutBand(i)=d.pop.outBandFraction; S.lowKDensity(i)=d.tg.densityLowKEnergy; S.densityRelRMS(i)=d.tg.densityRelRMS;
S.massFluxResidual(i)=get_nested(lastStepDiag,'rmsMassFluxDivResidual',NaN); S.divParticleAfter(i)=get_nested(lastStepDiag,'rmsDivParticleAfter',NaN);
S.momentumRawDV(i)=get_nested(lastStepDiag,'projectionMomentumTotalRawMeanDVNorm',NaN); S.momentumResidualDV(i)=get_nested(lastStepDiag,'projectionMomentumTotalResidualMeanDVNorm',NaN);
end
function S=trim_series(S,n)
f=fieldnames(S); for i=1:numel(f), S.(f{i})=S.(f{i})(1:n); end
end
function S=add_forced_balance_estimates(S,params)
A=S.tgAmplitude; t=S.t; F=params.taylorGreenForceAmplitude; kx=2*pi*params.taylorGreenModeX/params.Lx; ky=2*pi*params.taylorGreenModeY/params.Ly; k2=kx^2+ky^2;
dA=nan(size(A));
if numel(A)>2
    dA(2:end-1)=(A(3:end)-A(1:end-2))./(t(3:end)-t(1:end-2));
    dA(1)=(A(2)-A(1))/max(t(2)-t(1),eps);
    dA(end)=(A(end)-A(end-1))/max(t(end)-t(end-1),eps);
end
S.dAmplitudeDt=dA;
S.nuEffForced=(F-dA)./max(k2*A,eps);
% Robust run-level fit is stored separately in summarize_run().
end
function meanAccum=init_mean_accumulator(params)
meanAccum=struct('n',0,'sumN',zeros(params.Nx,params.Ny),'sumUx',zeros(params.Nx,params.Ny),'sumUy',zeros(params.Nx,params.Ny));
end
function A=update_mean_accumulator(A,G)
A.n=A.n+1; A.sumN=A.sumN+double(G.N); A.sumUx=A.sumUx+double(G.Ux); A.sumUy=A.sumUy+double(G.Uy);
end
function G=mean_accumulator_to_grid(A,params)
G=struct(); G.N=A.sumN/max(A.n,1); G.Ux=A.sumUx/max(A.n,1); G.Uy=A.sumUy/max(A.n,1); G.Px=G.N.*G.Ux; G.Py=G.N.*G.Uy; G.dx=params.Lx/params.Nx; G.dy=params.Ly/params.Ny; G.Nx=params.Nx; G.Ny=params.Ny; G.Lx=params.Lx; G.Ly=params.Ly; G.rho=G.N/(G.dx*G.dy); G.valid=true(params.Nx,params.Ny);
end
function [unstable,reason]=check_unstable(S,i,params)
unstable=false; reasons={};
if S.kBTCell(i)>params.maxStableKBTCell, unstable=true; reasons{end+1}='kBT'; end
if S.highKVelocityFraction(i)>params.maxStableHighKFraction, unstable=true; reasons{end+1}='highK'; end
if S.densityRelRMS(i)>params.maxStableDensityRelRMS, unstable=true; reasons{end+1}='densityRelRMS'; end
reason=strjoin(reasons,', ');
end
function summary=summarize_run(S,params,lastStep,elapsed,stopReason)
valid=isfinite(S.t); tail=valid & S.t >= (0.5*max(S.t(valid)));
summary=struct(); summary.method=params.method; summary.actualLastStep=lastStep; summary.stoppedEarly=lastStep<params.nSteps; summary.stopReason=stopReason; summary.elapsedWallClock=elapsed;
summary.meanAmplitude=mean(S.tgAmplitude(tail),'omitnan'); summary.finalAmplitude=lastfinite(S.tgAmplitude);
summary.meanCoherence=mean(S.tgCoherence(tail),'omitnan'); summary.finalCoherence=lastfinite(S.tgCoherence);
summary.meanModeEnergy=mean(S.tgModeEnergy(tail),'omitnan'); summary.finalModeEnergy=lastfinite(S.tgModeEnergy);
summary.meanEnstrophy=mean(S.enstrophy(tail),'omitnan'); summary.finalEnstrophy=lastfinite(S.enstrophy);
summary.meanHighKFraction=mean(S.highKVelocityFraction(tail),'omitnan'); summary.finalHighKFraction=lastfinite(S.highKVelocityFraction);
summary.meanLowKDensity=mean(S.lowKDensity(tail),'omitnan'); summary.finalLowKDensity=lastfinite(S.lowKDensity);
summary.meanDensityRelRMS=mean(S.densityRelRMS(tail),'omitnan'); summary.finalDensityRelRMS=lastfinite(S.densityRelRMS);
summary.meanKBTCell=mean(S.kBTCell(tail),'omitnan'); summary.finalKBTCell=lastfinite(S.kBTCell);
summary.meanDivGrid=mean(S.rmsDivGrid(tail),'omitnan'); summary.finalDivGrid=lastfinite(S.rmsDivGrid);
summary.meanDivParticleAfter=mean(S.divParticleAfter(tail),'omitnan'); summary.finalDivParticleAfter=lastfinite(S.divParticleAfter);
summary.meanMassFluxResidual=mean(S.massFluxResidual(tail),'omitnan'); summary.finalMassFluxResidual=lastfinite(S.massFluxResidual);
summary.meanMomentumRawDV=mean(S.momentumRawDV(tail),'omitnan'); summary.meanMomentumResidualDV=mean(S.momentumResidualDV(tail),'omitnan');
summary.meanNuEffForced=mean(S.nuEffForced(tail & isfinite(S.nuEffForced)),'omitnan'); summary.finalNuEffForced=lastfinite(S.nuEffForced);
fitFraction = get_nested(params, 'taylorGreenViscosityFitFraction', 0.60);
try
    fit = projection_fit_forced_taylor_green_viscosity(S.t, S.tgAmplitude, params, 'fitFraction', fitFraction);
catch ME
    fit = struct('status','failed','message',ME.message,'nuPreferred',NaN,'preferredMethod','','nuDerivativeKnownF',NaN,'nuDerivativeFreeF',NaN,'forceDerivativeFit',NaN,'nuExpKnownF',NaN,'nuPlateau',NaN,'r2ExpKnownF',NaN,'r2DerivativeKnownF',NaN,'tStart',NaN,'tEnd',NaN,'nFit',NaN);
end
summary.viscosityFit = fit;
summary.nuEffFitPreferred = get_nested(fit,'nuPreferred',NaN);
summary.nuEffFitDerivativeKnownF = get_nested(fit,'nuDerivativeKnownF',NaN);
summary.nuEffFitDerivativeFreeF = get_nested(fit,'nuDerivativeFreeF',NaN);
summary.tgForceFit = get_nested(fit,'forceDerivativeFit',NaN);
summary.nuEffFitExpKnownF = get_nested(fit,'nuExpKnownF',NaN);
summary.nuEffFitPlateau = get_nested(fit,'nuPlateau',NaN);
summary.viscosityFitR2 = get_nested(fit,'r2ExpKnownF',NaN);
summary.viscosityFitMethod = get_nested(fit,'preferredMethod','');
summary.viscosityFitT0 = get_nested(fit,'tStart',NaN);
summary.viscosityFitT1 = get_nested(fit,'tEnd',NaN);
summary.viscosityFitNSamples = get_nested(fit,'nFit',NaN);
end
function x=lastfinite(v), idx=find(isfinite(v),1,'last'); if isempty(idx), x=NaN; else, x=v(idx); end; end
function params=set_default(params,name,value), if ~isfield(params,name)||isempty(params.(name)), params.(name)=value; end; end
function v=get_nested(s,name,defaultValue), if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=defaultValue; end; end
