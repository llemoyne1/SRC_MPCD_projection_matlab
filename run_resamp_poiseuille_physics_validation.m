function out = run_resamp_poiseuille_physics_validation(varargin)
%RUN_RESAMP_POISEUILLE_PHYSICS_VALIDATION Forced channel validation for weighted resampling.
%
% Minimal Poiseuille benchmark: periodic x, bounded y, bounceback walls and
% bodyForceX.  The runner compares a weighted classic reference and the
% weighted Q6/resampled formulation.
%
% Example:
%   out = run_resamp_poiseuille_physics_validation('steps',5000,'Nx',64,'Ny',32);

opts = parse_options(varargin{:});
if ~exist(opts.outputRoot,'dir'), mkdir(opts.outputRoot); end
caseDefs = build_cases(opts);
outs = struct(); rows = repmat(empty_summary_row(), numel(caseDefs), 1);
for ic=1:numel(caseDefs)
    c=caseDefs(ic);
    fprintf('\n=== Weighted Poiseuille physics case: %s ===\n', c.label);
    outCase = run_one_case(c, opts, ic);
    outs.(c.label)=outCase;
    rows(ic)=summarize_case(c,outCase,opts);
end
summary=struct2table(rows);
writetable(summary, fullfile(opts.outputRoot,'resamp_poiseuille_physics_summary.csv'));
out=struct('options',opts,'caseDefs',caseDefs,'cases',outs,'summary',summary,'outputRoot',opts.outputRoot);
save(fullfile(opts.outputRoot,'resamp_poiseuille_physics_validation.mat'),'out','-v7.3');
fprintf('\nWeighted Poiseuille summary:\n'); disp(summary);
fprintf('Wrote %s\n', fullfile(opts.outputRoot,'resamp_poiseuille_physics_summary.csv'));
end

function caseDefs = build_cases(opts)
caseDefs=struct([]);
caseDefs(1).label='classic_reference'; caseDefs(1).method='weighted_classic'; caseDefs(1).projectionStrength=0; caseDefs(1).extractEvery=0; caseDefs(1).insertEvery=0; caseDefs(1).remapEvery=0; caseDefs(1).thermostatAfterStep=true; caseDefs(1).thermostatAfterRemap=false; caseDefs(1).massSafetyEnable=false; caseDefs(1).NMin=0; caseDefs(1).NMax=2*opts.NTarget;
caseDefs(2).label='q6_resampled'; caseDefs(2).method='weighted_q6'; caseDefs(2).projectionStrength=1; caseDefs(2).extractEvery=1; caseDefs(2).insertEvery=1; caseDefs(2).remapEvery=1; caseDefs(2).thermostatAfterStep=false; caseDefs(2).thermostatAfterRemap=true; caseDefs(2).massSafetyEnable=true; caseDefs(2).NMin=opts.NMin; caseDefs(2).NMax=opts.NMax;
end

function outCase = run_one_case(c, opts, ic)
if ~isempty(opts.rngSeed), rng(opts.rngSeed+ic-1,'twister'); end
params = make_params(opts, c);
[state, initInfo] = resamp_initialize_particles_poiseuille_weighted(params);
[state, poolInfo0] = resamp_enable_particle_pool(state,'capacityFactor',opts.capacityFactor);
state.uMemUx=zeros(params.Nx,params.Ny); state.uMemUy=zeros(params.Nx,params.Ny); state.uMemValid=false(params.Nx,params.Ny);
state.cellWetMask=true(params.Nx,params.Ny); params.cellWetMask=state.cellWetMask;

nSamples=floor(opts.steps/opts.sampleEvery)+1;
UxProfiles=nan(params.Ny,nSamples); UyProfiles=nan(params.Ny,nSamples); NProfiles=nan(params.Ny,nSamples); MProfiles=nan(params.Ny,nSamples); sampleTimes=nan(nSamples,1); sampleSteps=nan(nSamples,1);
rows=repmat(empty_timeseries_row(), floor(opts.steps/opts.summaryEvery)+2, 1); irow=0; isamp=0;
insertStats=empty_insert_stats(); extractStats=empty_extract_stats(); lastStepDiag=struct(); lastInsert=struct(); lastExtract=struct(); lastRemap=struct(); lastThermostat=struct();

[prof, md] = sample_state(state, params); isamp=isamp+1; sampleSteps(isamp)=0; sampleTimes(isamp)=0; UxProfiles(:,isamp)=prof.Ux; UyProfiles(:,isamp)=prof.Uy; NProfiles(:,isamp)=prof.N; MProfiles(:,isamp)=prof.M;
irow=irow+1; rows(irow)=make_timeseries_row(0,0,md,prof,struct(),struct(),struct(),struct(),struct());

for step=1:opts.steps
    switch c.method
        case 'weighted_classic'
            [state, stepDiag]=resamp_step_classic_poiseuille_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag]=resamp_step_projection_poiseuille_weighted(state, params);
        otherwise
            error('Unknown method: %s', c.method);
    end
    preEditG=[];
    if opts.preservePreEditVelocity && c.remapEvery>0 && mod(step,c.remapEvery)==0
        preEditG=resamp_deposit_weighted_to_grid(state.x,state.v,state.m,params,'periodicX',true,'periodicY',false,'activeMask',resamp_active_mask(state),'cellWetMask',state.cellWetMask);
    end
    if c.extractEvery>0 && mod(step,c.extractEvery)==0
        [state, extractDiag]=resamp_extract_overpopulated_particles(state,params,'NTarget',opts.NTarget,'NMin',c.NMin,'NMax',c.NMax,'extractSelectionMode',opts.extractSelectionMode,'periodicX',true,'periodicY',false,'cellWetMask',state.cellWetMask,'computeDiagnostics',false);
    else
        extractDiag=empty_extract_diag();
    end
    extractStats=update_extract_stats(extractStats,extractDiag,step); extractDiag=attach_extract_cumulative(extractDiag,extractStats); lastExtract=extractDiag;
    if c.insertEvery>0 && mod(step,c.insertEvery)==0
        [state, insertDiag]=resamp_insert_underpopulated_particles(state,params,'NTarget',opts.NTarget,'NMin',c.NMin,'NMax',c.NMax,'memoryMinParticles',opts.memoryMinParticles,'particleMass',opts.particleMass,'kBT',opts.kBT,'insertVelocityMode',opts.insertVelocityMode,'periodicX',true,'periodicY',false,'cellWetMask',state.cellWetMask,'computeDiagnostics',false);
    else
        insertDiag=empty_insert_diag();
    end
    insertStats=update_insert_stats(insertStats,insertDiag,step); insertDiag=attach_insert_cumulative(insertDiag,insertStats); lastInsert=insertDiag;
    if c.remapEvery>0 && mod(step,c.remapEvery)==0
        remapArgs={'targetCellMass',params.resampTargetCellMass,'method',opts.remapMethod,'massMin',opts.massMinFactor*params.resampParticleMass,'massMax',opts.massMaxFactor*params.resampParticleMass,'periodicX',true,'periodicY',false,'cellWetMask',state.cellWetMask,'constraintTolerance',opts.constraintTolerance,'massSafetyEnable',c.massSafetyEnable,'massSafetyMode',opts.remapMassSafetyMode,'massSafetyMinFactor',opts.massSafetyMinFactor,'massSafetyMaxFactor',opts.massSafetyMaxFactor,'computeDiagnostics',false};
        didEdit=get_field(extractDiag,'nExtractedParticles',0)>0 || get_field(insertDiag,'nInsertedParticles',0)>0;
        if opts.preservePreEditVelocity && didEdit
            remapArgs=[remapArgs, {'targetVelocityMode','grid','targetUx',preEditG.Ux,'targetUy',preEditG.Uy}]; %#ok<AGROW>
        else
            remapArgs=[remapArgs, {'targetVelocityMode','preserve_cell_velocity'}]; %#ok<AGROW>
        end
        [state, remapDiag]=resamp_local_mass_moment_remap(state,params,remapArgs{:});
    else
        remapDiag=empty_remap_diag();
    end
    lastRemap=remapDiag;
    if c.thermostatAfterRemap
        [state.v, thermostatDiag]=resamp_apply_cell_thermostat_weighted(state.x,state.v,state.m,params,'periodicX',true,'periodicY',false,'activeMask',resamp_active_mask(state),'targetKBT',opts.kBT,'strength',opts.thermostatStrength);
    else
        thermostatDiag=struct('enabled',false,'meanKBTBefore',NaN,'meanKBTAfter',NaN);
    end
    lastThermostat=thermostatDiag; lastStepDiag=stepDiag;
    doSample=(mod(step,opts.sampleEvery)==0)||step==opts.steps; doSummary=(mod(step,opts.summaryEvery)==0)||step==opts.steps;
    if doSample || doSummary
        [prof, md] = sample_state(state, params);
    end
    if doSample
        isamp=isamp+1; sampleSteps(isamp)=step; sampleTimes(isamp)=step*params.dt; UxProfiles(:,isamp)=prof.Ux; UyProfiles(:,isamp)=prof.Uy; NProfiles(:,isamp)=prof.N; MProfiles(:,isamp)=prof.M;
    end
    if doSummary
        irow=irow+1; rows(irow)=make_timeseries_row(step,step*params.dt,md,prof,stepDiag,insertDiag,extractDiag,remapDiag,thermostatDiag);
        fprintf('poiseuille %-18s step=%6d t=%g Ucw=%.4g N[%g,%g] Mrel=%.2e kBT=%.4g ins=%g ext=%g\n', c.label, step, step*params.dt, prof.centerMinusWall, md.NMin, md.NMax, md.MRelRms, md.kBTWeighted, get_field(insertDiag,'nInsertedParticles',NaN), get_field(extractDiag,'nExtractedParticles',NaN));
    end
end
sampleSteps=sampleSteps(1:isamp); sampleTimes=sampleTimes(1:isamp); UxProfiles=UxProfiles(:,1:isamp); UyProfiles=UyProfiles(:,1:isamp); NProfiles=NProfiles(:,1:isamp); MProfiles=MProfiles(:,1:isamp); rows=rows(1:irow); timeseries=struct2table(rows);
caseDir=fullfile(opts.outputRoot,c.label); if ~exist(caseDir,'dir'), mkdir(caseDir); end
writetable(timeseries, fullfile(caseDir,'poiseuille_timeseries.csv'));
outCase=struct(); outCase.params=params; outCase.caseDef=c; outCase.initialInfo=initInfo; outCase.initialPoolInfo=poolInfo0; outCase.state=state; outCase.summary=timeseries; outCase.yCenters=prof.yCenters; outCase.sampleSteps=sampleSteps; outCase.sampleTimes=sampleTimes; outCase.UxProfiles=UxProfiles; outCase.UyProfiles=UyProfiles; outCase.NProfiles=NProfiles; outCase.MProfiles=MProfiles; outCase.finalMassDiagnostics=md; outCase.lastStepDiag=lastStepDiag; outCase.lastInsertDiag=lastInsert; outCase.lastExtractDiag=lastExtract; outCase.lastRemapDiag=lastRemap; outCase.lastThermostatDiag=lastThermostat;
try
    outCase.viscosity = analyze_resamp_poiseuille_viscosity(outCase,'excludeWallCells',opts.excludeWallCells,'fitStartFraction',opts.fitStartFraction);
catch ME
    outCase.viscosity = struct('nuEff',NaN,'R2',NaN,'signalToNoise',NaN,'errorMessage',ME.message);
end
save(fullfile(caseDir,'poiseuille_case.mat'),'outCase','-v7.3');
end

function params=make_params(opts,c)
params=struct(); params.Lx=opts.Nx; params.Ly=opts.Ny; params.Nx=opts.Nx; params.Ny=opts.Ny; params.gamma=opts.gamma; params.dt=opts.dt; params.alphaDeg=opts.alphaDeg; params.kBT=opts.kBT; params.resampParticleMass=opts.particleMass; params.resampTargetCellMass=opts.NTarget*opts.particleMass; params.resampNTarget=opts.NTarget; params.resampNMin=c.NMin; params.resampNMax=c.NMax; params.bodyForceX=opts.bodyForceX; params.bodyForceY=0; params.wallModeY=opts.wallModeY; params.Ubottom=0; params.Utop=0; params.useRandomGridShiftX=true; params.useRandomGridShiftY=false; params.useRandomGridShift=false; params.projectionEnable=~strcmp(c.method,'weighted_classic'); params.projectionStrength=c.projectionStrength; params.projectionInterpolationMethod=opts.projectionInterpolationMethod; params.projectionMomentumCorrectionEnable=true; params.projectionMomentumCorrectionMode='particle_global_exact'; params.thermostatAfterStep=c.thermostatAfterStep; params.thermostatAfterProjection=false; params.thermostatTargetKBT=opts.kBT; params.thermostatStrength=opts.thermostatStrength; params.thermostatMinParticlesPerCell=2; params.thermostatMaxScale=10; params.initialPopulationMode='exact_per_cell'; params.initialVelocityZeroGlobalMean=true; params.initialPoiseuilleProfileEnable=opts.initialPoiseuilleProfileEnable; params.initialPoiseuilleNuGuess=opts.initialPoiseuilleNuGuess; params.initialPoiseuilleScale=opts.initialPoiseuilleScale;
end

function [prof, md]=sample_state(state,params)
prof=resamp_poiseuille_profile_diagnostics(state,params);
md=resamp_population_mass_diagnostics(state,params,'periodicX',true,'periodicY',false,'cellWetMask',state.cellWetMask);
end

function row=make_timeseries_row(step,t,md,prof,stepDiag,insertDiag,extractDiag,remapDiag,thermostatDiag)
row=empty_timeseries_row(); row.step=step; row.t=t; row.NpActive=md.NpActive; row.Nfree=md.Nfree; row.NMin=md.NMin; row.NMax=md.NMax; row.NStd=md.NStd; row.MRelRms=md.MRelRms; row.mParticleRelStd=md.mParticleRelStd; row.mParticleMin=md.mParticleMin; row.mParticleMax=md.mParticleMax; row.kBTWeighted=md.kBTWeighted; row.meanVxWeighted=md.meanVxWeighted; row.centerMinusWall=prof.centerMinusWall; row.centerVelocity=prof.centerVelocity; row.wallMeanVelocity=prof.wallMeanVelocity; row.insertedParticles=get_field(insertDiag,'nInsertedParticles',NaN); row.insertedParticlesCumulative=get_field(insertDiag,'nInsertedParticlesCumulative',NaN); row.extractedParticles=get_field(extractDiag,'nExtractedParticles',NaN); row.extractedParticlesCumulative=get_field(extractDiag,'nExtractedParticlesCumulative',NaN); row.remapMassResidualRelRms=get_field(remapDiag,'massResidualRelRms',NaN); row.remapMomentumResidualRms=get_field(remapDiag,'momentumResidualRms',NaN); row.remapMassSafetyCells=get_field(remapDiag,'nCellsMassSafetyApplied',NaN); row.thermostatKBTAfter=get_field(thermostatDiag,'meanKBTAfter',NaN); row.rmsDivBefore=get_field(stepDiag,'rmsDivBefore',NaN); row.rmsDivAfter=get_field(stepDiag,'rmsDivParticleAfter',NaN);
end
function row=empty_timeseries_row()
row=struct('step',NaN,'t',NaN,'NpActive',NaN,'Nfree',NaN,'NMin',NaN,'NMax',NaN,'NStd',NaN,'MRelRms',NaN,'mParticleRelStd',NaN,'mParticleMin',NaN,'mParticleMax',NaN,'kBTWeighted',NaN,'meanVxWeighted',NaN,'centerMinusWall',NaN,'centerVelocity',NaN,'wallMeanVelocity',NaN,'insertedParticles',NaN,'insertedParticlesCumulative',NaN,'extractedParticles',NaN,'extractedParticlesCumulative',NaN,'remapMassResidualRelRms',NaN,'remapMomentumResidualRms',NaN,'remapMassSafetyCells',NaN,'thermostatKBTAfter',NaN,'rmsDivBefore',NaN,'rmsDivAfter',NaN);
end
function row=summarize_case(c,outCase,opts)
T=outCase.summary; last=T(end,:); visc=outCase.viscosity; row=empty_summary_row(); row.label={c.label}; row.method={c.method}; row.finalStep=last.step; row.finalTime=last.t; row.finalNMin=last.NMin; row.finalNMax=last.NMax; row.finalMRelRms=last.MRelRms; row.finalKBT=last.kBTWeighted; row.finalMassRelStd=last.mParticleRelStd; row.insertedCum=last.insertedParticlesCumulative; row.extractedCum=last.extractedParticlesCumulative; row.centerMinusWallFinal=last.centerMinusWall; row.nuEff=get_field(visc,'nuEff',NaN); row.R2=get_field(visc,'R2',NaN); row.signalToNoise=get_field(visc,'signalToNoise',NaN); row.physicalCandidate=double(get_field(visc,'physicalCandidate',false)); row.fitTMin=get_field(visc,'tMin',NaN); row.fitTMax=get_field(visc,'tMax',NaN); row.bodyForceX=opts.bodyForceX;
end
function row=empty_summary_row()
row=struct('label',{{''}},'method',{{''}},'finalStep',NaN,'finalTime',NaN,'finalNMin',NaN,'finalNMax',NaN,'finalMRelRms',NaN,'finalKBT',NaN,'finalMassRelStd',NaN,'insertedCum',NaN,'extractedCum',NaN,'centerMinusWallFinal',NaN,'nuEff',NaN,'R2',NaN,'signalToNoise',NaN,'physicalCandidate',NaN,'fitTMin',NaN,'fitTMax',NaN,'bodyForceX',NaN);
end

function s=empty_insert_stats(), s=struct('cumulativeParticles',0,'cumulativeCells',0,'maxParticlesPerStep',0,'lastStep',NaN); end
function s=update_insert_stats(s,d,step), n=get_field(d,'nInsertedParticles',0); c=get_field(d,'nCellsInserted',0); if isfinite(n)&&n>0, s.cumulativeParticles=s.cumulativeParticles+n; s.cumulativeCells=s.cumulativeCells+c; s.maxParticlesPerStep=max(s.maxParticlesPerStep,n); s.lastStep=step; end; end
function d=attach_insert_cumulative(d,s), d.nInsertedParticlesCumulative=s.cumulativeParticles; d.nInsertedCellsCumulative=s.cumulativeCells; d.maxInsertedParticlesPerStep=s.maxParticlesPerStep; d.lastInsertionStep=s.lastStep; end
function d=empty_insert_diag(), d=struct('nInsertedParticles',0,'nCellsInserted',0,'nPoorCellsAfter',NaN,'nOverCellsAfter',NaN); end
function s=empty_extract_stats(), s=struct('cumulativeParticles',0,'cumulativeCells',0,'maxParticlesPerStep',0,'lastStep',NaN); end
function s=update_extract_stats(s,d,step), n=get_field(d,'nExtractedParticles',0); c=get_field(d,'nCellsExtracted',0); if isfinite(n)&&n>0, s.cumulativeParticles=s.cumulativeParticles+n; s.cumulativeCells=s.cumulativeCells+c; s.maxParticlesPerStep=max(s.maxParticlesPerStep,n); s.lastStep=step; end; end
function d=attach_extract_cumulative(d,s), d.nExtractedParticlesCumulative=s.cumulativeParticles; d.nExtractedCellsCumulative=s.cumulativeCells; d.maxExtractedParticlesPerStep=s.maxParticlesPerStep; d.lastExtractionStep=s.lastStep; end
function d=empty_extract_diag(), d=struct('nExtractedParticles',0,'nCellsExtracted',0,'nPoorCellsAfter',NaN,'nOverCellsAfter',NaN); end
function d=empty_remap_diag(), d=struct('massResidualRelRms',NaN,'momentumResidualRms',NaN,'nCellsMassSafetyApplied',NaN); end
function v=get_field(s,name,defaultValue)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=defaultValue; end
end
function opts=parse_options(varargin)
opts=struct(); opts.outputRoot=fullfile('runs','resamp_poiseuille_physics_validation'); opts.Nx=64; opts.Ny=32; opts.gamma=20; opts.NTarget=20; opts.NMin=14; opts.NMax=26; opts.steps=5000; opts.sampleEvery=50; opts.summaryEvery=100; opts.dt=0.005; opts.alphaDeg=90; opts.kBT=0.01; opts.particleMass=1.0; opts.bodyForceX=0.005; opts.wallModeY='bounceback'; opts.capacityFactor=2.0; opts.thermostatStrength=0.25; opts.projectionInterpolationMethod='nearest'; opts.extractSelectionMode='closest_to_cell_mean'; opts.insertVelocityMode='current_or_memory_pairwise'; opts.memoryMinParticles=14; opts.remapMethod='scale_preserve_velocity'; opts.massMinFactor=0.05; opts.massMaxFactor=20.0; opts.remapMassSafetyMode='uniform_mass_velocity_shift'; opts.massSafetyMinFactor=0.25; opts.massSafetyMaxFactor=4.0; opts.constraintTolerance=1e-10; opts.preservePreEditVelocity=true; opts.initialPoiseuilleProfileEnable=false; opts.initialPoiseuilleNuGuess=0.05; opts.initialPoiseuilleScale=1.0; opts.excludeWallCells=2; opts.fitStartFraction=0.5; opts.rngSeed=12345;
if mod(numel(varargin),2)~=0, error('Options must be name/value pairs.'); end
for k=1:2:numel(varargin)
    key=lower(char(string(varargin{k}))); val=varargin{k+1};
    switch key
        case 'outputroot', opts.outputRoot=char(string(val));
        case 'nx', opts.Nx=val; case 'ny', opts.Ny=val; case 'gamma', opts.gamma=val; opts.NTarget=val; case 'ntarget', opts.NTarget=val; case 'nmin', opts.NMin=val; case 'nmax', opts.NMax=val; case 'steps', opts.steps=val; case 'sampleevery', opts.sampleEvery=val; case 'summaryevery', opts.summaryEvery=val; case 'dt', opts.dt=val; case 'alphadeg', opts.alphaDeg=val; case 'kbt', opts.kBT=val; case 'particlemass', opts.particleMass=val; case 'bodyforcex', opts.bodyForceX=val; case 'wallmodey', opts.wallModeY=char(string(val)); case 'capacityfactor', opts.capacityFactor=val; case 'thermostatstrength', opts.thermostatStrength=val; case 'projectioninterpolationmethod', opts.projectionInterpolationMethod=char(string(val)); case 'extractselectionmode', opts.extractSelectionMode=lower(char(string(val))); case 'insertvelocitymode', opts.insertVelocityMode=lower(char(string(val))); case 'memoryminparticles', opts.memoryMinParticles=val; case 'remapmethod', opts.remapMethod=lower(char(string(val))); case 'massminfactor', opts.massMinFactor=val; case 'massmaxfactor', opts.massMaxFactor=val; case 'remapmasssafetymode', opts.remapMassSafetyMode=lower(char(string(val))); case 'masssafetyminfactor', opts.massSafetyMinFactor=val; case 'masssafetymaxfactor', opts.massSafetyMaxFactor=val; case 'constrainttolerance', opts.constraintTolerance=val; case 'preservepreeditvelocity', opts.preservePreEditVelocity=logical(val); case 'initialpoiseuilleprofileenable', opts.initialPoiseuilleProfileEnable=logical(val); case 'initialpoiseuillenuguess', opts.initialPoiseuilleNuGuess=val; case 'initialpoiseuillescale', opts.initialPoiseuilleScale=val; case 'excludewallcells', opts.excludeWallCells=val; case 'fitstartfraction', opts.fitStartFraction=val; case 'rngseed', opts.rngSeed=val;
        otherwise, error('Unknown option: %s', key);
    end
end
end
