function out = run_resamp_poiseuille_wallvp_validation(varargin)
%RUN_RESAMP_POISEUILLE_WALLVP_VALIDATION Weighted-resampling Poiseuille with wall VPs.
%
% This runner is intentionally narrow: it validates Poiseuille only after
% adding weighted virtual wall particles to the SRC collision.  It compares a
% classic weighted wallVP reference and the Q6/resampled wallVP method.

opts = parse_options(varargin{:});
if ~exist(opts.outputRoot,'dir'), mkdir(opts.outputRoot); end
caseDefs = build_cases(opts);
outs = struct(); rows = repmat(empty_summary_row(), numel(caseDefs), 1);
for ic = 1:numel(caseDefs)
    c = caseDefs(ic);
    fprintf('\n=== Weighted Poiseuille wallVP case: %s ===\n', c.label);
    outCase = run_one_case(c, opts, ic);
    outs.(c.label) = outCase;
    rows(ic) = summarize_case(c, outCase, opts);
end
summary = struct2table(rows);
writetable(summary, fullfile(opts.outputRoot,'resamp_poiseuille_wallvp_summary.csv'));
out = struct('options',opts,'caseDefs',caseDefs,'cases',outs,'summary',summary,'outputRoot',opts.outputRoot);
save(fullfile(opts.outputRoot,'resamp_poiseuille_wallvp_validation.mat'),'out','-v7.3');
fprintf('\nWeighted Poiseuille wallVP summary:\n'); disp(summary);
fprintf('Wrote %s\n', fullfile(opts.outputRoot,'resamp_poiseuille_wallvp_summary.csv'));
end

function caseDefs = build_cases(opts)
allCases = struct([]);
allCases(1).label = 'classic_wallvp_reference';
allCases(1).method = 'weighted_classic';
allCases(1).projectionStrength = 0;
allCases(1).extractEvery = 0;
allCases(1).insertEvery = 0;
allCases(1).remapEvery = 0;
allCases(1).thermostatAfterStep = true;
allCases(1).thermostatAfterRemap = false;
allCases(1).massSafetyEnable = false;
allCases(1).massRegularizationEnable = false;
allCases(1).NMin = 0;
allCases(1).NMax = 2 * opts.NTarget;

allCases(2).label = 'q6_resampled_wallvp';
allCases(2).method = 'weighted_q6';
allCases(2).projectionStrength = 1;
allCases(2).extractEvery = 1;
allCases(2).insertEvery = 1;
allCases(2).remapEvery = 1;
allCases(2).thermostatAfterStep = false;
allCases(2).thermostatAfterRemap = true;
allCases(2).massSafetyEnable = true;
allCases(2).massRegularizationEnable = opts.massRegularizationEnable;
allCases(2).NMin = opts.NMin;
allCases(2).NMax = opts.NMax;

wanted = normalize_case_list(opts.cases);
caseDefs = allCases([]);
n = 0;
for i = 1:numel(allCases)
    if any(strcmp(allCases(i).label, wanted))
        n = n + 1;
        caseDefs(n) = allCases(i); %#ok<AGROW>
    end
end
if isempty(caseDefs)
    error('No Poiseuille wallVP cases selected. Valid cases are: classic_wallvp_reference, q6_resampled_wallvp.');
end
end

function wanted = normalize_case_list(cases)
if ischar(cases) || isstring(cases)
    cases = cellstr(string(cases));
end
if isempty(cases)
    cases = {'classic_wallvp_reference','q6_resampled_wallvp'};
end
wanted = cell(size(cases));
for i = 1:numel(cases)
    key = lower(strtrim(char(string(cases{i}))));
    switch key
        case {'classic_wallvp_reference','src_classic_wallvp_only','classic','weighted_classic'}
            wanted{i} = 'classic_wallvp_reference';
        case {'q6_resampled_wallvp','q6','weighted_q6','resampled'}
            wanted{i} = 'q6_resampled_wallvp';
        case {'all','both'}
            wanted = {'classic_wallvp_reference','q6_resampled_wallvp'};
            return;
        otherwise
            error('Unknown Poiseuille wallVP case: %s', key);
    end
end
wanted = unique(wanted, 'stable');
end

function outCase = run_one_case(c, opts, ic)
if ~isempty(opts.rngSeed), rng(opts.rngSeed + ic - 1, 'twister'); end
params = make_params(opts, c);
[state, initInfo] = resamp_initialize_particles_poiseuille_weighted(params);
[state, poolInfo0] = resamp_enable_particle_pool(state, 'capacityFactor', opts.capacityFactor);
state.uMemUx = zeros(params.Nx, params.Ny); state.uMemUy = zeros(params.Nx, params.Ny); state.uMemValid = false(params.Nx, params.Ny);
state.cellWetMask = true(params.Nx, params.Ny); params.cellWetMask = state.cellWetMask;

nSamples = floor(opts.steps / opts.sampleEvery) + 1;
UxProfiles = nan(params.Ny, nSamples); UyProfiles = nan(params.Ny, nSamples); NProfiles = nan(params.Ny, nSamples); MProfiles = nan(params.Ny, nSamples);
sampleTimes = nan(nSamples,1); sampleSteps = nan(nSamples,1);
rows = repmat(empty_timeseries_row(), floor(opts.steps / opts.summaryEvery) + 2, 1);
irow = 0; isamp = 0;
insertStats = empty_counter_stats(); extractStats = empty_counter_stats();
lastStepDiag = struct(); lastInsert = empty_insert_diag(); lastExtract = empty_extract_diag(); lastRemap = empty_remap_diag(); lastThermostat = empty_thermostat_diag(); lastMassReg = empty_mass_regularization_diag();
massRegStats = empty_counter_stats();

[prof, md] = sample_state(state, params);
isamp = isamp + 1; sampleSteps(isamp)=0; sampleTimes(isamp)=0; UxProfiles(:,isamp)=prof.Ux; UyProfiles(:,isamp)=prof.Uy; NProfiles(:,isamp)=prof.N; MProfiles(:,isamp)=prof.M;
irow = irow + 1; rows(irow) = make_timeseries_row(0, 0, md, prof, struct(), lastInsert, lastExtract, lastRemap, lastThermostat, lastMassReg);
if opts.visualEvery > 0
    visualize_case_frame(state, params, prof, md, c, opts, 0, 0, struct(), false);
end

for step = 1:opts.steps
    switch c.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_poiseuille_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag] = resamp_step_projection_poiseuille_weighted(state, params);
        otherwise
            error('Unknown method: %s', c.method);
    end
    lastStepDiag = stepDiag;

    preEditG = [];
    if opts.preservePreEditVelocity && c.remapEvery > 0 && mod(step, c.remapEvery) == 0
        preEditG = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
            'periodicX', true, 'periodicY', false, 'activeMask', resamp_active_mask(state), 'cellWetMask', state.cellWetMask);
    end

    if c.extractEvery > 0 && mod(step, c.extractEvery) == 0
        [state, lastExtract] = resamp_extract_overpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, 'NMin', c.NMin, 'NMax', c.NMax, ...
            'selectionMode', opts.extractSelectionMode, 'periodicX', true, 'periodicY', false, 'cellWetMask', state.cellWetMask);
        extractStats = update_counter_stats(extractStats, get_field(lastExtract,'nExtractedParticles',0), get_field(lastExtract,'nCellsExtracted',0), step);
        lastExtract = attach_extract_cumulative(lastExtract, extractStats);
    else
        lastExtract = attach_extract_cumulative(empty_extract_diag(), extractStats);
    end

    if c.insertEvery > 0 && mod(step, c.insertEvery) == 0
        [state, lastInsert] = resamp_insert_underpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, 'NMin', c.NMin, 'NMax', c.NMax, ...
            'memoryMinParticles', opts.memoryMinParticles, 'insertVelocityMode', opts.insertVelocityMode, ...
            'periodicX', true, 'periodicY', false, 'cellWetMask', state.cellWetMask);
        insertStats = update_counter_stats(insertStats, get_field(lastInsert,'nInsertedParticles',0), get_field(lastInsert,'nCellsInserted',0), step);
        lastInsert = attach_insert_cumulative(lastInsert, insertStats);
    else
        lastInsert = attach_insert_cumulative(empty_insert_diag(), insertStats);
    end

    if c.remapEvery > 0 && mod(step, c.remapEvery) == 0
        if ~isempty(preEditG) && (get_field(lastInsert,'nInsertedParticles',0) > 0 || get_field(lastExtract,'nExtractedParticles',0) > 0 || opts.alwaysUsePreEditVelocityForRemap)
            [state, lastRemap] = resamp_local_mass_moment_remap(state, params, ...
                'periodicX', true, 'periodicY', false, 'cellWetMask', state.cellWetMask, ...
                'targetVelocityMode', 'grid', 'targetUx', preEditG.Ux, 'targetUy', preEditG.Uy, ...
                'method', opts.remapMethod, 'massMin', opts.massMinFactor * opts.particleMass, 'massMax', opts.massMaxFactor * opts.particleMass, ...
                'massSafetyEnable', c.massSafetyEnable, 'massSafetyMode', opts.remapMassSafetyMode, ...
                'massSafetyMinFactor', opts.massSafetyMinFactor, 'massSafetyMaxFactor', opts.massSafetyMaxFactor, ...
                'constraintTolerance', opts.constraintTolerance);
        else
            [state, lastRemap] = resamp_local_mass_moment_remap(state, params, ...
                'periodicX', true, 'periodicY', false, 'cellWetMask', state.cellWetMask, ...
                'method', opts.remapMethod, 'massMin', opts.massMinFactor * opts.particleMass, 'massMax', opts.massMaxFactor * opts.particleMass, ...
                'massSafetyEnable', c.massSafetyEnable, 'massSafetyMode', opts.remapMassSafetyMode, ...
                'massSafetyMinFactor', opts.massSafetyMinFactor, 'massSafetyMaxFactor', opts.massSafetyMaxFactor, ...
                'constraintTolerance', opts.constraintTolerance);
        end
    else
        lastRemap = empty_remap_diag();
    end

    if c.thermostatAfterRemap
        [state.v, lastThermostat] = resamp_apply_cell_thermostat_weighted(state.x, state.v, state.m, params, ...
            'periodicX', true, 'periodicY', false, 'activeMask', resamp_active_mask(state), 'cellWetMask', state.cellWetMask);
    else
        lastThermostat = empty_thermostat_diag();
    end

    if c.massRegularizationEnable && should_apply_mass_regularization(step, opts)
        [state, lastMassReg] = resamp_regularize_particle_masses_preserve_mue(state, params, ...
            'periodicX', true, 'periodicY', false, 'activeMask', resamp_active_mask(state), 'cellWetMask', state.cellWetMask, ...
            'strength', opts.massRegularizationStrength, ...
            'triggerMode', opts.massRegularizationTriggerMode, ...
            'relStdTrigger', opts.massRegularizationRelStdTrigger, ...
            'minFactor', opts.massRegularizationMinFactor, ...
            'maxFactor', opts.massRegularizationMaxFactor, ...
            'preserveEnergy', opts.massRegularizationPreserveEnergy, ...
            'nominalParticleMass', opts.particleMass, ...
            'minParticlesPerCell', opts.massRegularizationMinParticlesPerCell, ...
            'tolerance', opts.constraintTolerance);
        lastMassReg.step = step;
        massRegStats = update_counter_stats(massRegStats, get_field(lastMassReg,'nParticlesRegularized',0), get_field(lastMassReg,'nCellsRegularized',0), step);
        lastMassReg = attach_mass_regularization_cumulative(lastMassReg, massRegStats);
    else
        lastMassReg = attach_mass_regularization_cumulative(empty_mass_regularization_diag(), massRegStats);
    end

    if mod(step, opts.sampleEvery) == 0
        [prof, md] = sample_state(state, params);
        isamp = isamp + 1; sampleSteps(isamp)=step; sampleTimes(isamp)=step*opts.dt; UxProfiles(:,isamp)=prof.Ux; UyProfiles(:,isamp)=prof.Uy; NProfiles(:,isamp)=prof.N; MProfiles(:,isamp)=prof.M;
    end
    if mod(step, opts.summaryEvery) == 0 || step == opts.steps
        [prof, md] = sample_state(state, params);
        irow = irow + 1; rows(irow) = make_timeseries_row(step, step*opts.dt, md, prof, lastStepDiag, lastInsert, lastExtract, lastRemap, lastThermostat, lastMassReg);
        fprintf('[%s] step=%d/%d N=[%g,%g] Mrel=%.3g kBT=%.4g wallVP=%.0f center-wall=%.4g\n', ...
            c.label, step, opts.steps, md.NMin, md.NMax, md.MRelRms, md.kBTWeighted, get_field(stepDiag,'nVirtualWallParticlesEquivalent',NaN), prof.centerMinusWall);
    end
    if opts.visualEvery > 0 && (mod(step, opts.visualEvery) == 0 || step == opts.steps)
        if ~(exist('prof','var') && isstruct(prof)) || get_field(prof,'lastVisualStep',NaN) ~= step
            [profVis, mdVis] = sample_state(state, params);
        else
            profVis = prof; mdVis = md;
        end
        visualize_case_frame(state, params, profVis, mdVis, c, opts, step, step*opts.dt, stepDiag, false);
    end
end

rows = rows(1:irow); summary = struct2table(rows);
UxProfiles = UxProfiles(:,1:isamp); UyProfiles = UyProfiles(:,1:isamp); NProfiles = NProfiles(:,1:isamp); MProfiles = MProfiles(:,1:isamp); sampleTimes = sampleTimes(1:isamp); sampleSteps = sampleSteps(1:isamp);

caseDir = fullfile(opts.outputRoot, c.label); if ~exist(caseDir,'dir'), mkdir(caseDir); end
writetable(summary, fullfile(caseDir,'timeseries.csv'));
profileTable = table((1:params.Ny)', mean(UxProfiles(:,max(1,ceil(0.5*isamp)):isamp),2,'omitnan'), mean(UyProfiles(:,max(1,ceil(0.5*isamp)):isamp),2,'omitnan'), mean(NProfiles(:,max(1,ceil(0.5*isamp)):isamp),2,'omitnan'), mean(MProfiles(:,max(1,ceil(0.5*isamp)):isamp),2,'omitnan'), ...
    'VariableNames', {'iy','UxMeanLate','UyMeanLate','NMeanLate','MMeanLate'});
writetable(profileTable, fullfile(caseDir,'profile_late_mean.csv'));

try
    fitOut = struct();
    fitOut.params = params;
    fitOut.yCenters = ((1:params.Ny)' - 0.5) .* (params.Ly / params.Ny);
    fitOut.UxProfiles = UxProfiles;
    fitOut.NProfiles = NProfiles;
    fitOut.sampleTimes = sampleTimes;
    visc = analyze_resamp_poiseuille_viscosity(fitOut, ...
        'excludeWallCells', opts.excludeWallCells, ...
        'fitStartFraction', opts.fitStartFraction, ...
        'fitModel', opts.poiseuilleFitModel);
catch ME
    warning('Viscosity analysis failed for %s: %s', c.label, ME.message);
    visc = struct('nuEff', NaN, 'R2', NaN, 'physicalCandidate', false, 'signalToNoise', NaN, 'errorMessage', ME.message);
end

profileTable = add_fit_columns_to_profile_table(profileTable, visc);
writetable(profileTable, fullfile(caseDir,'profile_late_mean.csv'));
if opts.saveFinalFigures || opts.visualEvery > 0
    [profFinal, mdFinal] = sample_state(state, params);
    visualize_case_frame(state, params, profFinal, mdFinal, c, opts, opts.steps, opts.steps*opts.dt, lastStepDiag, true, visc, caseDir);
end

outCase = struct('caseDef', c, 'params', params, 'initInfo', initInfo, 'poolInfo0', poolInfo0, 'stateFinal', state, 'summary', summary, ...
    'sampleSteps', sampleSteps, 'sampleTimes', sampleTimes, 'UxProfiles', UxProfiles, 'UyProfiles', UyProfiles, 'NProfiles', NProfiles, 'MProfiles', MProfiles, 'profileTable', profileTable, 'viscosity', visc);
save(fullfile(caseDir,'poiseuille_wallvp_case.mat'),'outCase','-v7.3');
end

function params = make_params(opts, c)
params = struct();
params.Lx=opts.Nx; params.Ly=opts.Ny; params.Nx=opts.Nx; params.Ny=opts.Ny; params.gamma=opts.gamma;
params.dt=opts.dt; params.alphaDeg=opts.alphaDeg; params.kBT=opts.kBT; params.resampParticleMass=opts.particleMass;
params.resampTargetCellMass=opts.NTarget*opts.particleMass; params.resampNTarget=opts.NTarget; params.resampNMin=c.NMin; params.resampNMax=c.NMax;
params.bodyForceX=opts.bodyForceX; params.bodyForceY=0; params.wallModeY=opts.wallModeY; params.Ubottom=0; params.Utop=0;
params.useRandomGridShiftX=true; params.useRandomGridShiftY=false; params.useRandomGridShift=false;
params.projectionEnable=~strcmp(c.method,'weighted_classic'); params.projectionStrength=c.projectionStrength; params.projectionInterpolationMethod=opts.projectionInterpolationMethod;
params.projectionMomentumCorrectionEnable=true; params.projectionMomentumCorrectionMode='particle_global_exact';
params.thermostatAfterStep=c.thermostatAfterStep; params.thermostatAfterProjection=false; params.thermostatTargetKBT=opts.kBT; params.thermostatStrength=opts.thermostatStrength; params.thermostatMinParticlesPerCell=2; params.thermostatMaxScale=10;
params.initialPopulationMode='exact_per_cell'; params.initialVelocityZeroGlobalMean=true; params.initialPoiseuilleProfileEnable=opts.initialPoiseuilleProfileEnable; params.initialPoiseuilleNuGuess=opts.initialPoiseuilleNuGuess; params.initialPoiseuilleScale=opts.initialPoiseuilleScale;
params.wallVirtualParticlesEnable=opts.wallVirtualParticlesEnable; params.wallVirtualParticlesGeometryMode=opts.wallVirtualParticlesGeometryMode; params.wallVirtualParticlesForceRandomShiftY=opts.wallVirtualParticlesForceRandomShiftY;
params.wallVirtualParticlesDensityFactor=opts.wallVirtualParticlesDensityFactor; params.wallVirtualParticlesPerCell=opts.wallVirtualParticlesPerCell; params.wallVirtualParticleMass=opts.wallVirtualParticleMass;
params.wallVirtualParticlesThermal=opts.wallVirtualParticlesThermal; params.wallVirtualParticlesKBT=opts.wallVirtualParticlesKBT; params.wallVirtualParticlesStochasticCount=opts.wallVirtualParticlesStochasticCount;
params.wallVirtualParticlesIncludeBottom=true; params.wallVirtualParticlesIncludeTop=true;
end

function [prof, md] = sample_state(state, params)
prof = resamp_poiseuille_profile_diagnostics(state, params);
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', false, 'cellWetMask', state.cellWetMask);
end

function row = make_timeseries_row(step,t,md,prof,stepDiag,insertDiag,extractDiag,remapDiag,thermostatDiag,massRegDiag)
row = empty_timeseries_row(); row.step=step; row.t=t; row.NpActive=md.NpActive; row.Nfree=md.Nfree; row.NMin=md.NMin; row.NMax=md.NMax; row.NStd=md.NStd; row.MRelRms=md.MRelRms; row.mParticleRelStd=md.mParticleRelStd; row.mParticleMin=md.mParticleMin; row.mParticleMax=md.mParticleMax; row.kBTWeighted=md.kBTWeighted; row.meanVxWeighted=md.meanVxWeighted; row.centerMinusWall=prof.centerMinusWall; row.centerVelocity=prof.centerVelocity; row.wallMeanVelocity=prof.wallMeanVelocity; row.insertedParticles=get_field(insertDiag,'nInsertedParticles',NaN); row.insertedParticlesCumulative=get_field(insertDiag,'nInsertedParticlesCumulative',NaN); row.extractedParticles=get_field(extractDiag,'nExtractedParticles',NaN); row.extractedParticlesCumulative=get_field(extractDiag,'nExtractedParticlesCumulative',NaN); row.remapMassResidualRelRms=get_field(remapDiag,'massResidualRelRms',NaN); row.remapMomentumResidualRms=get_field(remapDiag,'momentumResidualRms',NaN); row.remapMassSafetyCells=get_field(remapDiag,'nCellsMassSafetyApplied',NaN); row.thermostatKBTAfter=get_field(thermostatDiag,'meanKBTAfter',NaN); row.rmsDivBefore=get_field(stepDiag,'rmsDivBefore',NaN); row.rmsDivAfter=get_field(stepDiag,'rmsDivParticleAfter',NaN); row.nVirtualWallCells=get_field(stepDiag,'nVirtualWallCells',get_field(get_field(stepDiag,'classic',struct()),'nVirtualWallCells',NaN)); row.nVirtualWallParticlesEquivalent=get_field(stepDiag,'nVirtualWallParticlesEquivalent',get_field(get_field(stepDiag,'classic',struct()),'nVirtualWallParticlesEquivalent',NaN)); row.virtualWallMassTotal=get_field(stepDiag,'virtualWallMassTotal',get_field(get_field(stepDiag,'classic',struct()),'virtualWallMassTotal',NaN)); row.massRegCells=get_field(massRegDiag,'nCellsRegularized',NaN); row.massRegParticles=get_field(massRegDiag,'nParticlesRegularized',NaN); row.massRegCellsCumulative=get_field(massRegDiag,'nCellsRegularizedCumulative',NaN); row.massRegParticlesCumulative=get_field(massRegDiag,'nParticlesRegularizedCumulative',NaN); row.massRegRelStdBefore=get_field(massRegDiag,'meanRelStdBefore',NaN); row.massRegRelStdAfter=get_field(massRegDiag,'meanRelStdAfter',NaN); row.massRegEnergyResidualRelRms=get_field(massRegDiag,'energyResidualRelRms',NaN); row.massRegMomentumResidualRms=get_field(massRegDiag,'momentumResidualRms',NaN);
end
function row = empty_timeseries_row()
row = struct('step',NaN,'t',NaN,'NpActive',NaN,'Nfree',NaN,'NMin',NaN,'NMax',NaN,'NStd',NaN,'MRelRms',NaN,'mParticleRelStd',NaN,'mParticleMin',NaN,'mParticleMax',NaN,'kBTWeighted',NaN,'meanVxWeighted',NaN,'centerMinusWall',NaN,'centerVelocity',NaN,'wallMeanVelocity',NaN,'insertedParticles',NaN,'insertedParticlesCumulative',NaN,'extractedParticles',NaN,'extractedParticlesCumulative',NaN,'remapMassResidualRelRms',NaN,'remapMomentumResidualRms',NaN,'remapMassSafetyCells',NaN,'thermostatKBTAfter',NaN,'rmsDivBefore',NaN,'rmsDivAfter',NaN,'nVirtualWallCells',NaN,'nVirtualWallParticlesEquivalent',NaN,'virtualWallMassTotal',NaN,'massRegCells',NaN,'massRegParticles',NaN,'massRegCellsCumulative',NaN,'massRegParticlesCumulative',NaN,'massRegRelStdBefore',NaN,'massRegRelStdAfter',NaN,'massRegEnergyResidualRelRms',NaN,'massRegMomentumResidualRms',NaN);
end
function row = summarize_case(c,outCase,opts)
T = outCase.summary;
last = T(end,:);
visc = outCase.viscosity;
row = empty_summary_row();
row.label={c.label}; row.method={c.method}; row.finalStep=last.step; row.finalTime=last.t; row.finalNMin=last.NMin; row.finalNMax=last.NMax; row.finalMRelRms=last.MRelRms; row.finalKBT=last.kBTWeighted; row.finalMassRelStd=last.mParticleRelStd; row.insertedCum=last.insertedParticlesCumulative; row.extractedCum=last.extractedParticlesCumulative; row.centerMinusWallFinal=last.centerMinusWall; row.wallMeanVelocityFinal=last.wallMeanVelocity; row.nuEff=get_field(visc,'nuEff',NaN); row.R2=get_field(visc,'R2',NaN); row.signalToNoise=get_field(visc,'signalToNoise',NaN); row.physicalCandidate=double(get_field(visc,'physicalCandidate',false)); row.fitTMin=get_field(visc,'tMin',NaN); row.fitTMax=get_field(visc,'tMax',NaN); row.bodyForceX=opts.bodyForceX; row.wallVirtualParticlesEnable=double(opts.wallVirtualParticlesEnable); row.nVirtualWallParticlesEquivalentFinal=last.nVirtualWallParticlesEquivalent; row.virtualWallMassTotalFinal=last.virtualWallMassTotal; row.poiseuilleFitModel={get_field(visc,'fitModel','')}; row.nuEffSlip=get_field(visc,'nuEffSlip',NaN); row.R2Slip=get_field(visc,'R2Slip',NaN); row.nuEffNoSlip=get_field(visc,'nuEffNoSlip',NaN); row.R2NoSlip=get_field(visc,'R2NoSlip',NaN); row.slipVelocity=get_field(visc,'slipVelocity',NaN); row.slipRatio=get_field(visc,'slipRatio',NaN); row.wallSlipRatioObserved=get_field(visc,'wallSlipRatioObserved',NaN);
if any(strcmp('massRegCellsCumulative', T.Properties.VariableNames))
    row.massRegCellsCum = last.massRegCellsCumulative;
    row.massRegParticlesCum = last.massRegParticlesCumulative;
    idxReg = find(T.massRegCells > 0, 1, 'last');
    if ~isempty(idxReg)
        row.massRegRelStdBeforeLast = T.massRegRelStdBefore(idxReg);
        row.massRegRelStdAfterLast = T.massRegRelStdAfter(idxReg);
    end
end
end
function row = empty_summary_row()
row = struct('label',{{''}},'method',{{''}},'finalStep',NaN,'finalTime',NaN,'finalNMin',NaN,'finalNMax',NaN,'finalMRelRms',NaN,'finalKBT',NaN,'finalMassRelStd',NaN,'insertedCum',NaN,'extractedCum',NaN,'centerMinusWallFinal',NaN,'wallMeanVelocityFinal',NaN,'nuEff',NaN,'R2',NaN,'signalToNoise',NaN,'physicalCandidate',NaN,'fitTMin',NaN,'fitTMax',NaN,'bodyForceX',NaN,'wallVirtualParticlesEnable',NaN,'nVirtualWallParticlesEquivalentFinal',NaN,'virtualWallMassTotalFinal',NaN,'poiseuilleFitModel',{{''}},'nuEffSlip',NaN,'R2Slip',NaN,'nuEffNoSlip',NaN,'R2NoSlip',NaN,'slipVelocity',NaN,'slipRatio',NaN,'wallSlipRatioObserved',NaN,'massRegCellsCum',NaN,'massRegParticlesCum',NaN,'massRegRelStdBeforeLast',NaN,'massRegRelStdAfterLast',NaN);
end


function profileTable = add_fit_columns_to_profile_table(profileTable, visc)
if isstruct(visc) && isfield(visc,'y') && height(profileTable) == numel(visc.y)
    profileTable.fitMask = logical(visc.fitMask(:));
    profileTable.UxFitSlip = visc.UxFitSlip(:);
    profileTable.UxFitNoSlip = visc.UxFitNoSlip(:);
    profileTable.UxFitChosen = visc.UxFit(:);
else
    profileTable.fitMask = false(height(profileTable),1);
    profileTable.UxFitSlip = nan(height(profileTable),1);
    profileTable.UxFitNoSlip = nan(height(profileTable),1);
    profileTable.UxFitChosen = nan(height(profileTable),1);
end
end

function visualize_case_frame(state, params, prof, md, c, opts, step, t, stepDiag, isFinal, varargin)
fitInfo = [];
caseDir = '';
if nargin >= 11 && ~isempty(varargin)
    fitInfo = varargin{1};
end
if nargin >= 12 && numel(varargin) >= 2
    caseDir = varargin{2};
end
saveFrame = opts.saveFrames && step > 0 && mod(step, opts.saveFrameEvery) == 0;
if isFinal && opts.saveFinalFigures
    saveFrame = true;
end
if isempty(caseDir)
    frameDir = fullfile(opts.outputRoot, c.label, opts.frameDirName);
else
    frameDir = fullfile(caseDir, opts.frameDirName);
end
resamp_poiseuille_visualize_frame(state, params, prof, md, ...
    'caseLabel', c.label, ...
    'step', step, ...
    't', t, ...
    'figureId', opts.figureId + 10*(find_case_offset(c.label)), ...
    'profileFigureId', opts.profileFigureId + 10*(find_case_offset(c.label)), ...
    'showProfileFigure', opts.showProfileFigure, ...
    'fitInfo', fitInfo, ...
    'saveFrame', saveFrame, ...
    'frameDir', frameDir, ...
    'particleMarkerSize', opts.particleMarkerSize, ...
    'nVirtualWallParticlesEquivalent', get_field(stepDiag,'nVirtualWallParticlesEquivalent',get_field(get_field(stepDiag,'classic',struct()),'nVirtualWallParticlesEquivalent',NaN)), ...
    'pngResolution', opts.pngResolution);
end

function k = find_case_offset(label)
if strcmp(label,'classic_wallvp_reference')
    k = 0;
elseif strcmp(label,'q6_resampled_wallvp')
    k = 1;
else
    k = 2;
end
end

function s = empty_counter_stats(), s=struct('cumulativeParticles',0,'cumulativeCells',0,'maxParticlesPerStep',0,'lastStep',NaN); end
function s = update_counter_stats(s,n,c,step), if isfinite(n)&&n>0, s.cumulativeParticles=s.cumulativeParticles+n; s.cumulativeCells=s.cumulativeCells+c; s.maxParticlesPerStep=max(s.maxParticlesPerStep,n); s.lastStep=step; end; end
function d = attach_insert_cumulative(d,s), d.nInsertedParticlesCumulative=s.cumulativeParticles; d.nInsertedCellsCumulative=s.cumulativeCells; d.maxInsertedParticlesPerStep=s.maxParticlesPerStep; d.lastInsertionStep=s.lastStep; end
function d = attach_extract_cumulative(d,s), d.nExtractedParticlesCumulative=s.cumulativeParticles; d.nExtractedCellsCumulative=s.cumulativeCells; d.maxExtractedParticlesPerStep=s.maxParticlesPerStep; d.lastExtractionStep=s.lastStep; end
function d = attach_mass_regularization_cumulative(d,s), d.nParticlesRegularizedCumulative=s.cumulativeParticles; d.nCellsRegularizedCumulative=s.cumulativeCells; d.maxParticlesRegularizedPerStep=s.maxParticlesPerStep; d.lastMassRegularizationStep=s.lastStep; end
function d = empty_insert_diag(), d=struct('nInsertedParticles',0,'nCellsInserted',0,'nPoorCellsAfter',NaN,'nOverCellsAfter',NaN); end
function d = empty_extract_diag(), d=struct('nExtractedParticles',0,'nCellsExtracted',0,'nPoorCellsAfter',NaN,'nOverCellsAfter',NaN); end
function d = empty_remap_diag(), d=struct('massResidualRelRms',NaN,'momentumResidualRms',NaN,'nCellsMassSafetyApplied',NaN); end
function d = empty_thermostat_diag(), d=struct('enabled',false,'meanKBTAfter',NaN); end
function d = empty_mass_regularization_diag(), d=struct('enabled',false,'nCellsRegularized',0,'nParticlesRegularized',0,'meanRelStdBefore',NaN,'meanRelStdAfter',NaN,'energyResidualRelRms',NaN,'momentumResidualRms',NaN); end
function tf = should_apply_mass_regularization(step, opts)
tf = false;
if ~opts.massRegularizationEnable
    return;
end
steps = opts.massRegularizationSteps;
if ~isempty(steps) && any(step == steps)
    tf = true;
    return;
end
if isempty(steps) && opts.massRegularizationEvery > 0 && mod(step, opts.massRegularizationEvery) == 0
    tf = true;
end
end
function v = get_field(s,name,defaultValue)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), v=s.(name); else, v=defaultValue; end
end

function opts = parse_options(varargin)
opts=struct(); opts.outputRoot=fullfile('runs','resamp_poiseuille_wallvp_validation'); opts.cases={'classic_wallvp_reference','q6_resampled_wallvp'}; opts.massRegularizationEnable=false; opts.massRegularizationSteps=[]; opts.massRegularizationEvery=0; opts.massRegularizationStrength=1.0; opts.massRegularizationTriggerMode='always'; opts.massRegularizationRelStdTrigger=0.20; opts.massRegularizationMinFactor=0.25; opts.massRegularizationMaxFactor=4.0; opts.massRegularizationPreserveEnergy=true; opts.massRegularizationMinParticlesPerCell=2; opts.Nx=64; opts.Ny=32; opts.gamma=20; opts.NTarget=20; opts.NMin=14; opts.NMax=26; opts.steps=5000; opts.sampleEvery=50; opts.summaryEvery=100; opts.dt=0.005; opts.alphaDeg=90; opts.kBT=0.01; opts.particleMass=1.0; opts.bodyForceX=0.005; opts.wallModeY='bounceback'; opts.capacityFactor=2.0; opts.thermostatStrength=0.25; opts.projectionInterpolationMethod='nearest'; opts.extractSelectionMode='closest_to_cell_mean'; opts.insertVelocityMode='current_or_memory_pairwise'; opts.memoryMinParticles=14; opts.remapMethod='scale_preserve_velocity'; opts.massMinFactor=0.05; opts.massMaxFactor=20.0; opts.remapMassSafetyMode='uniform_mass_velocity_shift'; opts.massSafetyMinFactor=0.25; opts.massSafetyMaxFactor=4.0; opts.constraintTolerance=1e-10; opts.preservePreEditVelocity=true; opts.alwaysUsePreEditVelocityForRemap=false; opts.initialPoiseuilleProfileEnable=false; opts.initialPoiseuilleNuGuess=0.05; opts.initialPoiseuilleScale=1.0; opts.excludeWallCells=2; opts.fitStartFraction=0.5; opts.poiseuilleFitModel='slip'; opts.rngSeed=12345; opts.visualEvery=0; opts.figureId=720; opts.profileFigureId=721; opts.showProfileFigure=true; opts.saveFrames=false; opts.saveFrameEvery=100; opts.saveFinalFigures=true; opts.frameDirName='frames'; opts.particleMarkerSize=3; opts.pngResolution=150; opts.wallVirtualParticlesEnable=true; opts.wallVirtualParticlesGeometryMode='shifted_solid_fraction'; opts.wallVirtualParticlesForceRandomShiftY=true; opts.wallVirtualParticlesDensityFactor=1.0; opts.wallVirtualParticlesPerCell=[]; opts.wallVirtualParticleMass=1.0; opts.wallVirtualParticlesThermal=true; opts.wallVirtualParticlesKBT=0.01; opts.wallVirtualParticlesStochasticCount=true;
if mod(numel(varargin),2)~=0, error('Options must be name/value pairs.'); end
for k=1:2:numel(varargin)
    key=lower(char(string(varargin{k}))); val=varargin{k+1};
    switch key
        case 'outputroot', opts.outputRoot=char(string(val));
        case 'cases', opts.cases=val;
        case 'massregularizationenable', opts.massRegularizationEnable=logical(val);
        case {'massregularizationsteps','massregularizationstep'}, opts.massRegularizationSteps=val;
        case 'massregularizationevery', opts.massRegularizationEvery=val;
        case 'massregularizationstrength', opts.massRegularizationStrength=val;
        case 'massregularizationtriggermode', opts.massRegularizationTriggerMode=lower(char(string(val)));
        case 'massregularizationrelstdtrigger', opts.massRegularizationRelStdTrigger=val;
        case 'massregularizationminfactor', opts.massRegularizationMinFactor=val;
        case 'massregularizationmaxfactor', opts.massRegularizationMaxFactor=val;
        case 'massregularizationpreserveenergy', opts.massRegularizationPreserveEnergy=logical(val);
        case 'massregularizationminparticlespercell', opts.massRegularizationMinParticlesPerCell=val;
        case 'nx', opts.Nx=val; case 'ny', opts.Ny=val; case 'gamma', opts.gamma=val; opts.NTarget=val; case 'ntarget', opts.NTarget=val; case 'nmin', opts.NMin=val; case 'nmax', opts.NMax=val; case 'steps', opts.steps=val; case 'sampleevery', opts.sampleEvery=val; case 'summaryevery', opts.summaryEvery=val; case 'dt', opts.dt=val; case 'alphadeg', opts.alphaDeg=val; case 'kbt', opts.kBT=val; opts.wallVirtualParticlesKBT=val; case 'particlemass', opts.particleMass=val; opts.wallVirtualParticleMass=val; case 'bodyforcex', opts.bodyForceX=val; case 'wallmodey', opts.wallModeY=char(string(val)); case 'capacityfactor', opts.capacityFactor=val; case 'thermostatstrength', opts.thermostatStrength=val; case 'projectioninterpolationmethod', opts.projectionInterpolationMethod=char(string(val)); case 'extractselectionmode', opts.extractSelectionMode=lower(char(string(val))); case 'insertvelocitymode', opts.insertVelocityMode=lower(char(string(val))); case 'memoryminparticles', opts.memoryMinParticles=val; case 'remapmethod', opts.remapMethod=lower(char(string(val))); case 'massminfactor', opts.massMinFactor=val; case 'massmaxfactor', opts.massMaxFactor=val; case 'remapmasssafetymode', opts.remapMassSafetyMode=lower(char(string(val))); case 'masssafetyminfactor', opts.massSafetyMinFactor=val; case 'masssafetymaxfactor', opts.massSafetyMaxFactor=val; case 'constrainttolerance', opts.constraintTolerance=val; case 'preservepreeditvelocity', opts.preservePreEditVelocity=logical(val); case 'alwaysusepreeditvelocityforremap', opts.alwaysUsePreEditVelocityForRemap=logical(val); case 'initialpoiseuilleprofileenable', opts.initialPoiseuilleProfileEnable=logical(val); case 'initialpoiseuillenuguess', opts.initialPoiseuilleNuGuess=val; case 'initialpoiseuillescale', opts.initialPoiseuilleScale=val; case 'excludewallcells', opts.excludeWallCells=val; case 'fitstartfraction', opts.fitStartFraction=val; case 'poiseuillefitmodel', opts.poiseuilleFitModel=lower(char(string(val))); case 'visualevery', opts.visualEvery=val; case 'figureid', opts.figureId=val; case 'profilefigureid', opts.profileFigureId=val; case 'showprofilefigure', opts.showProfileFigure=logical(val); case 'saveframes', opts.saveFrames=logical(val); case 'saveframeevery', opts.saveFrameEvery=val; case 'savefinalfigures', opts.saveFinalFigures=logical(val); case 'framedirname', opts.frameDirName=char(string(val)); case 'particlemarkersize', opts.particleMarkerSize=val; case 'pngresolution', opts.pngResolution=val; case 'rngseed', opts.rngSeed=val;
        case 'wallvirtualparticlesenable', opts.wallVirtualParticlesEnable=logical(val); case 'wallvirtualparticlesgeometrymode', opts.wallVirtualParticlesGeometryMode=lower(char(string(val))); case 'wallvirtualparticlesforcerandomshifty', opts.wallVirtualParticlesForceRandomShiftY=logical(val); case 'wallvirtualparticlesdensityfactor', opts.wallVirtualParticlesDensityFactor=val; case 'wallvirtualparticlespercell', opts.wallVirtualParticlesPerCell=val; case 'wallvirtualparticlemass', opts.wallVirtualParticleMass=val; case 'wallvirtualparticlesthermal', opts.wallVirtualParticlesThermal=logical(val); case 'wallvirtualparticleskbt', opts.wallVirtualParticlesKBT=val; case 'wallvirtualparticlesstochasticcount', opts.wallVirtualParticlesStochasticCount=logical(val);
        otherwise, error('Unknown option: %s', key);
    end
end
if ischar(opts.massRegularizationSteps) || isstring(opts.massRegularizationSteps)
    if strlength(string(opts.massRegularizationSteps)) == 0
        opts.massRegularizationSteps = [];
    else
        opts.massRegularizationSteps = str2num(char(opts.massRegularizationSteps)); %#ok<ST2NM>
    end
end
opts.massRegularizationSteps = opts.massRegularizationSteps(:)';
end
