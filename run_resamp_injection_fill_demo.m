function out = run_resamp_injection_fill_demo(varargin)
%RUN_RESAMP_INJECTION_FILL_DEMO Latent-domain injection/fill test.
%
% Starts from gamma latent particles per cell (m=0, v=0, non-fluid), then
% activates gamma fluid particles per step in a source cell.  Wet cells evolve
% dynamically from deposited fluid mass.  Mass/remap/recycling act only on wet
% cells.
%
% Example:
%   out = run_resamp_injection_fill_demo('steps',200,'method','weighted_classic');

opts = parse_options(varargin{:});
if ~isempty(opts.rngSeed)
    rng(opts.rngSeed, 'twister');
end
params = default_params(opts);
if ~exist(opts.outputDir, 'dir'), mkdir(opts.outputDir); end
if opts.saveFrames && ~exist(opts.frameDir, 'dir'), mkdir(opts.frameDir); end

% Build a full latent marker grid: gamma storage slots per cell, but no fluid.

[state, initInfo] = resamp_initialize_particles_taylor_green_forced(params);
[state, poolInfo0] = resamp_enable_particle_pool(state, 'capacityFactor', opts.capacityFactor);
state = make_all_slots_latent(state);
state.cellWetMask = false(params.Nx, params.Ny);
params.cellWetMask = state.cellWetMask;
state.uMemUx = zeros(params.Nx, params.Ny);
state.uMemUy = zeros(params.Nx, params.Ny);
state.uMemValid = false(params.Nx, params.Ny);

rows = repmat(empty_row(), floor(opts.steps/max(opts.summaryEvery,1))+3, 1);
irow = 1;
insertStats = empty_cumulative_stats();
extractStats = empty_cumulative_stats();
activateStats = empty_cumulative_stats();
lastInjectDiag = empty_activation_diag();
lastInsertDiag = empty_insert_diag();
lastExtractDiag = empty_extract_diag();
lastRemapDiag = empty_remap_diag();
lastStepDiag = empty_step_diag();
lastTherm = empty_thermostat_diag();

[md, tg] = diagnostics(state, params);
rows(irow) = make_row(0, md, tg, lastInjectDiag, lastExtractDiag, lastInsertDiag, lastRemapDiag, lastTherm);
visualize_if_needed(state, params, opts, 0, lastInsertDiag, lastExtractDiag, lastRemapDiag, lastStepDiag);

for step = 1:opts.steps
    % Injection source: activate nominal support in one or several cells.
    if opts.injectionEvery > 0 && mod(step, opts.injectionEvery) == 0
        sourceCells = source_cell_ids(params, opts);
        [state, injectDiag] = resamp_activate_particles_in_cell(state, params, sourceCells, ...
            'nPerCell', opts.injectionParticlesPerCell, ...
            'particleMass', opts.particleMass, ...
            'velocity', [opts.inletVelocityX, opts.inletVelocityY], ...
            'kBT', opts.injectionKBT, ...
            'preferLatent', true);
        state.cellWetMask(sourceCells_to_sub(params, sourceCells)) = true;
        params.cellWetMask = state.cellWetMask;
    else
        injectDiag = empty_activation_diag();
    end
    activateStats = update_cumulative(activateStats, get_field(injectDiag,'nActivatedParticles',0), step);
    injectDiag = attach_activation_cumulative(injectDiag, activateStats);
    lastInjectDiag = injectDiag;

    % Evolve only fluid-active particles.  Latent markers are excluded by resamp_active_mask.
    switch opts.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_periodic_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag] = resamp_step_projection_periodic_weighted(state, params);
        otherwise
            error('Unknown method: %s', opts.method);
    end
    lastStepDiag = stepDiag;

    % Wet mask evolves from current deposited fluid mass.
    if opts.wetMaskUpdateEvery > 0 && mod(step, opts.wetMaskUpdateEvery) == 0
        [state, wetDiag] = resamp_update_cell_wet_mask(state, params, ...
            'mode', opts.wetMaskUpdateMode, ...
            'wetMassOnThreshold', opts.wetMassOnThreshold, ...
            'wetMassOffThreshold', opts.wetMassOffThreshold, ...
            'wetParticleOnThreshold', opts.wetParticleOnThreshold, ...
            'wetParticleOffThreshold', opts.wetParticleOffThreshold);
        params.cellWetMask = state.cellWetMask;
        if opts.latentizeDryCellsAfterWetUpdate
            [state, latentDiag] = resamp_set_dry_cell_particles_latent(state, params, 'cellWetMask', state.cellWetMask);
        else
            latentDiag = struct('nConvertedToLatent',0);
        end
    else
        wetDiag = struct('nWetCells', nnz(state.cellWetMask), 'nDryCells', numel(state.cellWetMask)-nnz(state.cellWetMask));
        latentDiag = struct('nConvertedToLatent',0);
    end

    % Save target field before population edits.
    preEditG = [];
    if opts.preservePreEditVelocity && opts.remapEvery > 0 && mod(step, opts.remapEvery) == 0
        preEditG = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
            'periodicX', true, 'periodicY', true, 'activeMask', resamp_active_mask(state), 'cellWetMask', state.cellWetMask);
    end

    if opts.extractEvery > 0 && mod(step, opts.extractEvery) == 0
        [state, extractDiag] = resamp_extract_overpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, 'NMin', opts.NMin, 'NMax', opts.NMax, ...
            'extractSelectionMode', opts.extractSelectionMode, 'computeDiagnostics', true);
    else
        extractDiag = empty_extract_diag();
    end
    extractStats = update_cumulative(extractStats, get_field(extractDiag,'nExtractedParticles',0), step);
    extractDiag.nExtractedParticlesCumulative = extractStats.cumulative;
    extractDiag.lastExtractionStep = extractStats.lastStep;
    lastExtractDiag = extractDiag;

    if opts.insertEvery > 0 && mod(step, opts.insertEvery) == 0
        [state, insertDiag] = resamp_insert_underpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, 'NMin', opts.NMin, 'NMax', opts.NMax, ...
            'memoryMinParticles', opts.memoryMinParticles, 'particleMass', opts.particleMass, ...
            'kBT', opts.insertKBT, 'insertVelocityMode', opts.insertVelocityMode, 'computeDiagnostics', true);
    else
        insertDiag = empty_insert_diag();
    end
    insertStats = update_cumulative(insertStats, get_field(insertDiag,'nInsertedParticles',0), step);
    insertDiag.nInsertedParticlesCumulative = insertStats.cumulative;
    insertDiag.lastInsertionStep = insertStats.lastStep;
    lastInsertDiag = insertDiag;

    if opts.remapEvery > 0 && mod(step, opts.remapEvery) == 0
        remapArgs = {'targetCellMass', params.resampTargetCellMass, 'method', opts.remapMethod, ...
            'massSafetyEnable', opts.remapMassSafetyEnable, 'massSafetyMode', opts.remapMassSafetyMode, ...
            'massSafetyMinFactor', opts.remapMassSafetyMinFactor, 'massSafetyMaxFactor', opts.remapMassSafetyMaxFactor, ...
            'massMin', opts.massMinFactor * opts.particleMass, 'massMax', opts.massMaxFactor * opts.particleMass, ...
            'constraintTolerance', opts.constraintTolerance, 'computeDiagnostics', true};
        didEdit = get_field(extractDiag,'nExtractedParticles',0)>0 || get_field(insertDiag,'nInsertedParticles',0)>0 || get_field(injectDiag,'nActivatedParticles',0)>0;
        if opts.preservePreEditVelocity && didEdit && ~isempty(preEditG)
            remapArgs = [remapArgs, {'targetVelocityMode','grid','targetUx',preEditG.Ux,'targetUy',preEditG.Uy}]; %#ok<AGROW>
        else
            remapArgs = [remapArgs, {'targetVelocityMode', opts.targetVelocityMode}]; %#ok<AGROW>
        end
        [state, remapDiag] = resamp_local_mass_moment_remap(state, params, remapArgs{:});
    else
        remapDiag = empty_remap_diag();
    end
    lastRemapDiag = remapDiag;

    if opts.thermostatAfterRemap
        [state.v, therm] = resamp_apply_cell_thermostat_weighted(state.x, state.v, state.m, params, ...
            'periodicX', true, 'periodicY', true, 'activeMask', resamp_active_mask(state), ...
            'targetKBT', opts.thermostatTargetKBT, 'strength', opts.thermostatStrength);
    else
        therm = empty_thermostat_diag();
    end
    therm.nConvertedDryToLatent = get_field(latentDiag, 'nConvertedToLatent', 0);
    lastTherm = therm;

    doSummary = (opts.summaryEvery > 0 && mod(step, opts.summaryEvery)==0) || step == opts.steps;
    doVisual = (opts.visualEvery > 0 && mod(step, opts.visualEvery)==0) || step == opts.steps;
    if doSummary || doVisual
        [md, tg] = diagnostics(state, params);
    end
    if doSummary
        irow = irow + 1;
        rows(irow) = make_row(step, md, tg, injectDiag, extractDiag, insertDiag, remapDiag, therm);
        fprintf('fill step=%5d wet=%4d/%4d fluid=%6d latent=%6d injCum=%6g extCum=%6g insCum=%6g N[min,max]=[%g,%g] Mrel=%.3e kBT=%.4g\n', ...
            step, md.nWetCells, md.NCells, md.NpActive, get_role_count(state,'latent'), ...
            get_field(injectDiag,'nActivatedParticlesCumulative',0), get_field(extractDiag,'nExtractedParticlesCumulative',0), ...
            get_field(insertDiag,'nInsertedParticlesCumulative',0), md.NMin, md.NMax, md.MRelRms, md.kBTWeighted);
    end
    if doVisual
        visualize_if_needed(state, params, opts, step, insertDiag, extractDiag, remapDiag, stepDiag);
    end
end

summary = struct2table(rows(1:irow));
out = struct();
out.params = params;
out.options = opts;
out.initialInfo = initInfo;
out.initialPoolInfo = poolInfo0;
out.state = state;
out.summary = summary;
out.finalMassDiagnostics = md;
out.finalTaylorGreenDiagnostics = tg;
out.lastInjectionDiag = lastInjectDiag;
out.lastExtractDiag = lastExtractDiag;
out.lastInsertDiag = lastInsertDiag;
out.lastRemapDiag = lastRemapDiag;
out.lastThermostatAfterRemapDiag = lastTherm;
if opts.writeCsv
    csvPath = fullfile(opts.outputDir, sprintf('resamp_injection_fill_%s.csv', opts.method));
    writetable(summary, csvPath);
    out.csvPath = csvPath;
    fprintf('Wrote %s\n', csvPath);
end
end

function params = default_params(opts)
params = struct();
params.Lx = opts.Nx; params.Ly = opts.Ny; params.Nx = opts.Nx; params.Ny = opts.Ny;
params.gamma = opts.gamma; params.dt = opts.dt; params.alphaDeg = opts.alphaDeg;
params.kBT = opts.kBT; params.initialPopulationMode = 'exact_per_cell';
params.taylorGreenAmplitude = 0; params.taylorGreenInitialAmplitude = 0;
params.taylorGreenThermalNoise = false; params.taylorGreenForcingEnable = false;
params.bodyForceX = 0; params.bodyForceY = 0; params.useRandomGridShift = true;
params.projectionEnable = true; params.projectionStrength = opts.projectionStrength;
params.projectionInterpolationMethod = opts.projectionInterpolationMethod;
params.projectionMomentumCorrectionEnable = true; params.projectionMomentumCorrectionMode = 'particle_global_exact';
params.thermostatAfterStep = false; params.thermostatAfterProjection = false;
params.thermostatTargetKBT = opts.thermostatTargetKBT; params.thermostatStrength = opts.thermostatStrength;
params.resampParticleMass = opts.particleMass; params.resampTargetCellMass = opts.NTarget * opts.particleMass;
params.resampNTarget = opts.NTarget; params.resampNMin = opts.NMin; params.resampNMax = opts.NMax;
params.resampMemoryMinParticles = opts.memoryMinParticles;
params.resampWetMassOnThreshold = opts.wetMassOnThreshold; params.resampWetMassOffThreshold = opts.wetMassOffThreshold;
params.resampWetParticleOnThreshold = opts.wetParticleOnThreshold; params.resampWetParticleOffThreshold = opts.wetParticleOffThreshold;
params.visualMaxParticles = opts.visualMaxParticles; params.visualQuiverScale = opts.visualQuiverScale;
params.computeDiagnostics = true;
end

function opts = parse_options(varargin)
opts = struct();
opts.method = 'weighted_classic';
opts.steps = 200; opts.summaryEvery = 10; opts.visualEvery = 5;
opts.Nx = 64; opts.Ny = 32; opts.gamma = 20; opts.dt = 0.01; opts.alphaDeg = 90;
opts.kBT = 0; opts.particleMass = 1.0; opts.NTarget = []; opts.NMin = []; opts.NMax = [];
opts.capacityFactor = 2.0; opts.rngSeed = 12345;
opts.injectionCell = []; opts.injectionPatchSize = [1 1]; opts.injectionEvery = 1;
opts.injectionParticlesPerCell = []; opts.inletVelocityX = 0.15; opts.inletVelocityY = 0.0; opts.injectionKBT = 0.01;
opts.wetMaskUpdateEvery = 1; opts.wetMaskUpdateMode = 'mass_hysteresis';
opts.wetMassOnThreshold = []; opts.wetMassOffThreshold = []; opts.wetParticleOnThreshold = []; opts.wetParticleOffThreshold = [];
opts.latentizeDryCellsAfterWetUpdate = true;
opts.extractEvery = 1; opts.insertEvery = 1; opts.remapEvery = 1; opts.extractSelectionMode = 'closest_to_cell_mean';
opts.preservePreEditVelocity = true; opts.memoryMinParticles = []; opts.insertVelocityMode = 'current_or_memory_pairwise'; opts.insertKBT = [];
opts.remapMethod = 'scale_preserve_velocity'; opts.targetVelocityMode = 'preserve_cell_velocity'; opts.massMinFactor = 0.05; opts.massMaxFactor = 20.0;
opts.remapMassSafetyEnable = true; opts.remapMassSafetyMode = 'uniform_mass_velocity_shift'; opts.remapMassSafetyMinFactor = 0.25; opts.remapMassSafetyMaxFactor = 4.0; opts.constraintTolerance = 1e-10;
opts.projectionStrength = 1.0; opts.projectionInterpolationMethod = 'nearest';
opts.thermostatAfterRemap = true; opts.thermostatStrength = 0.25; opts.thermostatTargetKBT = 0.01;
opts.figureId = 720; opts.debugFigureId = 721; opts.visualMaxParticles = 12000; opts.visualQuiverScale = 1.2;
opts.saveFrames = false; opts.frameDir = fullfile('runs','resamp_injection_fill','frames'); opts.framePrefix = 'resamp_fill';
opts.outputDir = fullfile('runs','resamp_injection_fill'); opts.writeCsv = true;
if mod(numel(varargin),2) ~= 0, error('Options must be name/value pairs.'); end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k}))); val = varargin{k+1};
    switch key
        case 'method', opts.method = lower(strrep(char(string(val)),'-','_'));
        case 'steps', opts.steps = val;
        case 'summaryevery', opts.summaryEvery = val;
        case 'visualevery', opts.visualEvery = val;
        case 'nx', opts.Nx = val;
        case 'ny', opts.Ny = val;
        case 'gamma', opts.gamma = val;
        case 'dt', opts.dt = val;
        case 'alphadeg', opts.alphaDeg = val;
        case 'particlemass', opts.particleMass = val;
        case 'capacityfactor', opts.capacityFactor = val;
        case 'rngseed', opts.rngSeed = val;
        case 'ntarget', opts.NTarget = val;
        case 'nmin', opts.NMin = val;
        case 'nmax', opts.NMax = val;
        case {'injectioncell','sourcecell'}, opts.injectionCell = val;
        case 'injectionpatchsize', opts.injectionPatchSize = val;
        case 'injectionevery', opts.injectionEvery = val;
        case {'injectionparticlespercell','injectionn'}, opts.injectionParticlesPerCell = val;
        case {'inletvelocityx','uin','uinx'}, opts.inletVelocityX = val;
        case {'inletvelocityy','uiny'}, opts.inletVelocityY = val;
        case 'injectionkbt', opts.injectionKBT = val;
        case 'wetmaskupdateevery', opts.wetMaskUpdateEvery = val;
        case 'wetmaskupdatemode', opts.wetMaskUpdateMode = lower(char(string(val)));
        case 'wetmassonthreshold', opts.wetMassOnThreshold = val;
        case 'wetmassoffthreshold', opts.wetMassOffThreshold = val;
        case 'wetparticleonthreshold', opts.wetParticleOnThreshold = val;
        case 'wetparticleoffthreshold', opts.wetParticleOffThreshold = val;
        case 'latentizedrycellsafterwetupdate', opts.latentizeDryCellsAfterWetUpdate = logical(val);
        case 'extractevery', opts.extractEvery = val;
        case 'insertevery', opts.insertEvery = val;
        case 'remapevery', opts.remapEvery = val;
        case 'extractselectionmode', opts.extractSelectionMode = lower(char(string(val)));
        case 'preservepreeditvelocity', opts.preservePreEditVelocity = logical(val);
        case 'memoryminparticles', opts.memoryMinParticles = val;
        case 'insertvelocitymode', opts.insertVelocityMode = lower(char(string(val)));
        case 'insertkbt', opts.insertKBT = val;
        case 'remapmethod', opts.remapMethod = lower(char(string(val)));
        case 'targetvelocitymode', opts.targetVelocityMode = lower(char(string(val)));
        case 'massminfactor', opts.massMinFactor = val;
        case 'massmaxfactor', opts.massMaxFactor = val;
        case 'remapmasssafetyenable', opts.remapMassSafetyEnable = logical(val);
        case 'remapmasssafetymode', opts.remapMassSafetyMode = lower(char(string(val)));
        case 'remapmasssafetyminfactor', opts.remapMassSafetyMinFactor = val;
        case 'remapmasssafetymaxfactor', opts.remapMassSafetyMaxFactor = val;
        case 'thermostatafterremap', opts.thermostatAfterRemap = logical(val);
        case 'thermostatstrength', opts.thermostatStrength = val;
        case 'thermostattargetkbt', opts.thermostatTargetKBT = val;
        case 'projectionstrength', opts.projectionStrength = val;
        case 'projectioninterpolationmethod', opts.projectionInterpolationMethod = char(string(val));
        case 'figureid', opts.figureId = val;
        case 'visualmaxparticles', opts.visualMaxParticles = val;
        case 'visualquiverscale', opts.visualQuiverScale = val;
        case 'saveframes', opts.saveFrames = logical(val);
        case 'framedir', opts.frameDir = char(string(val));
        case 'frameprefix', opts.framePrefix = char(string(val));
        case 'outputdir', opts.outputDir = char(string(val));
        case 'writecsv', opts.writeCsv = logical(val);
        otherwise, error('Unknown option: %s', key);
    end
end
if isempty(opts.NTarget), opts.NTarget = opts.gamma; end
if isempty(opts.NMin), opts.NMin = ceil(0.5*opts.gamma); end
if isempty(opts.NMax), opts.NMax = ceil(1.5*opts.gamma); end
if isempty(opts.memoryMinParticles), opts.memoryMinParticles = opts.NMin; end
if isempty(opts.injectionParticlesPerCell), opts.injectionParticlesPerCell = opts.gamma; end
if isempty(opts.insertKBT), opts.insertKBT = opts.injectionKBT; end
if isempty(opts.wetMassOnThreshold), opts.wetMassOnThreshold = 0.05 * opts.NTarget * opts.particleMass; end
if isempty(opts.wetMassOffThreshold), opts.wetMassOffThreshold = 0.005 * opts.NTarget * opts.particleMass; end
if isempty(opts.wetParticleOnThreshold), opts.wetParticleOnThreshold = 1; end
if isempty(opts.wetParticleOffThreshold), opts.wetParticleOffThreshold = 0; end
end

function state = make_all_slots_latent(state)
Np = size(state.x,1);
state.active = false(Np,1);
state.particleRole = zeros(Np,1);
state.particleRole(1:get_field(state,'Ninitial',Np)) = 2;
state.m(:) = 0;
state.v(:,:) = 0;
state.Nactive = 0;
state.Ncapacity = Np;
end

function ids = source_cell_ids(params, opts)
Nx=params.Nx; Ny=params.Ny;
if isempty(opts.injectionCell)
    cx = max(1, round(0.1*Nx)); cy = max(1, round(0.5*Ny));
else
    cx = opts.injectionCell(1); cy = opts.injectionCell(2);
end
sx = max(1, round(opts.injectionPatchSize(1))); sy = max(1, round(opts.injectionPatchSize(2)));
% Build exactly sx-by-sy periodic cell lists.  The previous expression
% placed the final '-1' on the colon upper bound and returned an empty list
% for sx=1/sy=1, which silently disabled injection.
ix0 = cx - floor((sx - 1) / 2);
iy0 = cy - floor((sy - 1) / 2);
ixList = mod((ix0:(ix0 + sx - 1)) - 1, Nx) + 1;
iyList = mod((iy0:(iy0 + sy - 1)) - 1, Ny) + 1;
[IX,IY]=ndgrid(ixList,iyList);
ids=(IX(:)-1)*Ny+IY(:);
ids=unique(ids(:),'stable');
end

function mask = sourceCells_to_sub(params, ids)
mask = false(params.Nx, params.Ny);
Ny=params.Ny;
for k=1:numel(ids)
    ix=floor((ids(k)-1)/Ny)+1; iy=ids(k)-Ny*(ix-1);
    mask(ix,iy)=true;
end
end

function visualize_if_needed(state, params, opts, step, insertDiag, extractDiag, remapDiag, stepDiag)
if opts.visualEvery <= 0 && step ~= opts.steps, return; end
if step ~= 0 && mod(step, opts.visualEvery) ~= 0 && step ~= opts.steps, return; end
resamp_pool_visualize_frame(state, params, step, 'insertDiag', insertDiag, 'extractDiag', extractDiag, ...
    'remapDiag', remapDiag, 'stepDiag', stepDiag, 'figureId', opts.figureId, ...
    'saveFrame', opts.saveFrames, 'frameDir', opts.frameDir, 'framePrefix', opts.framePrefix, ...
    'titleSuffix', 'latent injection fill');
end

function [md,tg]=diagnostics(state,params)
md=resamp_population_mass_diagnostics(state,params,'periodicX',true,'periodicY',true,'cellWetMask',state.cellWetMask);
md.dt = params.dt;
if exist('projection_taylor_green_diagnostics','file')==2
    tg=projection_taylor_green_diagnostics(md.G,params);
else
    tg=struct('modeAmplitude',NaN,'modeCoherence',NaN,'enstrophy',NaN);
end
end

function row = make_row(step, md, tg, injectDiag, extractDiag, insertDiag, remapDiag, therm)
roleCounts = get_role_count(md.G, 'none'); %#ok<NASGU>
row=empty_row();
row.step=step; row.t=step * get_field(md, 'dt', NaN);
row.NpActive=md.NpActive; row.Ncapacity=md.Ncapacity; row.Nfree=md.Nfree;
row.nWetCells=md.nWetCells; row.nDryCells=md.nDryCells; row.wetFraction=md.wetFraction;
row.NMin=md.NMin; row.NMax=md.NMax; row.NStd=md.NStd; row.MRelRms=md.MRelRms;
row.mParticleRelStd=md.mParticleRelStd; row.mParticleMin=md.mParticleMin; row.mParticleMax=md.mParticleMax;
row.kBTWeighted=md.kBTWeighted; row.rmsDivParticleAfter=get_field(tg,'rmsDivParticleAfter',NaN);
row.activatedParticles=get_field(injectDiag,'nActivatedParticles',NaN); row.activatedParticlesCumulative=get_field(injectDiag,'nActivatedParticlesCumulative',NaN);
row.extractedParticles=get_field(extractDiag,'nExtractedParticles',NaN); row.extractedParticlesCumulative=get_field(extractDiag,'nExtractedParticlesCumulative',NaN);
row.insertedParticles=get_field(insertDiag,'nInsertedParticles',NaN); row.insertedParticlesCumulative=get_field(insertDiag,'nInsertedParticlesCumulative',NaN);
row.remapMassSafetyCells=get_field(remapDiag,'nCellsMassSafetyApplied',NaN);
row.thermostatKBTAfter=get_field(therm,'meanKBTAfter',NaN);
row.nConvertedDryToLatent=get_field(therm,'nConvertedDryToLatent',NaN);
end

function row=empty_row()
row=struct('step',NaN,'t',NaN,'NpActive',NaN,'Ncapacity',NaN,'Nfree',NaN,'nWetCells',NaN,'nDryCells',NaN,'wetFraction',NaN, ...
    'NMin',NaN,'NMax',NaN,'NStd',NaN,'MRelRms',NaN,'mParticleRelStd',NaN,'mParticleMin',NaN,'mParticleMax',NaN, ...
    'kBTWeighted',NaN,'rmsDivParticleAfter',NaN,'activatedParticles',NaN,'activatedParticlesCumulative',NaN, ...
    'extractedParticles',NaN,'extractedParticlesCumulative',NaN,'insertedParticles',NaN,'insertedParticlesCumulative',NaN, ...
    'remapMassSafetyCells',NaN,'thermostatKBTAfter',NaN,'nConvertedDryToLatent',NaN);
end

function s=empty_activation_diag(), s=struct('nActivatedParticles',0,'nActivatedParticlesCumulative',0,'lastActivationStep',NaN); end
function s=empty_insert_diag(), s=struct('nInsertedParticles',0,'nInsertedParticlesCumulative',0,'lastInsertionStep',NaN); end
function s=empty_extract_diag(), s=struct('nExtractedParticles',0,'nExtractedParticlesCumulative',0,'lastExtractionStep',NaN); end
function s=empty_remap_diag(), s=struct('nCellsMassSafetyApplied',NaN); end
function s=empty_step_diag(), s=struct(); end
function s=empty_thermostat_diag(), s=struct('meanKBTAfter',NaN,'nConvertedDryToLatent',NaN); end
function st=empty_cumulative_stats(), st=struct('cumulative',0,'lastStep',NaN); end
function st=update_cumulative(st,n,step), if ~isfinite(n), n=0; end, if n>0, st.cumulative=st.cumulative+n; st.lastStep=step; end, end
function d=attach_activation_cumulative(d,st), d.nActivatedParticlesCumulative=st.cumulative; d.lastActivationStep=st.lastStep; end
function n=get_role_count(stateOrG,what), n=NaN; %#ok<INUSD>
end
function v=get_field(s,name,default), if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=default; end, end
