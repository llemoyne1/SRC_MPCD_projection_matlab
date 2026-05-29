function out = run_resamp_pool_insertion_visual_demo(varargin)
%RUN_RESAMP_POOL_INSERTION_VISUAL_DEMO Live visual demo for weighted pool insertion.
%
% Examples:
%   out = run_resamp_pool_insertion_visual_demo('method','weighted_q6');
%   out = run_resamp_pool_insertion_visual_demo('method','weighted_classic','initialDepletion','patch');
%   out = run_resamp_pool_insertion_visual_demo('showDebugFigure',true,'saveFrames',true);

opts = parse_options(varargin{:});
if ~isempty(opts.rngSeed)
    rng(opts.rngSeed, 'twister');
end

params = default_params(opts);
params.method = opts.method;
params.visualFigureId = opts.figureId;
params.visualFigureName = sprintf('resamp pool %s', opts.method);
params.visualMaxParticles = opts.visualMaxParticles;

[state, initInfo] = resamp_initialize_particles_taylor_green_forced(params);
[state, poolInfo0] = resamp_enable_particle_pool(state, 'capacityFactor', opts.capacityFactor);
state = initialize_velocity_memory(state, params);
[state, initialPopulationEditInfo] = apply_initial_depletion(state, params, opts);
[state, initialWetMaskInfo] = resamp_update_cell_wet_mask(state, params, ...
    'mode', opts.initialWetMaskMode, ...
    'cellWetMask', opts.cellWetMask, ...
    'wetMassOnThreshold', opts.wetMassOnThreshold, ...
    'wetMassOffThreshold', opts.wetMassOffThreshold, ...
    'wetParticleOnThreshold', opts.wetParticleOnThreshold, ...
    'wetParticleOffThreshold', opts.wetParticleOffThreshold);
params.cellWetMask = state.cellWetMask;

if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir);
end
if opts.saveFrames && ~exist(opts.frameDir, 'dir')
    mkdir(opts.frameDir);
end

nRows = floor(opts.steps / opts.summaryEvery) + 2;
rows = repmat(empty_row(), nRows, 1);
irow = 1;
insertStats = empty_insert_stats();
extractStats = empty_extract_stats();
thermostatAfterRemapDiag = empty_thermostat_after_remap_diag();
insertDiag0 = attach_insert_cumulative(empty_insert_diag(), insertStats);
extractDiag0 = attach_extract_cumulative(empty_extract_diag(), extractStats);
stepDiag0 = attach_thermostat_after_remap(empty_step_diag(), thermostatAfterRemapDiag);
[md, tg] = diagnostics(state, params);
rows(irow) = make_row(0, 0.0, md, tg, stepDiag0, insertDiag0, extractDiag0, empty_remap_diag());

resamp_pool_visualize_frame(state, params, 0, ...
    'insertDiag', insertDiag0, ...
    'extractDiag', extractDiag0, ...
    'remapDiag', empty_remap_diag(), ...
    'stepDiag', stepDiag0, ...
    'figureId', opts.figureId, ...
    'showDebugFigure', opts.showDebugFigure, ...
    'debugFigureId', opts.debugFigureId, ...
    'saveFrame', opts.saveFrames, ...
    'frameDir', opts.frameDir, ...
    'framePrefix', opts.framePrefix, ...
    'titleSuffix', opts.initialDepletion);

lastStepDiag = stepDiag0;
lastInsertDiag = insertDiag0;
lastExtractDiag = extractDiag0;
lastRemapDiag = empty_remap_diag();
lastThermostatAfterRemapDiag = thermostatAfterRemapDiag;

for step = 1:opts.steps
    if opts.wetMaskUpdateEvery > 0 && mod(step, opts.wetMaskUpdateEvery) == 0
        [state, wetMaskInfo] = resamp_update_cell_wet_mask(state, params, ...
            'mode', opts.wetMaskUpdateMode, ...
            'cellWetMask', opts.cellWetMask, ...
            'wetMassOnThreshold', opts.wetMassOnThreshold, ...
            'wetMassOffThreshold', opts.wetMassOffThreshold, ...
            'wetParticleOnThreshold', opts.wetParticleOnThreshold, ...
            'wetParticleOffThreshold', opts.wetParticleOffThreshold);
        params.cellWetMask = state.cellWetMask;
    else
        wetMaskInfo = empty_wet_mask_update_diag();
    end

    switch opts.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_periodic_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag] = resamp_step_projection_periodic_weighted(state, params);
        otherwise
            error('Unknown method: %s', opts.method);
    end
    preEditG = [];
    if opts.preservePreEditVelocity && opts.remapEvery > 0 && mod(step, opts.remapEvery) == 0
        activePreEdit = resamp_active_mask(state);
        preEditG = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
            'periodicX', true, 'periodicY', true, 'activeMask', activePreEdit, 'cellWetMask', state.cellWetMask);
    end

    if opts.extractEvery > 0 && mod(step, opts.extractEvery) == 0
        [state, extractDiag] = resamp_extract_overpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, ...
            'NMin', opts.NMin, ...
            'NMax', opts.NMax, ...
            'extractSelectionMode', opts.extractSelectionMode, ...
            'computeDiagnostics', true);
    else
        extractDiag = empty_extract_diag();
    end
    extractStats = update_extract_stats(extractStats, extractDiag, step);
    extractDiag = attach_extract_cumulative(extractDiag, extractStats);
    lastExtractDiag = extractDiag;

    if opts.insertEvery > 0 && mod(step, opts.insertEvery) == 0
        [state, insertDiag] = resamp_insert_underpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, ...
            'NMin', opts.NMin, ...
            'NMax', opts.NMax, ...
            'memoryMinParticles', opts.memoryMinParticles, ...
            'particleMass', opts.particleMass, ...
            'kBT', opts.insertKBT, ...
            'insertVelocityMode', opts.insertVelocityMode, ...
            'computeDiagnostics', true);
    else
        insertDiag = empty_insert_diag();
    end
    insertStats = update_insert_stats(insertStats, insertDiag, step);
    insertDiag = attach_insert_cumulative(insertDiag, insertStats);
    lastInsertDiag = insertDiag;

    if opts.remapEvery > 0 && mod(step, opts.remapEvery) == 0
        remapTargetVelocityMode = opts.targetVelocityMode;
        remapArgs = {'targetCellMass', params.resampTargetCellMass, ...
            'method', opts.remapMethod, ...
            'massMin', opts.massMinFactor * params.resampParticleMass, ...
            'massMax', opts.massMaxFactor * params.resampParticleMass, ...
            'constraintTolerance', opts.constraintTolerance, ...
            'massSafetyEnable', opts.remapMassSafetyEnable, ...
            'massSafetyMode', opts.remapMassSafetyMode, ...
            'massSafetyMinFactor', opts.remapMassSafetyMinFactor, ...
            'massSafetyMaxFactor', opts.remapMassSafetyMaxFactor, ...
            'computeDiagnostics', true};
        didPopulationEdit = get_field(extractDiag, 'nExtractedParticles', 0) > 0 || ...
            get_field(insertDiag, 'nInsertedParticles', 0) > 0;
        if opts.preservePreEditVelocity && didPopulationEdit
            if isempty(preEditG)
                activePreEdit = resamp_active_mask(state);
                preEditG = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
                    'periodicX', true, 'periodicY', true, 'activeMask', activePreEdit, 'cellWetMask', state.cellWetMask);
            end
            remapTargetVelocityMode = 'grid';
            remapArgs = [remapArgs, {'targetVelocityMode', remapTargetVelocityMode, ...
                'targetUx', preEditG.Ux, 'targetUy', preEditG.Uy}]; %#ok<AGROW>
        else
            remapArgs = [remapArgs, {'targetVelocityMode', remapTargetVelocityMode}]; %#ok<AGROW>
        end
        [state, remapDiag] = resamp_local_mass_moment_remap(state, params, remapArgs{:});
        remapDiag.populationEditDidModifyParticles = didPopulationEdit;
        remapDiag.populationEditPreservePreEditVelocity = opts.preservePreEditVelocity;
        remapDiag.populationEditTargetVelocityMode = remapTargetVelocityMode;
    else
        remapDiag = empty_remap_diag();
    end
    lastRemapDiag = remapDiag;

    if opts.thermostatAfterRemap
        [state.v, thermostatAfterRemapDiag] = resamp_apply_cell_thermostat_weighted(state.x, state.v, state.m, params, ...
            'periodicX', true, 'periodicY', true, ...
            'activeMask', resamp_active_mask(state), ...
            'targetKBT', opts.thermostatTargetKBT, ...
            'strength', opts.thermostatStrength);
    else
        thermostatAfterRemapDiag = empty_thermostat_after_remap_diag();
    end
    lastThermostatAfterRemapDiag = thermostatAfterRemapDiag;
    stepDiag = attach_thermostat_after_remap(stepDiag, thermostatAfterRemapDiag);
    lastStepDiag = stepDiag;

    doVisual = (opts.visualEvery > 0 && mod(step, opts.visualEvery) == 0) || step == opts.steps;
    doSummary = (opts.summaryEvery > 0 && mod(step, opts.summaryEvery) == 0) || step == opts.steps;
    if doSummary || doVisual
        [md, tg] = diagnostics(state, params);
    end

    if doSummary
        irow = irow + 1;
        rows(irow) = make_row(step, step * params.dt, md, tg, stepDiag, insertDiag, extractDiag, remapDiag);
        fprintf(['visual resamp %-12s step=%6d t=%.4g Nact=%7d free=%7d ', ...
                 'N[min,max]=[%3g,%3g] Mrel=%.3e insertedNow=%s insertedCum=%s lastInsert=%s ', ...
                 'poor=%4g over=%4g mRelStd=%.3e safetyCells=%s kBT=%.5g thermAfter=%.5g\n'], ...
            opts.method, step, step*params.dt, md.NpActive, md.Nfree, md.NMin, md.NMax, md.MRelRms, ...
            format_scalar(get_field(insertDiag, 'nInsertedParticles', NaN)), ...
            format_scalar(get_field(insertDiag, 'nInsertedParticlesCumulative', NaN)), ...
            format_scalar(get_field(insertDiag, 'lastInsertionStep', NaN)), ...
            get_field(insertDiag, 'nPoorCellsAfter', NaN), ...
            get_field(insertDiag, 'nOverCellsAfter', NaN), md.mParticleRelStd, ...
            format_scalar(get_field(remapDiag, 'nCellsMassSafetyApplied', NaN)), md.kBTWeighted, ...
            get_nested(stepDiag, {'thermostatAfterRemap','meanKBTAfter'}, NaN));
    end

    if doVisual
        resamp_pool_visualize_frame(state, params, step, ...
            'insertDiag', insertDiag, ...
            'extractDiag', extractDiag, ...
            'remapDiag', remapDiag, ...
            'stepDiag', stepDiag, ...
            'figureId', opts.figureId, ...
            'showDebugFigure', opts.showDebugFigure, ...
            'debugFigureId', opts.debugFigureId, ...
            'saveFrame', opts.saveFrames, ...
            'frameDir', opts.frameDir, ...
            'framePrefix', opts.framePrefix, ...
            'titleSuffix', opts.initialDepletion);
    end
end

rows = rows(1:irow);
summary = struct2table(rows);

out = struct();
out.params = params;
out.options = opts;
out.initialInfo = initInfo;
out.initialPopulationEditInfo = initialPopulationEditInfo;
out.initialPoolInfo = poolInfo0;
out.initialWetMaskInfo = initialWetMaskInfo;
out.finalWetMask = state.cellWetMask;
[~, out.finalWetMaskInfo] = resamp_cell_wet_mask(state, params, 'mode', 'auto');
out.state = state;
out.summary = summary;
out.finalMassDiagnostics = md;
out.finalTaylorGreenDiagnostics = tg;
out.finalPoolInfo = resamp_particle_pool_info(state);
out.lastStepDiag = lastStepDiag;
out.lastInsertDiag = lastInsertDiag;
out.lastExtractDiag = lastExtractDiag;
out.lastRemapDiag = lastRemapDiag;
out.lastThermostatAfterRemapDiag = lastThermostatAfterRemapDiag;
out.insertStats = insertStats;
out.extractStats = extractStats;

if opts.writeCsv
    csvPath = fullfile(opts.outputDir, sprintf('resamp_pool_insertion_visual_%s_%s.csv', opts.method, opts.initialDepletion));
    writetable(summary, csvPath);
    out.csvPath = csvPath;
    fprintf('Wrote %s\n', csvPath);
end
end

function params = default_params(opts)
params = struct();
params.Lx = opts.Nx;
params.Ly = opts.Ny;
params.Nx = opts.Nx;
params.Ny = opts.Ny;
params.gamma = opts.gamma;
params.dt = opts.dt;
params.alphaDeg = opts.alphaDeg;
params.kBT = opts.kBT;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;
params.taylorGreenAmplitude = opts.taylorGreenAmplitude;
params.taylorGreenInitialAmplitude = opts.taylorGreenInitialAmplitude;
params.taylorGreenModeX = 1;
params.taylorGreenModeY = 1;
params.taylorGreenThermalNoise = true;
params.taylorGreenForcingEnable = opts.taylorGreenForcingEnable;
params.taylorGreenForcingAmplitude = opts.taylorGreenForcingAmplitude;
params.bodyForceX = 0;
params.bodyForceY = 0;
params.useRandomGridShift = true;
params.projectionEnable = true;
params.projectionStrength = opts.projectionStrength;
params.projectionInterpolationMethod = opts.projectionInterpolationMethod;
params.projectionMomentumCorrectionEnable = true;
params.projectionMomentumCorrectionMode = 'particle_global_exact';
params.thermostatAfterStep = opts.thermostatAfterStep;
params.thermostatAfterProjection = opts.thermostatAfterProjection;
params.thermostatAfterRemap = opts.thermostatAfterRemap;
params.thermostatTargetKBT = opts.thermostatTargetKBT;
params.thermostatStrength = opts.thermostatStrength;
params.thermostatMinParticlesPerCell = 2;
params.thermostatMaxScale = 10;
params.resampParticleMass = opts.particleMass;
params.resampTargetCellMass = opts.gamma * opts.particleMass;
params.resampNTarget = opts.NTarget;
params.resampNMin = opts.NMin;
params.resampNMax = opts.NMax;
params.resampMemoryMinParticles = opts.memoryMinParticles;
params.resampWetMassOnThreshold = opts.wetMassOnThreshold;
params.resampWetMassOffThreshold = opts.wetMassOffThreshold;
params.resampWetParticleOnThreshold = opts.wetParticleOnThreshold;
params.resampWetParticleOffThreshold = opts.wetParticleOffThreshold;
params.resampExtractEvery = opts.extractEvery;
params.resampExtractSelectionMode = opts.extractSelectionMode;
params.resampPreservePreEditVelocity = opts.preservePreEditVelocity;
params.resampMassSafetyEnable = opts.remapMassSafetyEnable;
params.resampMassSafetyMode = opts.remapMassSafetyMode;
params.resampMassSafetyMinFactor = opts.remapMassSafetyMinFactor;
params.resampMassSafetyMaxFactor = opts.remapMassSafetyMaxFactor;
params.computeDiagnostics = true;
params.visualMaxParticles = opts.visualMaxParticles;
params.visualQuiverScale = opts.visualQuiverScale;
end

function opts = parse_options(varargin)
opts = struct();
opts.method = 'weighted_q6';
opts.steps = 300;
opts.summaryEvery = 25;
opts.visualEvery = 5;
opts.Nx = 32;
opts.Ny = 32;
opts.gamma = 20;
opts.dt = 0.001;
opts.alphaDeg = 90;
opts.kBT = 0.01;
opts.particleMass = 1.0;
opts.taylorGreenAmplitude = 0.10;
opts.taylorGreenInitialAmplitude = 0.10;
opts.taylorGreenForcingEnable = false;
opts.taylorGreenForcingAmplitude = 0.0;
opts.projectionStrength = 1.0;
opts.projectionInterpolationMethod = 'nearest';
opts.thermostatAfterStep = false;
opts.thermostatAfterProjection = false;
opts.thermostatAfterRemap = false;
opts.thermostatTargetKBT = [];
opts.thermostatStrength = 1.0;
opts.rngSeed = 12345;
opts.writeCsv = true;
opts.outputDir = fullfile('runs', 'resamp_pool_insertion_visual');
opts.capacityFactor = 2.0;
opts.NTarget = [];
opts.NMin = [];
opts.NMax = [];
opts.memoryMinParticles = [];
opts.initialWetMaskMode = 'all';
opts.wetMaskUpdateMode = 'none';
opts.wetMaskUpdateEvery = 0;
opts.cellWetMask = [];
opts.wetMassOnThreshold = [];
opts.wetMassOffThreshold = [];
opts.wetParticleOnThreshold = [];
opts.wetParticleOffThreshold = [];
opts.extractEvery = 0;
opts.extractSelectionMode = 'closest_to_cell_mean';
opts.preservePreEditVelocity = true;
opts.insertEvery = 1;
opts.insertVelocityMode = 'current_or_memory_pairwise';
opts.insertKBT = [];
opts.remapEvery = 1;
opts.remapMethod = 'scale_preserve_velocity';
opts.targetVelocityMode = 'preserve_cell_velocity';
opts.massMinFactor = 0.05;
opts.massMaxFactor = 20.0;
opts.remapMassSafetyEnable = true;
opts.remapMassSafetyMode = 'uniform_mass_velocity_shift';
opts.remapMassSafetyMinFactor = 0.25;
opts.remapMassSafetyMaxFactor = 4.0;
opts.constraintTolerance = 1e-10;
opts.initialDepletion = 'patch';
opts.depletionPatchSize = [6 6];
opts.depletionPatchCenter = [];
opts.overpopulationPatchSize = [];
opts.overpopulationPatchCenter = [];
opts.populationHeterogeneityStd = 4.0;
opts.populationHeterogeneityMin = [];
opts.populationHeterogeneityMax = [];
opts.populationHeterogeneityVelocityMode = 'taylor_green_at_new_position';
opts.populationHeterogeneityThermalNoise = true;
opts.populationHeterogeneityZeroGlobalMean = true;
opts.figureId = 620;
opts.debugFigureId = 621;
opts.showDebugFigure = false;
opts.saveFrames = false;
opts.frameDir = fullfile('runs', 'resamp_pool_insertion_visual', 'frames');
opts.framePrefix = 'resamp_pool';
opts.visualMaxParticles = 12000;
opts.visualQuiverScale = 1.2;

if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    key = char(string(varargin{k}));
    val = varargin{k+1};
    switch lower(key)
        case 'method'
            opts.method = lower(strrep(char(string(val)), '-', '_'));
        case 'steps'
            opts.steps = val;
        case 'summaryevery'
            opts.summaryEvery = val;
        case 'visualevery'
            opts.visualEvery = val;
        case 'nx'
            opts.Nx = val;
        case 'ny'
            opts.Ny = val;
        case 'gamma'
            opts.gamma = val;
        case 'dt'
            opts.dt = val;
        case 'alphadeg'
            opts.alphaDeg = val;
        case 'kbt'
            opts.kBT = val;
        case {'particlemass','resampparticlemass'}
            opts.particleMass = val;
        case 'taylorgreenamplitude'
            opts.taylorGreenAmplitude = val;
        case 'taylorgreeninitialamplitude'
            opts.taylorGreenInitialAmplitude = val;
        case 'taylorgreenforcingenable'
            opts.taylorGreenForcingEnable = logical(val);
        case 'taylorgreenforcingamplitude'
            opts.taylorGreenForcingAmplitude = val;
        case 'projectionstrength'
            opts.projectionStrength = val;
        case 'projectioninterpolationmethod'
            opts.projectionInterpolationMethod = char(string(val));
        case 'thermostatafterstep'
            opts.thermostatAfterStep = logical(val);
        case 'thermostatafterprojection'
            opts.thermostatAfterProjection = logical(val);
        case 'thermostatafterremap'
            opts.thermostatAfterRemap = logical(val);
        case 'thermostattargetkbt'
            opts.thermostatTargetKBT = val;
        case 'thermostatstrength'
            opts.thermostatStrength = val;
        case 'rngseed'
            opts.rngSeed = val;
        case 'writecsv'
            opts.writeCsv = logical(val);
        case 'outputdir'
            opts.outputDir = char(string(val));
        case 'capacityfactor'
            opts.capacityFactor = val;
        case 'ntarget'
            opts.NTarget = val;
        case 'nmin'
            opts.NMin = val;
        case 'nmax'
            opts.NMax = val;
        case 'memoryminparticles'
            opts.memoryMinParticles = val;
        case {'initialwetmaskmode','cellwetmaskmode'}
            opts.initialWetMaskMode = lower(char(string(val)));
        case 'wetmaskupdatemode'
            opts.wetMaskUpdateMode = lower(char(string(val)));
        case 'wetmaskupdateevery'
            opts.wetMaskUpdateEvery = val;
        case {'cellwetmask','wetmask'}
            opts.cellWetMask = logical(val);
        case {'wetmassonthreshold','massonthreshold'}
            opts.wetMassOnThreshold = val;
        case {'wetmassoffthreshold','massoffthreshold'}
            opts.wetMassOffThreshold = val;
        case {'wetparticleonthreshold','particleonthreshold'}
            opts.wetParticleOnThreshold = val;
        case {'wetparticleoffthreshold','particleoffthreshold'}
            opts.wetParticleOffThreshold = val;
        case 'extractevery'
            opts.extractEvery = val;
        case {'extractselectionmode','selectionmode'}
            opts.extractSelectionMode = lower(char(string(val)));
        case {'preservepreeditvelocity','populationeditpreservepreeditvelocity'}
            opts.preservePreEditVelocity = logical(val);
        case 'insertevery'
            opts.insertEvery = val;
        case 'insertvelocitymode'
            opts.insertVelocityMode = lower(char(string(val)));
        case 'insertkbt'
            opts.insertKBT = val;
        case 'remapevery'
            opts.remapEvery = val;
        case 'remapmethod'
            opts.remapMethod = lower(char(string(val)));
        case 'targetvelocitymode'
            opts.targetVelocityMode = lower(char(string(val)));
        case 'massminfactor'
            opts.massMinFactor = val;
        case 'massmaxfactor'
            opts.massMaxFactor = val;
        case {'remapmasssafetyenable','masssafetyenable'}
            opts.remapMassSafetyEnable = logical(val);
        case {'remapmasssafetymode','masssafetymode'}
            opts.remapMassSafetyMode = lower(char(string(val)));
        case {'remapmasssafetyminfactor','masssafetyminfactor'}
            opts.remapMassSafetyMinFactor = val;
        case {'remapmasssafetymaxfactor','masssafetymaxfactor'}
            opts.remapMassSafetyMaxFactor = val;
        case 'constrainttolerance'
            opts.constraintTolerance = val;
        case 'initialdepletion'
            opts.initialDepletion = lower(char(string(val)));
        case 'depletionpatchsize'
            opts.depletionPatchSize = val;
        case 'depletionpatchcenter'
            opts.depletionPatchCenter = val;
        case {'overpopulationpatchsize','richpatchsize'}
            opts.overpopulationPatchSize = val;
        case {'overpopulationpatchcenter','richpatchcenter'}
            opts.overpopulationPatchCenter = val;
        case {'populationheterogeneitystd','heterogeneitystd','initialpopulationstd'}
            opts.populationHeterogeneityStd = val;
        case {'populationheterogeneitymin','heterogeneitymin','initialpopulationmin'}
            opts.populationHeterogeneityMin = val;
        case {'populationheterogeneitymax','heterogeneitymax','initialpopulationmax'}
            opts.populationHeterogeneityMax = val;
        case {'populationheterogeneityvelocitymode','heterogeneityvelocitymode'}
            opts.populationHeterogeneityVelocityMode = lower(char(string(val)));
        case {'populationheterogeneitythermalnoise','heterogeneitythermalnoise'}
            opts.populationHeterogeneityThermalNoise = logical(val);
        case {'populationheterogeneityzeroglobalmean','heterogeneityzeroglobalmean'}
            opts.populationHeterogeneityZeroGlobalMean = logical(val);
        case 'figureid'
            opts.figureId = val;
        case 'debugfigureid'
            opts.debugFigureId = val;
        case 'showdebugfigure'
            opts.showDebugFigure = logical(val);
        case {'saveframes','saveframe'}
            opts.saveFrames = logical(val);
        case 'framedir'
            opts.frameDir = char(string(val));
        case 'frameprefix'
            opts.framePrefix = char(string(val));
        case 'visualmaxparticles'
            opts.visualMaxParticles = val;
        case 'visualquiverscale'
            opts.visualQuiverScale = val;
        otherwise
            error('Unknown option: %s', key);
    end
end
if isempty(opts.NTarget)
    opts.NTarget = opts.gamma;
end
if isempty(opts.NMin)
    opts.NMin = ceil(0.5 * opts.gamma);
end
if isempty(opts.NMax)
    opts.NMax = ceil(1.5 * opts.gamma);
end
if isempty(opts.memoryMinParticles)
    opts.memoryMinParticles = opts.NMin;
end
if isempty(opts.insertKBT)
    opts.insertKBT = opts.kBT;
end
if isempty(opts.thermostatTargetKBT)
    opts.thermostatTargetKBT = opts.kBT;
end
if isempty(opts.wetMassOnThreshold)
    opts.wetMassOnThreshold = 0.5 * opts.NTarget * opts.particleMass;
end
if isempty(opts.wetMassOffThreshold)
    opts.wetMassOffThreshold = 0.05 * opts.NTarget * opts.particleMass;
end
if isempty(opts.wetParticleOnThreshold)
    opts.wetParticleOnThreshold = max(1, ceil(0.5 * opts.NTarget));
end
if isempty(opts.wetParticleOffThreshold)
    opts.wetParticleOffThreshold = 0;
end
end

function [md, tg] = diagnostics(state, params)
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
if exist('projection_taylor_green_diagnostics', 'file') == 2
    tg = projection_taylor_green_diagnostics(md.G, params);
else
    tg = struct('modeAmplitude', NaN, 'modeCoherence', NaN, 'enstrophy', NaN);
end
end

function state = initialize_velocity_memory(state, params)
activeMask = resamp_active_mask(state);
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', true, 'periodicY', true, 'activeMask', activeMask);
state.uMemUx = G.Ux;
state.uMemUy = G.Uy;
state.uMemValid = G.N > 0;
end

function [state, info] = apply_initial_depletion(state, params, opts)
%APPLY_INITIAL_DEPLETION Deterministic initial population edits for pool tests.
%
% Modes:
%   none/off             leave exact-per-cell initialization unchanged;
%   patch                deactivate all particles in one patch (legacy empty pocket);
%   paired_pockets       move particles from an empty patch into a rich patch;
%   heterogeneous_population distributed random target N_c with conserved Nactive.
%
% The paired-pockets mode is designed to exercise the complete recycling loop:
% extraction from overpopulated cells -> free pool -> insertion into empty cells -> remap.
% It keeps Nactive unchanged at t=0 while creating exactly one poor pocket and one rich pocket.

mode = lower(char(string(opts.initialDepletion)));
info = empty_initial_population_edit_info(mode);
if strcmp(mode, 'none') || strcmp(mode, 'off')
    return;
end
Nx = params.Nx;
Ny = params.Ny;
emptyCells = patch_cell_ids(Nx, Ny, opts.depletionPatchSize, opts.depletionPatchCenter, [floor(Nx/2), floor(Ny/2)]);

switch mode
    case 'patch'
        [state, killedIds] = deactivate_cells(state, params, emptyCells);
        info.mode = mode;
        info.nEmptyPatchCells = numel(emptyCells);
        info.nParticlesRemoved = numel(killedIds);
        info.nParticlesMovedToRichPatch = 0;
        info.NactiveAfter = nnz(resamp_active_mask(state));
        return;

    case {'paired_pockets','empty_overpop','empty_and_overpop','empty_and_overpopulated','empty_overpopulated_patches'}
        if isempty(opts.overpopulationPatchSize)
            richSize = opts.depletionPatchSize;
        else
            richSize = opts.overpopulationPatchSize;
        end
        richDefaultCenter = [max(1, floor(Nx/4)), max(1, floor(Ny/4))];
        richCells = patch_cell_ids(Nx, Ny, richSize, opts.overpopulationPatchCenter, richDefaultCenter);
        if any(ismember(emptyCells, richCells))
            error(['Initial paired pockets overlap. Choose non-overlapping depletionPatchCenter ', ...
                   'and overpopulationPatchCenter.']);
        end
        [state, movedIds, movedPerCell] = move_particles_from_empty_to_rich_patch(state, params, emptyCells, richCells);
        info.mode = mode;
        info.nEmptyPatchCells = numel(emptyCells);
        info.nRichPatchCells = numel(richCells);
        info.nParticlesRemoved = 0;
        info.nParticlesMovedToRichPatch = numel(movedIds);
        info.NactiveAfter = nnz(resamp_active_mask(state));
        info.emptyPatchCells = emptyCells(:).';
        info.richPatchCells = richCells(:).';
        info.movedPerRichCell = movedPerCell(:).';
        return;

    case {'heterogeneous_population','population_heterogeneous','heterogeneous','hetero','random_population','population_noise'}
        popMin = opts.populationHeterogeneityMin;
        popMax = opts.populationHeterogeneityMax;
        if isempty(popMin)
            popMin = max(0, floor(params.gamma - 3 * opts.populationHeterogeneityStd));
        end
        if isempty(popMax)
            popMax = ceil(params.gamma + 3 * opts.populationHeterogeneityStd);
        end
        [state, heteroInfo] = resamp_apply_population_heterogeneity(state, params, ...
            'populationStd', opts.populationHeterogeneityStd, ...
            'populationMin', popMin, ...
            'populationMax', popMax, ...
            'velocityMode', opts.populationHeterogeneityVelocityMode, ...
            'thermalNoise', opts.populationHeterogeneityThermalNoise, ...
            'zeroGlobalMean', opts.populationHeterogeneityZeroGlobalMean, ...
            'computeDiagnostics', true);
        info = merge_initial_info(info, heteroInfo);
        return;

    otherwise
        error('Unknown initialDepletion mode: %s', mode);
end
end

function info = empty_initial_population_edit_info(mode)
info = struct();
info.mode = mode;
info.nEmptyPatchCells = 0;
info.nRichPatchCells = 0;
info.nParticlesRemoved = 0;
info.nParticlesMovedToRichPatch = 0;
info.NactiveAfter = NaN;
info.emptyPatchCells = [];
info.richPatchCells = [];
info.movedPerRichCell = [];
info.populationStdRequested = NaN;
info.populationMinClamp = NaN;
info.populationMaxClamp = NaN;
info.nMovedParticles = NaN;
info.nDonorCells = NaN;
info.nReceiverCells = NaN;
info.targetNMin = NaN;
info.targetNMax = NaN;
info.targetNStd = NaN;
info.NMinAfter = NaN;
info.NMaxAfter = NaN;
info.NStdAfter = NaN;
info.MRelRmsAfter = NaN;
end


function info = merge_initial_info(info, extra)
if ~isstruct(extra)
    return;
end
names = fieldnames(extra);
for kk = 1:numel(names)
    info.(names{kk}) = extra.(names{kk});
end
end

function cells = patch_cell_ids(Nx, Ny, patchSize, patchCenter, defaultCenter)
if isempty(patchSize)
    patchSize = [6 6];
end
sx = patchSize(1);
sy = patchSize(2);
if isempty(patchCenter)
    cx = defaultCenter(1);
    cy = defaultCenter(2);
else
    cx = patchCenter(1);
    cy = patchCenter(2);
end
ixList = mod((cx - floor(sx/2)):(cx - floor(sx/2) + sx - 1) - 1, Nx) + 1;
iyList = mod((cy - floor(sy/2)):(cy - floor(sy/2) + sy - 1) - 1, Ny) + 1;
[IX, IY] = ndgrid(ixList, iyList);
cells = (IX(:) - 1) * Ny + IY(:);
cells = unique(cells(:), 'stable');
end

function [state, idsKill] = deactivate_cells(state, params, cellsToDeactivate)
activeMask = resamp_active_mask(state);
cellId = resamp_cell_ids_periodic(state.x(activeMask,:), params, 'periodicX', true, 'periodicY', true);
activeIds = find(activeMask);
kill = ismember(cellId(:), cellsToDeactivate(:));
idsKill = activeIds(kill);
state.active(idsKill) = false;
state.m(idsKill) = 0;
state.v(idsKill,:) = 0;
state.x(idsKill,:) = 0;
state.Nactive = nnz(resamp_active_mask(state));
end

function [state, movedIds, movedPerRichCell] = move_particles_from_empty_to_rich_patch(state, params, emptyCells, richCells)
% Move active particles from one patch into another one without changing Nactive.
% This creates a true empty pocket and a true overpopulated pocket at t=0.
Nx = params.Nx;
Ny = params.Ny;
dx = params.Lx / Nx;
dy = params.Ly / Ny;
activeMask = resamp_active_mask(state);
cellId = resamp_cell_ids_periodic(state.x(activeMask,:), params, 'periodicX', true, 'periodicY', true);
activeIds = find(activeMask);
moveMask = ismember(cellId(:), emptyCells(:));
movedIds = activeIds(moveMask);
if isempty(movedIds)
    movedPerRichCell = zeros(numel(richCells), 1);
    state.Nactive = nnz(resamp_active_mask(state));
    return;
end
nMove = numel(movedIds);
nRich = numel(richCells);
richAssign = richCells(mod((0:nMove-1)', nRich) + 1);
movedPerRichCell = accumarray(mod((0:nMove-1)', nRich) + 1, 1, [nRich, 1], @sum, 0);
for kk = 1:nMove
    c = richAssign(kk);
    ix = floor((c - 1) / Ny) + 1;
    iy = c - Ny * (ix - 1);
    state.x(movedIds(kk),1) = (ix - 1 + rand()) * dx;
    state.x(movedIds(kk),2) = (iy - 1 + rand()) * dy;
end
state.Nactive = nnz(resamp_active_mask(state));
end

function row = empty_row()
row = struct();
row.step = NaN;
row.t = NaN;
row.NpActive = NaN;
row.Ncapacity = NaN;
row.Nfree = NaN;
row.activeFraction = NaN;
row.nWetCells = NaN;
row.nDryCells = NaN;
row.wetFraction = NaN;
row.NMean = NaN;
row.NStd = NaN;
row.NMin = NaN;
row.NMax = NaN;
row.nEmptyCells = NaN;
row.MMean = NaN;
row.MStd = NaN;
row.MRelRms = NaN;
row.mParticleMean = NaN;
row.mParticleStd = NaN;
row.mParticleRelStd = NaN;
row.mParticleMin = NaN;
row.mParticleMax = NaN;
row.totalMass = NaN;
row.totalMomentumX = NaN;
row.totalMomentumY = NaN;
row.kBTWeighted = NaN;
row.tgAmplitude = NaN;
row.rmsDivParticleAfter = NaN;
row.dvAppliedRms = NaN;
row.weightedMomentumCorrectionResidual = NaN;
row.collisionDeltaPNorm = NaN;
row.extractedParticles = NaN;
row.extractCells = NaN;
row.poorCellsAfterExtract = NaN;
row.emptyCellsAfterExtract = NaN;
row.overCellsAfterExtract = NaN;
row.extractedParticlesCumulative = NaN;
row.extractCellsCumulative = NaN;
row.maxExtractedParticlesPerStep = NaN;
row.lastExtractionStep = NaN;
row.insertedParticles = NaN;
row.insertCells = NaN;
row.poorCellsAfterInsert = NaN;
row.emptyCellsAfterInsert = NaN;
row.overCellsAfterInsert = NaN;
row.capacityHit = NaN;
row.insertedParticlesCumulative = NaN;
row.insertCellsCumulative = NaN;
row.maxInsertedParticlesPerStep = NaN;
row.lastInsertionStep = NaN;
row.thermostatAfterRemapEnabled = NaN;
row.thermostatAfterRemapMeanKBTBefore = NaN;
row.thermostatAfterRemapMeanKBTAfter = NaN;
row.thermostatAfterRemapMeanScale = NaN;
row.thermostatAfterRemapRmsVelocityChange = NaN;
row.remapMassResidualRelRms = NaN;
row.remapMassSafetyCells = NaN;
row.remapMassSafetyParticles = NaN;
row.remapMassSafetyInfeasibleCells = NaN;
row.remapMassSafetyVelocityShiftRms = NaN;
row.remapMassSafetyVelocityShiftMax = NaN;
row.remapMassSafetyCandidateMinFactor = NaN;
row.remapMassSafetyCandidateMaxFactor = NaN;
row.remapMomentumResidualRms = NaN;
row.remapSuccessFractionNonEmpty = NaN;
end

function row = make_row(step, t, md, tg, stepDiag, insertDiag, extractDiag, remapDiag)
row = empty_row();
row.step = step;
row.t = t;
row.NpActive = md.NpActive;
row.Ncapacity = md.Ncapacity;
row.Nfree = md.Nfree;
row.activeFraction = md.activeFraction;
row.nWetCells = get_field(md, 'nWetCells', NaN);
row.nDryCells = get_field(md, 'nDryCells', NaN);
row.wetFraction = get_field(md, 'wetFraction', NaN);
row.NMean = md.NMean;
row.NStd = md.NStd;
row.NMin = md.NMin;
row.NMax = md.NMax;
row.nEmptyCells = md.nEmptyCells;
row.MMean = md.MMean;
row.MStd = md.MStd;
row.MRelRms = md.MRelRms;
row.mParticleMean = md.mParticleMean;
row.mParticleStd = md.mParticleStd;
row.mParticleRelStd = md.mParticleRelStd;
row.mParticleMin = md.mParticleMin;
row.mParticleMax = md.mParticleMax;
row.totalMass = md.totalMass;
row.totalMomentumX = md.totalMomentum(1);
row.totalMomentumY = md.totalMomentum(2);
row.kBTWeighted = md.kBTWeighted;
row.tgAmplitude = get_field(tg, 'modeAmplitude', NaN);
row.rmsDivParticleAfter = get_nested(stepDiag, {'rmsDivParticleAfter'}, NaN);
row.dvAppliedRms = get_nested(stepDiag, {'dvAppliedRms'}, NaN);
row.weightedMomentumCorrectionResidual = get_nested(stepDiag, {'momentumCorrection','residualDeltaPNorm'}, NaN);
row.collisionDeltaPNorm = get_nested(stepDiag, {'classic','collisionDeltaPNorm'}, get_nested(stepDiag, {'collisionDeltaPNorm'}, NaN));
row.extractedParticles = get_field(extractDiag, 'nExtractedParticles', NaN);
row.extractCells = get_field(extractDiag, 'nCellsExtracted', NaN);
row.poorCellsAfterExtract = get_field(extractDiag, 'nPoorCellsAfter', NaN);
row.emptyCellsAfterExtract = get_field(extractDiag, 'nEmptyCellsAfter', NaN);
row.overCellsAfterExtract = get_field(extractDiag, 'nOverCellsAfter', NaN);
row.extractedParticlesCumulative = get_field(extractDiag, 'nExtractedParticlesCumulative', NaN);
row.extractCellsCumulative = get_field(extractDiag, 'nExtractedCellsCumulative', NaN);
row.maxExtractedParticlesPerStep = get_field(extractDiag, 'maxExtractedParticlesPerStep', NaN);
row.lastExtractionStep = get_field(extractDiag, 'lastExtractionStep', NaN);
row.insertedParticles = get_field(insertDiag, 'nInsertedParticles', NaN);
row.insertCells = get_field(insertDiag, 'nCellsInserted', NaN);
row.poorCellsAfterInsert = get_field(insertDiag, 'nPoorCellsAfter', NaN);
row.emptyCellsAfterInsert = get_field(insertDiag, 'nEmptyCellsAfter', NaN);
row.overCellsAfterInsert = get_field(insertDiag, 'nOverCellsAfter', NaN);
row.capacityHit = double(get_field(insertDiag, 'capacityHit', NaN));
row.insertedParticlesCumulative = get_field(insertDiag, 'nInsertedParticlesCumulative', NaN);
row.insertCellsCumulative = get_field(insertDiag, 'nInsertedCellsCumulative', NaN);
row.maxInsertedParticlesPerStep = get_field(insertDiag, 'maxInsertedParticlesPerStep', NaN);
row.lastInsertionStep = get_field(insertDiag, 'lastInsertionStep', NaN);
row.thermostatAfterRemapEnabled = double(get_nested(stepDiag, {'thermostatAfterRemap','enabled'}, false));
row.thermostatAfterRemapMeanKBTBefore = get_nested(stepDiag, {'thermostatAfterRemap','meanKBTBefore'}, NaN);
row.thermostatAfterRemapMeanKBTAfter = get_nested(stepDiag, {'thermostatAfterRemap','meanKBTAfter'}, NaN);
row.thermostatAfterRemapMeanScale = get_nested(stepDiag, {'thermostatAfterRemap','meanScale'}, NaN);
row.thermostatAfterRemapRmsVelocityChange = get_nested(stepDiag, {'thermostatAfterRemap','rmsVelocityChange'}, NaN);
row.remapMassResidualRelRms = get_field(remapDiag, 'massResidualRelRms', NaN);
row.remapMassSafetyCells = get_field(remapDiag, 'nCellsMassSafetyApplied', NaN);
row.remapMassSafetyParticles = get_field(remapDiag, 'nParticlesMassSafetyApplied', NaN);
row.remapMassSafetyInfeasibleCells = get_field(remapDiag, 'nCellsMassSafetyInfeasible', NaN);
row.remapMassSafetyVelocityShiftRms = get_field(remapDiag, 'massSafetyVelocityShiftRms', NaN);
row.remapMassSafetyVelocityShiftMax = get_field(remapDiag, 'massSafetyVelocityShiftMax', NaN);
row.remapMassSafetyCandidateMinFactor = get_field(remapDiag, 'massSafetyCandidateMinFactor', NaN);
row.remapMassSafetyCandidateMaxFactor = get_field(remapDiag, 'massSafetyCandidateMaxFactor', NaN);
row.remapMomentumResidualRms = get_field(remapDiag, 'momentumResidualRms', NaN);
row.remapSuccessFractionNonEmpty = get_field(remapDiag, 'successFractionNonEmpty', NaN);
end

function d = empty_wet_mask_update_diag()
d = struct('kind','update_cell_wet_mask','mode','none','nWetCellsBefore',NaN, ...
    'nWetCellsAfter',NaN,'nDryCellsAfter',NaN,'wetFractionAfter',NaN, ...
    'nCellsBecameWet',0,'nCellsBecameDry',0);
end

function d = empty_step_diag()
d = struct();
end

function d = empty_insert_diag()
d = struct('nInsertedParticles', NaN, 'nCellsInserted', NaN, 'nPoorCellsAfter', NaN, ...
    'nEmptyCellsAfter', NaN, 'nOverCellsAfter', NaN, 'capacityHit', false, 'insertedPerCellGrid', [], ...
    'nInsertedParticlesCumulative', NaN, 'nInsertedCellsCumulative', NaN, ...
    'maxInsertedParticlesPerStep', NaN, 'lastInsertionStep', NaN);
end

function d = empty_extract_diag()
d = struct('nExtractedParticles', NaN, 'nCellsExtracted', NaN, 'nPoorCellsAfter', NaN, ...
    'nEmptyCellsAfter', NaN, 'nOverCellsAfter', NaN, 'extractedPerCellGrid', [], ...
    'nExtractedParticlesCumulative', NaN, 'nExtractedCellsCumulative', NaN, ...
    'maxExtractedParticlesPerStep', NaN, 'lastExtractionStep', NaN);
end

function d = empty_remap_diag()
d = struct('massResidualRelRms', NaN, 'momentumResidualRms', NaN, 'successFractionNonEmpty', NaN, ...
    'nCellsMassSafetyApplied', NaN, 'nParticlesMassSafetyApplied', NaN, ...
    'nCellsMassSafetyInfeasible', NaN, 'massSafetyVelocityShiftRms', NaN, ...
    'massSafetyVelocityShiftMax', NaN, 'massSafetyCandidateMinFactor', NaN, ...
    'massSafetyCandidateMaxFactor', NaN);
end


function stats = empty_insert_stats()
stats = struct();
stats.nInsertedParticlesCumulative = 0;
stats.nInsertedCellsCumulative = 0;
stats.maxInsertedParticlesPerStep = 0;
stats.lastInsertionStep = NaN;
end

function stats = update_insert_stats(stats, insertDiag, step)
nNow = get_field(insertDiag, 'nInsertedParticles', 0);
cNow = get_field(insertDiag, 'nCellsInserted', 0);
if ~isfinite(nNow)
    nNow = 0;
end
if ~isfinite(cNow)
    cNow = 0;
end
stats.nInsertedParticlesCumulative = stats.nInsertedParticlesCumulative + nNow;
stats.nInsertedCellsCumulative = stats.nInsertedCellsCumulative + cNow;
stats.maxInsertedParticlesPerStep = max(stats.maxInsertedParticlesPerStep, nNow);
if nNow > 0
    stats.lastInsertionStep = step;
end
end

function d = attach_insert_cumulative(d, stats)
d.nInsertedParticlesCumulative = stats.nInsertedParticlesCumulative;
d.nInsertedCellsCumulative = stats.nInsertedCellsCumulative;
d.maxInsertedParticlesPerStep = stats.maxInsertedParticlesPerStep;
d.lastInsertionStep = stats.lastInsertionStep;
end

function stats = empty_extract_stats()
stats = struct();
stats.nExtractedParticlesCumulative = 0;
stats.nExtractedCellsCumulative = 0;
stats.maxExtractedParticlesPerStep = 0;
stats.lastExtractionStep = NaN;
end

function stats = update_extract_stats(stats, extractDiag, step)
nNow = get_field(extractDiag, 'nExtractedParticles', 0);
cNow = get_field(extractDiag, 'nCellsExtracted', 0);
if ~isfinite(nNow)
    nNow = 0;
end
if ~isfinite(cNow)
    cNow = 0;
end
stats.nExtractedParticlesCumulative = stats.nExtractedParticlesCumulative + nNow;
stats.nExtractedCellsCumulative = stats.nExtractedCellsCumulative + cNow;
stats.maxExtractedParticlesPerStep = max(stats.maxExtractedParticlesPerStep, nNow);
if nNow > 0
    stats.lastExtractionStep = step;
end
end

function d = attach_extract_cumulative(d, stats)
d.nExtractedParticlesCumulative = stats.nExtractedParticlesCumulative;
d.nExtractedCellsCumulative = stats.nExtractedCellsCumulative;
d.maxExtractedParticlesPerStep = stats.maxExtractedParticlesPerStep;
d.lastExtractionStep = stats.lastExtractionStep;
end

function info = empty_thermostat_after_remap_diag()
info = struct('enabled', false, 'targetKBT', NaN, 'strength', NaN, ...
    'minParticlesPerCell', NaN, 'maxScale', NaN, 'nThermostattedCells', 0, ...
    'meanScale', NaN, 'minScale', NaN, 'maxScaleApplied', NaN, ...
    'meanKBTBefore', NaN, 'meanKBTAfter', NaN, 'rmsVelocityChange', 0.0);
end

function d = attach_thermostat_after_remap(d, thermostatInfo)
d.thermostatAfterRemap = thermostatInfo;
end

function s = format_scalar(v)
if isempty(v) || ~isfinite(v)
    s = 'n/a';
elseif abs(v - round(v)) < 10*eps(max(1, abs(v)))
    s = sprintf('%d', round(v));
else
    s = sprintf('%.4g', v);
end
end

function v = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end

function v = get_nested(s, names, defaultValue)
v = s;
for k = 1:numel(names)
    if ~isstruct(v) || ~isfield(v, names{k}) || isempty(v.(names{k}))
        v = defaultValue;
        return;
    end
    v = v.(names{k});
end
end
