function out = run_resamp_pool_insertion_smoke(varargin)
%RUN_RESAMP_POOL_INSERTION_SMOKE Smoke run for preallocated pool + poor-cell insertion.
%
% This is the third resampling prototype layer:
%   1) convert the weighted state to a preallocated pool;
%   2) optionally create an initial poor-cell patch by deactivating particles;
%   3) after each SRC/Q6 step, insert particles into cells with N < NMin;
%   4) apply the existing local mass/moment remap.
%
% No overpopulation removal or fusion is performed in this patch.
%
% Examples:
%   out = run_resamp_pool_insertion_smoke('method','weighted_classic');
%   out = run_resamp_pool_insertion_smoke('method','weighted_q6','initialDepletion','patch');

opts = parse_options(varargin{:});
if ~isempty(opts.rngSeed)
    rng(opts.rngSeed, 'twister');
end

params = default_params(opts);
[state, initInfo] = resamp_initialize_particles_taylor_green_forced(params);
[state, poolInfo0] = resamp_enable_particle_pool(state, 'capacityFactor', opts.capacityFactor);
state = initialize_velocity_memory(state, params);
state = apply_initial_depletion(state, params, opts);

nRows = floor(opts.steps / opts.summaryEvery) + 2;
rows = repmat(empty_row(), nRows, 1);
irow = 0;

md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
G = md.G;
tg = projection_taylor_green_diagnostics(G, params);
irow = irow + 1;
rows(irow) = make_row(0, 0.0, md, tg, empty_step_diag(), empty_insert_diag(), empty_remap_diag());

for step = 1:opts.steps
    switch opts.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_periodic_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag] = resamp_step_projection_periodic_weighted(state, params);
        otherwise
            error('Unknown method: %s', opts.method);
    end

    if opts.insertEvery > 0 && mod(step, opts.insertEvery) == 0
        [state, insertDiag] = resamp_insert_underpopulated_particles(state, params, ...
            'NTarget', opts.NTarget, ...
            'NMin', opts.NMin, ...
            'NMax', opts.NMax, ...
            'memoryMinParticles', opts.memoryMinParticles, ...
            'particleMass', opts.particleMass, ...
            'kBT', opts.insertKBT, ...
            'insertVelocityMode', opts.insertVelocityMode, ...
            'computeDiagnostics', opts.computeDiagnostics);
    else
        insertDiag = empty_insert_diag();
    end

    if opts.remapEvery > 0 && mod(step, opts.remapEvery) == 0
        [state, remapDiag] = resamp_local_mass_moment_remap(state, params, ...
            'targetCellMass', params.resampTargetCellMass, ...
            'targetVelocityMode', opts.targetVelocityMode, ...
            'method', opts.remapMethod, ...
            'massMin', opts.massMinFactor * params.resampParticleMass, ...
            'massMax', opts.massMaxFactor * params.resampParticleMass, ...
            'constraintTolerance', opts.constraintTolerance, ...
            'computeDiagnostics', opts.computeDiagnostics);
    else
        remapDiag = empty_remap_diag();
    end

    if mod(step, opts.summaryEvery) == 0 || step == opts.steps
        md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
        G = md.G;
        tg = projection_taylor_green_diagnostics(G, params);
        irow = irow + 1;
        rows(irow) = make_row(step, step * params.dt, md, tg, stepDiag, insertDiag, remapDiag);
        fprintf(['resamp pool %-12s step=%6d t=%.4g Nact=%7d free=%7d ', ...
                 'Nmin=%3g Nmax=%3g Nstd=%.4g Mrel=%.3e inserted=%5d poor=%4d capHit=%d ', ...
                 'mRelStd=%.3e remapMassRel=%.3e kBT=%.5g\n'], ...
            opts.method, step, step*params.dt, md.NpActive, md.Nfree, md.NMin, md.NMax, md.NStd, md.MRelRms, ...
            get_field(insertDiag, 'nInsertedParticles', 0), get_field(insertDiag, 'nPoorCellsAfter', NaN), ...
            double(get_field(insertDiag, 'capacityHit', false)), md.mParticleRelStd, ...
            get_field(remapDiag, 'massResidualRelRms', NaN), md.kBTWeighted);
    end
end
rows = rows(1:irow);
summary = struct2table(rows);

out = struct();
out.params = params;
out.options = opts;
out.initialInfo = initInfo;
out.initialPoolInfo = poolInfo0;
out.state = state;
out.summary = summary;
out.finalMassDiagnostics = md;
out.finalTaylorGreenDiagnostics = tg;
out.finalPoolInfo = resamp_particle_pool_info(state);

if opts.writeCsv
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end
    csvPath = fullfile(opts.outputDir, sprintf('resamp_pool_insertion_smoke_%s_%s.csv', opts.method, opts.initialDepletion));
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
params.thermostatTargetKBT = opts.kBT;
params.thermostatStrength = opts.thermostatStrength;
params.thermostatMinParticlesPerCell = 2;
params.thermostatMaxScale = 10;
params.resampParticleMass = opts.particleMass;
params.resampTargetCellMass = opts.gamma * opts.particleMass;
params.resampNTarget = opts.NTarget;
params.resampNMin = opts.NMin;
params.resampNMax = opts.NMax;
params.resampMemoryMinParticles = opts.memoryMinParticles;
params.computeDiagnostics = opts.computeDiagnostics;
end

function opts = parse_options(varargin)
opts = struct();
opts.method = 'weighted_q6';
opts.steps = 200;
opts.summaryEvery = 20;
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
opts.thermostatStrength = 1.0;
opts.computeDiagnostics = true;
opts.rngSeed = 12345;
opts.writeCsv = true;
opts.outputDir = fullfile('runs', 'resamp_pool_insertion_smoke');
opts.capacityFactor = 2.0;
opts.NTarget = [];
opts.NMin = [];
opts.NMax = [];
opts.memoryMinParticles = [];
opts.insertEvery = 1;
opts.insertVelocityMode = 'current_or_memory_pairwise';
opts.insertKBT = [];
opts.remapEvery = 1;
opts.remapMethod = 'scale_preserve_velocity';
opts.targetVelocityMode = 'preserve_cell_velocity';
opts.massMinFactor = 0.05;
opts.massMaxFactor = 20.0;
opts.constraintTolerance = 1e-10;
opts.initialDepletion = 'patch';
opts.depletionPatchSize = [6 6];
opts.depletionPatchCenter = [];

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
        case 'thermostatstrength'
            opts.thermostatStrength = val;
        case 'computediagnostics'
            opts.computeDiagnostics = logical(val);
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
        case 'constrainttolerance'
            opts.constraintTolerance = val;
        case 'initialdepletion'
            opts.initialDepletion = lower(char(string(val)));
        case 'depletionpatchsize'
            opts.depletionPatchSize = val;
        case 'depletionpatchcenter'
            opts.depletionPatchCenter = val;
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
end

function state = initialize_velocity_memory(state, params)
activeMask = resamp_active_mask(state);
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', true, 'periodicY', true, 'activeMask', activeMask);
state.uMemUx = G.Ux;
state.uMemUy = G.Uy;
state.uMemValid = G.N > 0;
end

function state = apply_initial_depletion(state, params, opts)
mode = lower(char(string(opts.initialDepletion)));
if strcmp(mode, 'none') || strcmp(mode, 'off')
    return;
end
if ~strcmp(mode, 'patch')
    error('Unknown initialDepletion mode: %s', mode);
end
Nx = params.Nx;
Ny = params.Ny;
if isempty(opts.depletionPatchCenter)
    cx = floor(Nx/2);
    cy = floor(Ny/2);
else
    cx = opts.depletionPatchCenter(1);
    cy = opts.depletionPatchCenter(2);
end
sx = opts.depletionPatchSize(1);
sy = opts.depletionPatchSize(2);
ixList = mod((cx - floor(sx/2)):(cx - floor(sx/2) + sx - 1) - 1, Nx) + 1;
iyList = mod((cy - floor(sy/2)):(cy - floor(sy/2) + sy - 1) - 1, Ny) + 1;
activeMask = resamp_active_mask(state);
cellId = resamp_cell_ids_periodic(state.x(activeMask,:), params, 'periodicX', true, 'periodicY', true);
activeIds = find(activeMask);
kill = false(size(activeIds));
for k = 1:numel(activeIds)
    c = cellId(k);
    ix = floor((c - 1) / Ny) + 1;
    iy = c - Ny * (ix - 1);
    kill(k) = any(ixList == ix) && any(iyList == iy);
end
idsKill = activeIds(kill);
state.active(idsKill) = false;
state.m(idsKill) = 0;
state.v(idsKill,:) = 0;
state.x(idsKill,:) = 0;
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
row.insertedParticles = NaN;
row.insertCells = NaN;
row.poorCellsAfterInsert = NaN;
row.emptyCellsAfterInsert = NaN;
row.overCellsAfterInsert = NaN;
row.capacityHit = NaN;
row.remapMassResidualRelRms = NaN;
row.remapMomentumResidualRms = NaN;
row.remapSuccessFractionNonEmpty = NaN;
end

function row = make_row(step, t, md, tg, stepDiag, insertDiag, remapDiag)
row = empty_row();
row.step = step;
row.t = t;
row.NpActive = md.NpActive;
row.Ncapacity = md.Ncapacity;
row.Nfree = md.Nfree;
row.activeFraction = md.activeFraction;
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
row.tgAmplitude = tg.modeAmplitude;
row.rmsDivParticleAfter = get_nested(stepDiag, {'rmsDivParticleAfter'}, NaN);
row.dvAppliedRms = get_nested(stepDiag, {'dvAppliedRms'}, NaN);
row.weightedMomentumCorrectionResidual = get_nested(stepDiag, {'momentumCorrection','residualDeltaPNorm'}, NaN);
row.collisionDeltaPNorm = get_nested(stepDiag, {'classic','collisionDeltaPNorm'}, get_nested(stepDiag, {'collisionDeltaPNorm'}, NaN));
row.insertedParticles = get_field(insertDiag, 'nInsertedParticles', NaN);
row.insertCells = get_field(insertDiag, 'nCellsInserted', NaN);
row.poorCellsAfterInsert = get_field(insertDiag, 'nPoorCellsAfter', NaN);
row.emptyCellsAfterInsert = get_field(insertDiag, 'nEmptyCellsAfter', NaN);
row.overCellsAfterInsert = get_field(insertDiag, 'nOverCellsAfter', NaN);
row.capacityHit = double(get_field(insertDiag, 'capacityHit', NaN));
row.remapMassResidualRelRms = get_field(remapDiag, 'massResidualRelRms', NaN);
row.remapMomentumResidualRms = get_field(remapDiag, 'momentumResidualRms', NaN);
row.remapSuccessFractionNonEmpty = get_field(remapDiag, 'successFractionNonEmpty', NaN);
end

function d = empty_step_diag()
d = struct();
end

function d = empty_insert_diag()
d = struct('nInsertedParticles', NaN, 'nCellsInserted', NaN, 'nPoorCellsAfter', NaN, ...
    'nEmptyCellsAfter', NaN, 'nOverCellsAfter', NaN, 'capacityHit', false);
end

function d = empty_remap_diag()
d = struct('massResidualRelRms', NaN, 'momentumResidualRms', NaN, 'successFractionNonEmpty', NaN);
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
