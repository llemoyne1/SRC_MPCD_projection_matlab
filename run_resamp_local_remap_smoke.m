function out = run_resamp_local_remap_smoke(varargin)
%RUN_RESAMP_LOCAL_REMAP_SMOKE Smoke run for local mass/moment remapping.
%
% This script exercises the first conservative remap layer:
% variable particle masses are adjusted locally, but particles are not yet
% created, deleted, split or merged.
%
% Examples:
%   out = run_resamp_local_remap_smoke();
%   out = run_resamp_local_remap_smoke('method','weighted_q6','steps',500);
%   out = run_resamp_local_remap_smoke('remapMethod','min_change');

opts = parse_options(varargin{:});
if ~isempty(opts.rngSeed)
    rng(opts.rngSeed, 'twister');
end

params = default_params(opts);
[state, initInfo] = resamp_initialize_particles_taylor_green_forced(params);
state = apply_initial_mass_perturbation(state, params, opts);

nRows = floor(opts.steps / opts.summaryEvery) + 2;
rows = repmat(empty_row(), nRows, 1);
irow = 0;

md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
G = md.G;
tg = projection_taylor_green_diagnostics(G, params);
irow = irow + 1;
rows(irow) = make_row(0, 0.0, md, tg, empty_step_diag(), empty_remap_diag());

for step = 1:opts.steps
    switch opts.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_periodic_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag] = resamp_step_projection_periodic_weighted(state, params);
        otherwise
            error('Unknown method: %s', opts.method);
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
        rows(irow) = make_row(step, step * params.dt, md, tg, stepDiag, remapDiag);
        fprintf(['resamp remap %-16s step=%6d t=%.4g Nstd=%.4g Mrel=%.3e ', ...
                 'mRelStd=%.3e remapMassRel=%.3e remapMom=%.3e solved=%.3f kBT=%.5g\n'], ...
            opts.method, step, step*params.dt, md.NStd, md.MRelRms, md.mParticleRelStd, ...
            get_field(remapDiag, 'massResidualRelRms', NaN), ...
            get_field(remapDiag, 'momentumResidualRms', NaN), ...
            get_field(remapDiag, 'successFractionNonEmpty', NaN), md.kBTWeighted);
    end
end
rows = rows(1:irow);
summary = struct2table(rows);

out = struct();
out.params = params;
out.options = opts;
out.initialInfo = initInfo;
out.state = state;
out.summary = summary;
out.finalMassDiagnostics = md;
out.finalTaylorGreenDiagnostics = tg;

if opts.writeCsv
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end
    csvPath = fullfile(opts.outputDir, sprintf('resamp_local_remap_smoke_%s_%s.csv', opts.method, opts.remapMethod));
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
params.computeDiagnostics = opts.computeDiagnostics;
end

function opts = parse_options(varargin)
opts = struct();
opts.method = 'weighted_q6';
opts.steps = 500;
opts.summaryEvery = 50;
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
opts.outputDir = fullfile('runs', 'resamp_local_remap_smoke');
opts.remapEvery = 1;
opts.remapMethod = 'scale_preserve_velocity';
opts.targetVelocityMode = 'preserve_cell_velocity';
opts.massMinFactor = 0.05;
opts.massMaxFactor = 20.0;
opts.constraintTolerance = 1e-10;
opts.initialMassPerturbation = 'cell_log_random';
opts.initialMassPerturbationAmplitude = 0.20;

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
        case {'constrainttolerance','tol'}
            opts.constraintTolerance = val;
        case 'initialmassperturbation'
            opts.initialMassPerturbation = lower(char(string(val)));
        case 'initialmassperturbationamplitude'
            opts.initialMassPerturbationAmplitude = val;
        otherwise
            error('Unknown option: %s', key);
    end
end
end

function state = apply_initial_mass_perturbation(state, params, opts)
mode = lower(opts.initialMassPerturbation);
amp = opts.initialMassPerturbationAmplitude;
if strcmp(mode, 'none') || amp == 0
    return;
end
cellId = resamp_cell_ids_periodic(state.x, params, 'periodicX', true, 'periodicY', true);
Nc = params.Nx * params.Ny;
switch mode
    case {'cell_log_random', 'cell'}
        % One smooth-ish multiplicative factor per cell.  This creates a
        % controlled mass imbalance while leaving population unchanged.
        logScale = amp * randn(Nc, 1);
        scale = exp(logScale);
        scale = scale / mean(scale, 'omitnan');
        state.m = state.m(:) .* scale(cellId);
    case {'particle_log_random', 'particle'}
        scale = exp(amp * randn(numel(state.m), 1));
        scale = scale / mean(scale, 'omitnan');
        state.m = state.m(:) .* scale;
    otherwise
        error('Unknown initialMassPerturbation: %s', opts.initialMassPerturbation);
end
end

function row = empty_row()
row = struct();
row.step = NaN;
row.t = NaN;
row.NMean = NaN;
row.NStd = NaN;
row.NMin = NaN;
row.NMax = NaN;
row.nEmptyCells = NaN;
row.MMean = NaN;
row.MStd = NaN;
row.MMin = NaN;
row.MMax = NaN;
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
row.kineticEnergyMeanWeighted = NaN;
row.tgAmplitude = NaN;
row.tgCoherence = NaN;
row.rmsDivParticleAfter = NaN;
row.dvAppliedRms = NaN;
row.weightedMomentumCorrectionResidual = NaN;
row.remapMassResidualRelRms = NaN;
row.remapMassResidualMaxAbs = NaN;
row.remapMomentumResidualRms = NaN;
row.remapMomentumResidualMaxNorm = NaN;
row.remapDeltaMassRelRms = NaN;
row.remapSuccessFractionNonEmpty = NaN;
row.remapCellsEmpty = NaN;
row.remapCellsUnresolved = NaN;
row.remapCellsUsedSolver = NaN;
row.remapCellsBounded = NaN;
end

function row = make_row(step, t, md, tg, stepDiag, remapDiag)
row = empty_row();
row.step = step;
row.t = t;
row.NMean = md.NMean;
row.NStd = md.NStd;
row.NMin = md.NMin;
row.NMax = md.NMax;
row.nEmptyCells = md.nEmptyCells;
row.MMean = md.MMean;
row.MStd = md.MStd;
row.MMin = md.MMin;
row.MMax = md.MMax;
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
row.kineticEnergyMeanWeighted = md.kineticEnergyMeanWeighted;
row.tgAmplitude = tg.modeAmplitude;
row.tgCoherence = tg.modeCoherence;
row.rmsDivParticleAfter = get_nested(stepDiag, {'rmsDivParticleAfter'}, NaN);
row.dvAppliedRms = get_nested(stepDiag, {'dvAppliedRms'}, NaN);
row.weightedMomentumCorrectionResidual = get_nested(stepDiag, {'momentumCorrection','residualDeltaPNorm'}, NaN);
row.remapMassResidualRelRms = get_field(remapDiag, 'massResidualRelRms', NaN);
row.remapMassResidualMaxAbs = get_field(remapDiag, 'massResidualMaxAbs', NaN);
row.remapMomentumResidualRms = get_field(remapDiag, 'momentumResidualRms', NaN);
row.remapMomentumResidualMaxNorm = get_field(remapDiag, 'momentumResidualMaxNorm', NaN);
row.remapDeltaMassRelRms = get_field(remapDiag, 'deltaMassRelRms', NaN);
row.remapSuccessFractionNonEmpty = get_field(remapDiag, 'successFractionNonEmpty', NaN);
row.remapCellsEmpty = get_field(remapDiag, 'nCellsEmpty', NaN);
row.remapCellsUnresolved = get_field(remapDiag, 'nCellsUnresolved', NaN);
row.remapCellsUsedSolver = get_field(remapDiag, 'nCellsUsedSolver', NaN);
row.remapCellsBounded = get_field(remapDiag, 'nCellsBounded', NaN);
end

function diag = empty_step_diag()
diag = struct();
end

function diag = empty_remap_diag()
diag = struct();
diag.massResidualRelRms = NaN;
diag.massResidualMaxAbs = NaN;
diag.momentumResidualRms = NaN;
diag.momentumResidualMaxNorm = NaN;
diag.deltaMassRelRms = NaN;
diag.successFractionNonEmpty = NaN;
diag.nCellsEmpty = NaN;
diag.nCellsUnresolved = NaN;
diag.nCellsUsedSolver = NaN;
diag.nCellsBounded = NaN;
end

function val = get_nested(s, path, defaultValue)
val = defaultValue;
if ~isstruct(s)
    return;
end
cur = s;
for k = 1:numel(path)
    name = path{k};
    if isstruct(cur) && isfield(cur, name)
        cur = cur.(name);
    else
        return;
    end
end
val = cur;
end

function val = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name)
    val = s.(name);
else
    val = defaultValue;
end
end
