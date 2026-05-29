function out = run_resamp_periodic_step_visual_demo(varargin)
%RUN_RESAMP_PERIODIC_STEP_VISUAL_DEMO Periodic step/shear visual demo for weighted resampling.
%
% This run is the first dynamic stress case after the depleted-patch smoke.
% It is fully periodic: no wall, no solid mask, no inlet/outlet.  A staircase
% velocity field creates strong shear and local population transport while the
% weighted pool/remap machinery maintains mass and support.
%
% Example:
%   out = run_resamp_periodic_step_visual_demo('method','weighted_q6');

opts = parse_options(varargin{:});
if ~isempty(opts.rngSeed)
    rng(opts.rngSeed, 'twister');
end

params = default_params(opts);
params.method = opts.method;
params.visualFigureId = opts.figureId;
params.visualMaxParticles = opts.visualMaxParticles;
params.visualQuiverScale = opts.visualQuiverScale;

[state, initInfo] = resamp_initialize_particles_periodic_step(params);
[state, poolInfo0] = resamp_enable_particle_pool(state, 'capacityFactor', opts.capacityFactor);
state = initialize_velocity_memory(state, params);

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
thermostatAfterRemapDiag = empty_thermostat_after_remap_diag();
insertDiag0 = attach_insert_cumulative(empty_insert_diag(), insertStats);
stepDiag0 = attach_thermostat_after_remap(empty_step_diag(), thermostatAfterRemapDiag);
[md, stepMetrics] = diagnostics(state, params);
rows(irow) = make_row(0, 0.0, md, stepMetrics, stepDiag0, insertDiag0, empty_remap_diag());

resamp_periodic_step_visualize_frame(state, params, 0, ...
    'insertDiag', insertDiag0, ...
    'remapDiag', empty_remap_diag(), ...
    'stepDiag', stepDiag0, ...
    'figureId', opts.figureId, ...
    'saveFrame', opts.saveFrames, ...
    'frameDir', opts.frameDir, ...
    'framePrefix', opts.framePrefix);

lastStepDiag = stepDiag0;
lastInsertDiag = insertDiag0;
lastRemapDiag = empty_remap_diag();
lastThermostatAfterRemapDiag = thermostatAfterRemapDiag;

for step = 1:opts.steps
    switch opts.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_periodic_weighted(state, params);
            [state, forceDiag] = apply_periodic_step_forcing(state, params, opts);
            stepDiag.kolmogorovForce = forceDiag;
        case 'weighted_q6'
            % Split the Q6 step here so that the structured body force is applied
            % after SRC streaming/collision and before the incompressible projection.
            [stateClassic, classicDiag] = resamp_step_classic_periodic_weighted(state, params);
            [stateForced, forceDiag] = apply_periodic_step_forcing(stateClassic, params, opts);
            [state, stepDiag] = resamp_apply_q6_periodic_weighted(stateForced, params);
            stepDiag.classic = classicDiag;
            stepDiag.kolmogorovForce = forceDiag;
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
            'computeDiagnostics', true);
    else
        insertDiag = empty_insert_diag();
    end
    insertStats = update_insert_stats(insertStats, insertDiag, step);
    insertDiag = attach_insert_cumulative(insertDiag, insertStats);
    lastInsertDiag = insertDiag;

    if opts.remapEvery > 0 && mod(step, opts.remapEvery) == 0
        [state, remapDiag] = resamp_local_mass_moment_remap(state, params, ...
            'targetCellMass', params.resampTargetCellMass, ...
            'targetVelocityMode', opts.targetVelocityMode, ...
            'method', opts.remapMethod, ...
            'massMin', opts.massMinFactor * params.resampParticleMass, ...
            'massMax', opts.massMaxFactor * params.resampParticleMass, ...
            'constraintTolerance', opts.constraintTolerance, ...
            'computeDiagnostics', true);
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
        [md, stepMetrics] = diagnostics(state, params);
    end
    if doSummary
        irow = irow + 1;
        rows(irow) = make_row(step, step * params.dt, md, stepMetrics, stepDiag, insertDiag, remapDiag);
        fprintf(['periodic-step %-12s step=%6d t=%.4g Nact=%7d free=%7d ', ...
                 'N[min,max]=[%3g,%3g] Mrel=%.3e insertedNow=%s insertedCum=%s lastInsert=%s ', ...
                 'poor=%4g over=%4g mRelStd=%.3e kBT=%.5g thermAfter=%.5g omegaRms=%.4g\n'], ...
            opts.method, step, step*params.dt, md.NpActive, md.Nfree, md.NMin, md.NMax, md.MRelRms, ...
            format_scalar(get_field(insertDiag, 'nInsertedParticles', NaN)), ...
            format_scalar(get_field(insertDiag, 'nInsertedParticlesCumulative', NaN)), ...
            format_scalar(get_field(insertDiag, 'lastInsertionStep', NaN)), ...
            get_field(insertDiag, 'nPoorCellsAfter', NaN), ...
            get_field(insertDiag, 'nOverCellsAfter', NaN), md.mParticleRelStd, md.kBTWeighted, ...
            get_nested(stepDiag, {'thermostatAfterRemap','meanKBTAfter'}, NaN), stepMetrics.omegaRms);
    end
    if doVisual
        resamp_periodic_step_visualize_frame(state, params, step, ...
            'insertDiag', insertDiag, ...
            'remapDiag', remapDiag, ...
            'stepDiag', stepDiag, ...
            'figureId', opts.figureId, ...
            'saveFrame', opts.saveFrames, ...
            'frameDir', opts.frameDir, ...
            'framePrefix', opts.framePrefix);
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
out.finalStepMetrics = stepMetrics;
out.finalPoolInfo = resamp_particle_pool_info(state);
out.lastStepDiag = lastStepDiag;
out.lastInsertDiag = lastInsertDiag;
out.lastRemapDiag = lastRemapDiag;
out.lastThermostatAfterRemapDiag = lastThermostatAfterRemapDiag;
out.insertStats = insertStats;

if opts.writeCsv
    csvPath = fullfile(opts.outputDir, sprintf('resamp_periodic_step_visual_%s.csv', opts.method));
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
params.bodyForceX = 0;
params.bodyForceY = 0;
params.forceMode = opts.forceMode;
params.kolmogorovForceAmplitude = opts.kolmogorovForceAmplitude;
params.kolmogorovForceWaveNumber = opts.kolmogorovForceWaveNumber;
params.kolmogorovForcePhase = opts.kolmogorovForcePhase;
params.kolmogorovForceDirection = opts.kolmogorovForceDirection;
params.kolmogorovForceZeroMeanKick = opts.kolmogorovForceZeroMeanKick;
params.kolmogorovForceApplyDt = opts.kolmogorovForceApplyDt;
params.resampParticleMass = opts.particleMass;
params.resampTargetCellMass = opts.gamma * opts.particleMass;
params.resampNTarget = opts.NTarget;
params.resampNMin = opts.NMin;
params.resampNMax = opts.NMax;
params.resampMemoryMinParticles = opts.memoryMinParticles;
params.computeDiagnostics = true;
params.visualMaxParticles = opts.visualMaxParticles;
params.visualQuiverScale = opts.visualQuiverScale;
params.periodicStepULower = opts.ULower;
params.periodicStepUUpper = opts.UUpper;
params.periodicStepVLower = opts.VLower;
params.periodicStepVUpper = opts.VUpper;
params.periodicStepXFraction = opts.stepXFraction;
params.periodicStepYLowFraction = opts.stepYLowFraction;
params.periodicStepYHighFraction = opts.stepYHighFraction;
params.periodicStepTransitionWidth = opts.stepTransitionWidth;
params.periodicStepThermalNoise = opts.stepThermalNoise;
params.periodicStepSubtractMeanVelocity = opts.stepSubtractMeanVelocity;
end

function opts = parse_options(varargin)
opts = struct();
opts.method = 'weighted_q6';
opts.steps = 1000;
opts.summaryEvery = 25;
opts.visualEvery = 10;
opts.Nx = 64;
opts.Ny = 32;
opts.gamma = 20;
opts.dt = 0.01;
opts.alphaDeg = 90;
opts.kBT = 0.01;
opts.particleMass = 1.0;
opts.ULower = 0.20;
opts.UUpper = -0.05;
opts.VLower = 0.00;
opts.VUpper = 0.00;
opts.stepXFraction = 0.35;
opts.stepYLowFraction = 0.35;
opts.stepYHighFraction = 0.65;
opts.stepTransitionWidth = 0.0;
opts.stepThermalNoise = true;
opts.stepSubtractMeanVelocity = true;
opts.projectionStrength = 1.0;
opts.projectionInterpolationMethod = 'nearest';
opts.thermostatAfterStep = false;
opts.thermostatAfterProjection = false;
opts.thermostatAfterRemap = true;
opts.thermostatTargetKBT = [];
opts.thermostatStrength = 0.25;
opts.forceMode = 'none';
opts.kolmogorovForceAmplitude = 0.0;
opts.kolmogorovForceWaveNumber = 1;
opts.kolmogorovForcePhase = 0.0;
opts.kolmogorovForceDirection = 'x';
opts.kolmogorovForceZeroMeanKick = true;
opts.kolmogorovForceApplyDt = true;
opts.rngSeed = 12345;
opts.writeCsv = true;
opts.outputDir = fullfile('runs', 'resamp_periodic_step_visual');
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
opts.figureId = 640;
opts.saveFrames = false;
opts.frameDir = fullfile('runs', 'resamp_periodic_step_visual', 'frames');
opts.framePrefix = 'resamp_periodic_step';
opts.visualMaxParticles = 15000;
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
        case {'ulower','periodicstepulower'}
            opts.ULower = val;
        case {'uupper','periodicstepuupper'}
            opts.UUpper = val;
        case {'vlower','periodicstepvlower'}
            opts.VLower = val;
        case {'vupper','periodicstepvupper'}
            opts.VUpper = val;
        case {'stepxfraction','periodicstepxfraction'}
            opts.stepXFraction = val;
        case {'stepylowfraction','periodicstepylowfraction'}
            opts.stepYLowFraction = val;
        case {'stepyhighfraction','periodicstepyhighfraction'}
            opts.stepYHighFraction = val;
        case {'steptransitionwidth','periodicsteptransitionwidth'}
            opts.stepTransitionWidth = val;
        case {'stepthermalnoise','periodicstepthermalnoise'}
            opts.stepThermalNoise = logical(val);
        case {'stepsubtractmeanvelocity','periodicstepsubtractmeanvelocity'}
            opts.stepSubtractMeanVelocity = logical(val);
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
        case {'forcemode','bodyforcemode'}
            opts.forceMode = lower(strrep(char(string(val)), '-', '_'));
        case {'kolmogorovforceamplitude','forceamplitude','bodyforceamplitude','kolmogorovamplitude'}
            opts.kolmogorovForceAmplitude = val;
        case {'kolmogorovforcewavenumber','kolmogorovwavenumber','forcemodenumber','kolmogorovk'}
            opts.kolmogorovForceWaveNumber = val;
        case {'kolmogorovforcephase','kolmogorovphase','forcephase'}
            opts.kolmogorovForcePhase = val;
        case {'kolmogorovforcedirection','forcedirection'}
            opts.kolmogorovForceDirection = lower(char(string(val)));
        case {'kolmogorovforcezeromeankick','forcemeanmomentumcorrection','zeromeanforcekick'}
            opts.kolmogorovForceZeroMeanKick = logical(val);
        case {'kolmogorovforceapplydt','forceapplydt'}
            opts.kolmogorovForceApplyDt = logical(val);
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
        case 'figureid'
            opts.figureId = val;
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
end

function [stateOut, forceDiag] = apply_periodic_step_forcing(state, params, opts)
mode = lower(strrep(char(string(opts.forceMode)), '-', '_'));
stateOut = state;
switch mode
    case {'none','off','disabled','no',''}
        forceDiag = empty_force_diag(mode);
    case {'kolmogorov','kolmogorov_x','kolmogorov_body_force','bodyforce_kolmogorov'}
        [stateOut, forceDiag] = resamp_apply_kolmogorov_body_force(stateOut, params, ...
            'amplitude', opts.kolmogorovForceAmplitude, ...
            'waveNumber', opts.kolmogorovForceWaveNumber, ...
            'phase', opts.kolmogorovForcePhase, ...
            'direction', opts.kolmogorovForceDirection, ...
            'zeroMeanKick', opts.kolmogorovForceZeroMeanKick, ...
            'applyDt', opts.kolmogorovForceApplyDt);
    otherwise
        error('Unknown ForceMode: %s', opts.forceMode);
end
end

function d = empty_force_diag(mode)
if nargin < 1 || isempty(mode)
    mode = 'none';
end
d = struct();
d.enabled = false;
d.mode = char(string(mode));
d.amplitude = 0;
d.waveNumber = NaN;
d.phase = NaN;
d.direction = '';
d.zeroMeanKick = false;
d.applyDt = true;
d.nActive = 0;
d.kickRms = 0;
d.kickMax = 0;
d.meanKick = [0 0];
d.momentumBefore = [NaN NaN];
d.momentumAfter = [NaN NaN];
d.momentumDelta = [0 0];
d.momentumDeltaNorm = 0;
end

function [md, metrics] = diagnostics(state, params)
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
G = md.G;
omega = periodic_vorticity(G.Ux, G.Uy, G.dx, G.dy);
metrics = struct();
metrics.omegaRms = sqrt(mean(omega(:).^2, 'omitnan'));
metrics.speedRms = sqrt(mean(G.Ux(:).^2 + G.Uy(:).^2, 'omitnan'));
metrics.meanUx = mean(G.Ux(:), 'omitnan');
metrics.meanUy = mean(G.Uy(:), 'omitnan');
metrics.omega = omega;
end

function state = initialize_velocity_memory(state, params)
activeMask = resamp_active_mask(state);
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', true, 'periodicY', true, 'activeMask', activeMask);
state.uMemUx = G.Ux;
state.uMemUy = G.Uy;
state.uMemValid = G.N > 0;
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
row.speedRms = NaN;
row.omegaRms = NaN;
row.meanUx = NaN;
row.meanUy = NaN;
row.forceEnabled = NaN;
row.forceMode = "";
row.kolmogorovForceAmplitude = NaN;
row.kolmogorovForceWaveNumber = NaN;
row.kolmogorovKickRms = NaN;
row.kolmogorovKickMax = NaN;
row.kolmogorovMomentumDeltaNorm = NaN;
row.kolmogorovMeanKickX = NaN;
row.kolmogorovMeanKickY = NaN;
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
row.remapMomentumResidualRms = NaN;
row.remapSuccessFractionNonEmpty = NaN;
end

function row = make_row(step, t, md, metrics, stepDiag, insertDiag, remapDiag)
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
row.speedRms = metrics.speedRms;
row.omegaRms = metrics.omegaRms;
row.meanUx = metrics.meanUx;
row.meanUy = metrics.meanUy;
forceDiag = get_nested(stepDiag, {'kolmogorovForce'}, struct());
row.forceEnabled = double(get_field(forceDiag, 'enabled', false));
row.forceMode = string(get_field(forceDiag, 'mode', 'none'));
row.kolmogorovForceAmplitude = get_field(forceDiag, 'amplitude', NaN);
row.kolmogorovForceWaveNumber = get_field(forceDiag, 'waveNumber', NaN);
row.kolmogorovKickRms = get_field(forceDiag, 'kickRms', NaN);
row.kolmogorovKickMax = get_field(forceDiag, 'kickMax', NaN);
row.kolmogorovMomentumDeltaNorm = get_field(forceDiag, 'momentumDeltaNorm', NaN);
meanKick = get_field(forceDiag, 'meanKick', [NaN NaN]);
if numel(meanKick) >= 2
    row.kolmogorovMeanKickX = meanKick(1);
    row.kolmogorovMeanKickY = meanKick(2);
end
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
row.remapMomentumResidualRms = get_field(remapDiag, 'momentumResidualRms', NaN);
row.remapSuccessFractionNonEmpty = get_field(remapDiag, 'successFractionNonEmpty', NaN);
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

function d = empty_remap_diag()
d = struct('massResidualRelRms', NaN, 'momentumResidualRms', NaN, 'successFractionNonEmpty', NaN);
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
if ~isfinite(nNow), nNow = 0; end
if ~isfinite(cNow), cNow = 0; end
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

function info = empty_thermostat_after_remap_diag()
info = struct('enabled', false, 'targetKBT', NaN, 'strength', NaN, ...
    'minParticlesPerCell', NaN, 'maxScale', NaN, 'nThermostattedCells', 0, ...
    'meanScale', NaN, 'minScale', NaN, 'maxScaleApplied', NaN, ...
    'meanKBTBefore', NaN, 'meanKBTAfter', NaN, 'rmsVelocityChange', 0.0);
end

function d = attach_thermostat_after_remap(d, thermostatInfo)
d.thermostatAfterRemap = thermostatInfo;
end

function omega = periodic_vorticity(Ux, Uy, dx, dy)
dUyDx = (circshift(Uy, [-1,0]) - circshift(Uy, [1,0])) / (2*dx);
dUxDy = (circshift(Ux, [0,-1]) - circshift(Ux, [0,1])) / (2*dy);
omega = dUyDx - dUxDy;
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
