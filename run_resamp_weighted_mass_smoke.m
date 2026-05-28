function out = run_resamp_weighted_mass_smoke(varargin)
%RUN_RESAMP_WEIGHTED_MASS_SMOKE Smoke run for weighted-particle SRC/Q6 infrastructure.
%
% This run is intentionally periodic and Q6-only.  It checks that the new
% state.m path can run without resampling, with m(:)=1 by default.
%
% Examples:
%   out = run_resamp_weighted_mass_smoke();
%   out = run_resamp_weighted_mass_smoke('method','weighted_q6','steps',2000);
%   out = run_resamp_weighted_mass_smoke('method','weighted_classic','steps',500);

opts = parse_options(varargin{:});
if ~isempty(opts.rngSeed)
    rng(opts.rngSeed, 'twister');
end

params = default_params(opts);
[state, initInfo] = resamp_initialize_particles_taylor_green_forced(params);

nRows = floor(opts.steps / opts.summaryEvery) + 1;
rows = repmat(empty_row(), nRows, 1);
irow = 0;

md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
G = md.G;
tg = projection_taylor_green_diagnostics(G, params);
irow = irow + 1;
rows(irow) = make_row(0, 0.0, md, tg, empty_step_diag());

for step = 1:opts.steps
    switch opts.method
        case 'weighted_classic'
            [state, stepDiag] = resamp_step_classic_periodic_weighted(state, params);
        case 'weighted_q6'
            [state, stepDiag] = resamp_step_projection_periodic_weighted(state, params);
        otherwise
            error('Unknown method: %s', opts.method);
    end

    if mod(step, opts.summaryEvery) == 0 || step == opts.steps
        md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
        G = md.G;
        tg = projection_taylor_green_diagnostics(G, params);
        irow = irow + 1;
        rows(irow) = make_row(step, step * params.dt, md, tg, stepDiag);
        fprintf(['resamp smoke %-16s step=%6d t=%.4g Nstd=%.4g Mstd=%.4g ', ...
                 'mRelStd=%.3g div=%.3e TGamp=%.5g kBT=%.5g\n'], ...
            opts.method, step, step*params.dt, md.NStd, md.MStd, md.mParticleRelStd, ...
            get_nested(stepDiag, {'rmsDivParticleAfter'}, NaN), tg.modeAmplitude, md.kBTWeighted);
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
    csvPath = fullfile(opts.outputDir, sprintf('resamp_weighted_mass_smoke_%s.csv', opts.method));
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
opts.outputDir = fullfile('runs', 'resamp_weighted_mass_smoke');

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
        otherwise
            error('Unknown option: %s', key);
    end
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
row.totalMass = NaN;
row.totalMomentumX = NaN;
row.totalMomentumY = NaN;
row.kBTWeighted = NaN;
row.kineticEnergyMeanWeighted = NaN;
row.tgAmplitude = NaN;
row.tgCoherence = NaN;
row.tgModeEnergy = NaN;
row.rmsDivParticleAfter = NaN;
row.dvAppliedRms = NaN;
row.weightedMomentumCorrectionResidual = NaN;
row.collisionDeltaPNorm = NaN;
end

function row = make_row(step, t, md, tg, stepDiag)
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
row.totalMass = md.totalMass;
row.totalMomentumX = md.totalMomentum(1);
row.totalMomentumY = md.totalMomentum(2);
row.kBTWeighted = md.kBTWeighted;
row.kineticEnergyMeanWeighted = md.kineticEnergyMeanWeighted;
row.tgAmplitude = tg.modeAmplitude;
row.tgCoherence = tg.modeCoherence;
row.tgModeEnergy = tg.modeEnergy;
row.rmsDivParticleAfter = get_nested(stepDiag, {'rmsDivParticleAfter'}, NaN);
row.dvAppliedRms = get_nested(stepDiag, {'dvAppliedRms'}, NaN);
row.weightedMomentumCorrectionResidual = get_nested(stepDiag, {'momentumCorrection','residualDeltaPNorm'}, NaN);
row.collisionDeltaPNorm = get_nested(stepDiag, {'classic','collisionDeltaPNorm'}, ...
    get_nested(stepDiag, {'collisionDeltaPNorm'}, NaN));
end

function val = get_nested(s, path, defaultValue)
val = defaultValue;
if ~isstruct(s)
    return;
end
cur = s;
for i = 1:numel(path)
    key = path{i};
    if isstruct(cur) && isfield(cur, key)
        cur = cur.(key);
    else
        return;
    end
end
if isnumeric(cur) && isscalar(cur)
    val = cur;
elseif islogical(cur) && isscalar(cur)
    val = double(cur);
end
end

function d = empty_step_diag()
d = struct();
d.rmsDivParticleAfter = NaN;
d.dvAppliedRms = NaN;
d.momentumCorrection = struct('residualDeltaPNorm', NaN);
d.collisionDeltaPNorm = NaN;
end
