function out = run_resamp_tg_physics_validation(varargin)
%RUN_RESAMP_TG_PHYSICS_VALIDATION Forced Taylor--Green validation for weighted resampling.
%
% This runner keeps the domain fully wet and periodic.  It is intended to
% calibrate the hydrodynamic behaviour of the weighted/resampled SRC/MPCD
% formulation before adding walls, interfaces or solids.
%
% Example:
%   out = run_resamp_tg_physics_validation('steps',3000,'Nx',64,'Ny',64);

opts = parse_options(varargin{:});
if ~exist(opts.outputRoot, 'dir'), mkdir(opts.outputRoot); end

caseDefs = build_cases(opts);
outs = struct();
rows = repmat(empty_summary_row(), numel(caseDefs), 1);

for ic = 1:numel(caseDefs)
    c = caseDefs(ic);
    fprintf('\n=== Forced TG resampling physics case: %s ===\n', c.label);
    outCaseDir = fullfile(opts.outputRoot, c.label);
    if ~exist(outCaseDir, 'dir'), mkdir(outCaseDir); end
    runArgs = { ...
        'method', c.method, ...
        'steps', opts.steps, ...
        'summaryEvery', opts.summaryEvery, ...
        'visualEvery', opts.visualEvery, ...
        'Nx', opts.Nx, 'Ny', opts.Ny, 'gamma', opts.gamma, ...
        'dt', opts.dt, 'alphaDeg', opts.alphaDeg, 'kBT', opts.kBT, ...
        'initialDepletion', 'none', ...
        'taylorGreenInitialAmplitude', opts.taylorGreenInitialAmplitude, ...
        'taylorGreenAmplitude', opts.taylorGreenInitialAmplitude, ...
        'taylorGreenForcingEnable', true, ...
        'taylorGreenForcingAmplitude', opts.taylorGreenForcingAmplitude, ...
        'projectionStrength', c.projectionStrength, ...
        'projectionInterpolationMethod', opts.projectionInterpolationMethod, ...
        'capacityFactor', opts.capacityFactor, ...
        'NMin', c.NMin, 'NTarget', opts.NTarget, 'NMax', c.NMax, ...
        'ExtractEvery', c.extractEvery, 'InsertEvery', c.insertEvery, ...
        'RemapEvery', c.remapEvery, ...
        'ThermostatAfterStep', c.thermostatAfterStep, ...
        'ThermostatAfterProjection', false, ...
        'ThermostatAfterRemap', c.thermostatAfterRemap, ...
        'ThermostatStrength', opts.thermostatStrength, ...
        'ThermostatTargetKBT', opts.kBT, ...
        'RemapMassSafetyEnable', c.massSafetyEnable, ...
        'RemapMassSafetyMinFactor', opts.massSafetyMinFactor, ...
        'RemapMassSafetyMaxFactor', opts.massSafetyMaxFactor, ...
        'writeCsv', true, ...
        'outputDir', outCaseDir, ...
        'figureId', opts.figureIdBase + ic - 1, ...
        'showDebugFigure', false, ...
        'saveFrames', opts.saveFrames, ...
        'frameDir', fullfile(outCaseDir, 'frames'), ...
        'framePrefix', sprintf('tg_%s', c.label), ...
        'rngSeed', opts.rngSeed + ic - 1};
    outCase = run_resamp_pool_insertion_visual_demo(runArgs{:});
    outs.(c.label) = outCase;
    rows(ic) = summarize_case(c, outCase, opts);
end

summary = struct2table(rows);
writetable(summary, fullfile(opts.outputRoot, 'resamp_tg_physics_summary.csv'));

out = struct();
out.options = opts;
out.caseDefs = caseDefs;
out.cases = outs;
out.summary = summary;
out.outputRoot = opts.outputRoot;
save(fullfile(opts.outputRoot, 'resamp_tg_physics_validation.mat'), 'out', '-v7.3');

fprintf('\nForced TG weighted-resampling summary:\n');
disp(summary);
fprintf('Wrote %s\n', fullfile(opts.outputRoot, 'resamp_tg_physics_summary.csv'));
end

function caseDefs = build_cases(opts)
caseDefs = struct([]);
caseDefs(1).label = 'classic_reference';
caseDefs(1).method = 'weighted_classic';
caseDefs(1).projectionStrength = 0;
caseDefs(1).extractEvery = 0;
caseDefs(1).insertEvery = 0;
caseDefs(1).remapEvery = 0;
caseDefs(1).thermostatAfterStep = true;
caseDefs(1).thermostatAfterRemap = false;
caseDefs(1).massSafetyEnable = false;
caseDefs(1).NMin = 0;
caseDefs(1).NMax = 2*opts.NTarget;

caseDefs(2).label = 'q6_resampled';
caseDefs(2).method = 'weighted_q6';
caseDefs(2).projectionStrength = 1;
caseDefs(2).extractEvery = 1;
caseDefs(2).insertEvery = 1;
caseDefs(2).remapEvery = 1;
caseDefs(2).thermostatAfterStep = false;
caseDefs(2).thermostatAfterRemap = true;
caseDefs(2).massSafetyEnable = true;
caseDefs(2).NMin = opts.NMin;
caseDefs(2).NMax = opts.NMax;
end

function row = summarize_case(c, outCase, opts)
T = outCase.summary;
last = T(end,:);
if isfield(outCase, 'finalTaylorGreenDiagnostics')
    tg = outCase.finalTaylorGreenDiagnostics;
else
    tg = struct();
end
fit = struct('status','not_run','nuPreferred',NaN,'preferredMethod','','r2DerivativeKnownF',NaN,'r2ExpKnownF',NaN,'tStart',NaN,'tEnd',NaN,'nFit',NaN);
try
    if ismember('tgAmplitude', T.Properties.VariableNames)
        pfit = outCase.params;
        pfit.taylorGreenForceAmplitude = opts.taylorGreenForcingAmplitude;
        fit = projection_fit_forced_taylor_green_viscosity(T.t, T.tgAmplitude, pfit, 'fitFraction', opts.fitFraction);
    end
catch ME
    fit.status = ['failed: ' ME.message];
end
row = empty_summary_row();
row.label = {c.label};
row.method = {c.method};
row.steps = opts.steps;
row.finalStep = last.step;
row.finalTime = last.t;
row.finalTGAmplitude = get_table_scalar(last, 'tgAmplitude');
row.finalTGCoherence = get_struct_field(tg, 'modeCoherence', NaN);
row.finalEnstrophy = get_struct_field(tg, 'enstrophy', NaN);
row.finalOmegaRms = get_struct_field(tg, 'omegaRms', NaN);
row.finalDensityRelRms = get_struct_field(tg, 'densityRelRMS', NaN);
row.finalRmsDiv = get_struct_field(tg, 'rmsDiv', get_table_scalar(last, 'rmsDivParticleAfter'));
row.finalNMin = get_table_scalar(last, 'NMin');
row.finalNMax = get_table_scalar(last, 'NMax');
row.finalMRelRms = get_table_scalar(last, 'MRelRms');
row.finalKBT = get_table_scalar(last, 'kBTWeighted');
row.finalMassRelStd = get_table_scalar(last, 'mParticleRelStd');
row.insertedCum = get_table_scalar(last, 'insertedParticlesCumulative');
row.extractedCum = get_table_scalar(last, 'extractedParticlesCumulative');
row.massSafetyCellsFinal = get_table_scalar(last, 'remapMassSafetyCells');
row.nuEffTG = get_struct_field(fit, 'nuPreferred', NaN);
row.nuEffTGMethod = {get_struct_field(fit, 'preferredMethod', '')};
row.nuEffTGStatus = {get_struct_field(fit, 'status', '')};
row.nuFitR2DerivativeKnownF = get_struct_field(fit, 'r2DerivativeKnownF', NaN);
row.nuFitR2ExpKnownF = get_struct_field(fit, 'r2ExpKnownF', NaN);
row.nuFitT0 = get_struct_field(fit, 'tStart', NaN);
row.nuFitT1 = get_struct_field(fit, 'tEnd', NaN);
row.finalModeEnergy = get_struct_field(tg, 'modeEnergy', NaN);
end

function row = empty_summary_row()
row = struct();
row.label = {''}; row.method = {''}; row.steps = NaN; row.finalStep = NaN; row.finalTime = NaN;
row.finalTGAmplitude = NaN; row.finalTGCoherence = NaN; row.finalEnstrophy = NaN; row.finalOmegaRms = NaN;
row.finalDensityRelRms = NaN; row.finalRmsDiv = NaN; row.finalNMin = NaN; row.finalNMax = NaN;
row.finalMRelRms = NaN; row.finalKBT = NaN; row.finalMassRelStd = NaN;
row.insertedCum = NaN; row.extractedCum = NaN; row.massSafetyCellsFinal = NaN;
row.nuEffTG = NaN; row.nuEffTGMethod = {''}; row.nuEffTGStatus = {''};
row.nuFitR2DerivativeKnownF = NaN; row.nuFitR2ExpKnownF = NaN; row.nuFitT0 = NaN; row.nuFitT1 = NaN;
row.finalModeEnergy = NaN;
end

function v = get_table_scalar(row, name)
if ismember(name, row.Properties.VariableNames)
    v = row.(name)(1);
else
    v = NaN;
end
end

function v = get_struct_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end

function opts = parse_options(varargin)
opts = struct();
opts.outputRoot = fullfile('runs','resamp_tg_physics_validation');
opts.Nx = 64; opts.Ny = 64; opts.gamma = 20;
opts.NTarget = 20; opts.NMin = 14; opts.NMax = 26;
opts.steps = 3000; opts.summaryEvery = 50; opts.visualEvery = 0;
opts.dt = 1.0e-3; opts.alphaDeg = 90; opts.kBT = 0.01;
opts.taylorGreenInitialAmplitude = 0.10;
opts.taylorGreenForcingAmplitude = 0.12;
opts.thermostatStrength = 0.25;
opts.capacityFactor = 2.0;
opts.massSafetyMinFactor = 0.25;
opts.massSafetyMaxFactor = 4.0;
opts.projectionInterpolationMethod = 'nearest';
opts.rngSeed = 12345;
opts.figureIdBase = 900;
opts.saveFrames = false;
opts.fitFraction = 0.60;
if mod(numel(varargin),2) ~= 0, error('Options must be name/value pairs.'); end
for k=1:2:numel(varargin)
    key = lower(char(string(varargin{k}))); val = varargin{k+1};
    switch key
        case 'outputroot', opts.outputRoot = char(string(val));
        case 'nx', opts.Nx = val;
        case 'ny', opts.Ny = val;
        case 'gamma', opts.gamma = val; opts.NTarget = val;
        case 'ntarget', opts.NTarget = val;
        case 'nmin', opts.NMin = val;
        case 'nmax', opts.NMax = val;
        case 'steps', opts.steps = val;
        case 'summaryevery', opts.summaryEvery = val;
        case 'visualevery', opts.visualEvery = val;
        case 'dt', opts.dt = val;
        case 'alphadeg', opts.alphaDeg = val;
        case 'kbt', opts.kBT = val;
        case 'taylorgreeninitialamplitude', opts.taylorGreenInitialAmplitude = val;
        case {'taylorgreenforcingamplitude','forceamplitude'}, opts.taylorGreenForcingAmplitude = val;
        case 'thermostatstrength', opts.thermostatStrength = val;
        case 'capacityfactor', opts.capacityFactor = val;
        case 'masssafetyminfactor', opts.massSafetyMinFactor = val;
        case 'masssafetymaxfactor', opts.massSafetyMaxFactor = val;
        case 'projectioninterpolationmethod', opts.projectionInterpolationMethod = char(string(val));
        case 'rngseed', opts.rngSeed = val;
        case 'figureidbase', opts.figureIdBase = val;
        case 'saveframes', opts.saveFrames = logical(val);
        case 'fitfraction', opts.fitFraction = val;
        otherwise, error('Unknown option: %s', key);
    end
end
end
