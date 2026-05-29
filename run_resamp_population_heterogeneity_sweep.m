function out = run_resamp_population_heterogeneity_sweep(varargin)
%RUN_RESAMP_POPULATION_HETEROGENEITY_SWEEP Compact sweep over initial heterogeneity levels.
%
% The sweep is intentionally small and non-visual.  It is meant to answer one
% question: for a controlled distributed initial N_c heterogeneity, does the
% closed-loop extraction/insertion/remap system converge to a bounded state?
%
% Example:
%   out = run_resamp_population_heterogeneity_sweep( ...
%       'PopulationHeterogeneityStdList',[0 2 4 6 8], 'steps',300);

opts = parse_options(varargin{:});
if ~exist(opts.outputRoot, 'dir')
    mkdir(opts.outputRoot);
end

nCases = numel(opts.populationHeterogeneityStdList) * numel(opts.seeds);
rows = repmat(empty_row(), nCases, 1);
icase = 0;
for iseed = 1:numel(opts.seeds)
    for ilev = 1:numel(opts.populationHeterogeneityStdList)
        level = opts.populationHeterogeneityStdList(ilev);
        seed = opts.seeds(iseed);
        icase = icase + 1;
        caseLabel = sprintf('hetero%g_seed%d', level, seed);
        caseLabel = regexprep(caseLabel, '[^A-Za-z0-9_\-]', 'p');
        caseDir = fullfile(opts.outputRoot, caseLabel);
        fprintf('Running heterogeneity case %d/%d: std=%g seed=%d\n', icase, nCases, level, seed);
        outCase = run_resamp_pool_insertion_visual_demo( ...
            'method', opts.method, ...
            'steps', opts.steps, ...
            'visualEvery', 0, ...
            'summaryEvery', opts.summaryEvery, ...
            'Nx', opts.Nx, ...
            'Ny', opts.Ny, ...
            'gamma', opts.gamma, ...
            'dt', opts.dt, ...
            'initialDepletion','heterogeneous_population', ...
            'PopulationHeterogeneityStd', level, ...
            'PopulationHeterogeneityMin', opts.populationHeterogeneityMin, ...
            'PopulationHeterogeneityMax', opts.populationHeterogeneityMax, ...
            'NMin', opts.NMin, ...
            'NTarget', opts.NTarget, ...
            'NMax', opts.NMax, ...
            'ExtractEvery', opts.extractEvery, ...
            'InsertEvery', opts.insertEvery, ...
            'ExtractSelectionMode', opts.extractSelectionMode, ...
            'RemapMassSafetyEnable', opts.remapMassSafetyEnable, ...
            'RemapMassSafetyMode', opts.remapMassSafetyMode, ...
            'RemapMassSafetyMinFactor', opts.remapMassSafetyMinFactor, ...
            'RemapMassSafetyMaxFactor', opts.remapMassSafetyMaxFactor, ...
            'ThermostatAfterRemap', opts.thermostatAfterRemap, ...
            'ThermostatStrength', opts.thermostatStrength, ...
            'ThermostatTargetKBT', opts.thermostatTargetKBT, ...
            'rngSeed', seed, ...
            'writeCsv', opts.writeCaseCsv, ...
            'outputDir', caseDir);
        rows(icase) = summarize_case(outCase, level, seed, caseLabel, caseDir);
    end
end

summary = struct2table(rows);
summaryPath = fullfile(opts.outputRoot, 'heterogeneity_sweep_summary.csv');
writetable(summary, summaryPath);
fprintf('Wrote %s\n', summaryPath);

out = struct();
out.options = opts;
out.summary = summary;
out.summaryPath = summaryPath;
end

function row = empty_row()
row = struct();
row.caseLabel = "";
row.outputDir = "";
row.populationHeterogeneityStd = NaN;
row.seed = NaN;
row.initialTargetNStd = NaN;
row.initialNMin = NaN;
row.initialNMax = NaN;
row.initialNStd = NaN;
row.initialMRelRms = NaN;
row.initialMovedParticles = NaN;
row.finalStep = NaN;
row.NpActiveFinal = NaN;
row.NfreeFinal = NaN;
row.NMinFinal = NaN;
row.NMaxFinal = NaN;
row.NStdFinal = NaN;
row.extractedCumulative = NaN;
row.insertedCumulative = NaN;
row.lastExtractionStep = NaN;
row.lastInsertionStep = NaN;
row.MRelRmsFinal = NaN;
row.mParticleMinFinal = NaN;
row.mParticleMaxFinal = NaN;
row.remapMassSafetyCellsFinal = NaN;
row.remapMassSafetyParticlesFinal = NaN;
row.remapMassSafetyInfeasibleCellsFinal = NaN;
row.remapMassSafetyCandidateMinFactorFinal = NaN;
row.remapMassSafetyCandidateMaxFactorFinal = NaN;
row.mParticleRelStdFinal = NaN;
row.kBTFinal = NaN;
row.rmsDivAfterFinal = NaN;
end

function row = summarize_case(outCase, level, seed, caseLabel, caseDir)
row = empty_row();
row.caseLabel = string(caseLabel);
row.outputDir = string(caseDir);
row.populationHeterogeneityStd = level;
row.seed = seed;
info = outCase.initialPopulationEditInfo;
row.initialTargetNStd = get_field(info, 'targetNStd', NaN);
row.initialNMin = get_field(info, 'NMinAfter', NaN);
row.initialNMax = get_field(info, 'NMaxAfter', NaN);
row.initialNStd = get_field(info, 'NStdAfter', NaN);
row.initialMRelRms = get_field(info, 'MRelRmsAfter', NaN);
row.initialMovedParticles = get_field(info, 'nMovedParticles', NaN);
T = outCase.summary;
last = T(end,:);
row.finalStep = last.step;
row.NpActiveFinal = last.NpActive;
row.NfreeFinal = last.Nfree;
row.NMinFinal = last.NMin;
row.NMaxFinal = last.NMax;
row.NStdFinal = last.NStd;
row.extractedCumulative = last.extractedParticlesCumulative;
row.insertedCumulative = last.insertedParticlesCumulative;
row.lastExtractionStep = last.lastExtractionStep;
row.lastInsertionStep = last.lastInsertionStep;
row.MRelRmsFinal = last.MRelRms;
row.mParticleRelStdFinal = last.mParticleRelStd;
row.mParticleMinFinal = last.mParticleMin;
row.mParticleMaxFinal = last.mParticleMax;
if ismember('remapMassSafetyCells', T.Properties.VariableNames)
    row.remapMassSafetyCellsFinal = last.remapMassSafetyCells;
    row.remapMassSafetyParticlesFinal = last.remapMassSafetyParticles;
    row.remapMassSafetyInfeasibleCellsFinal = last.remapMassSafetyInfeasibleCells;
    row.remapMassSafetyCandidateMinFactorFinal = last.remapMassSafetyCandidateMinFactor;
    row.remapMassSafetyCandidateMaxFactorFinal = last.remapMassSafetyCandidateMaxFactor;
end
row.kBTFinal = last.kBTWeighted;
row.rmsDivAfterFinal = last.rmsDivParticleAfter;
end

function opts = parse_options(varargin)
opts = struct();
opts.method = 'weighted_q6';
opts.steps = 300;
opts.summaryEvery = 25;
opts.Nx = 32;
opts.Ny = 32;
opts.gamma = 20;
opts.dt = 0.001;
opts.populationHeterogeneityStdList = [0 2 4 6 8];
opts.populationHeterogeneityMin = 0;
opts.populationHeterogeneityMax = 40;
opts.NMin = 14;
opts.NTarget = 20;
opts.NMax = 26;
opts.extractEvery = 1;
opts.insertEvery = 1;
opts.extractSelectionMode = 'closest_to_cell_mean';
opts.remapMassSafetyEnable = true;
opts.remapMassSafetyMode = 'uniform_mass_velocity_shift';
opts.remapMassSafetyMinFactor = 0.25;
opts.remapMassSafetyMaxFactor = 4.0;
opts.thermostatAfterRemap = true;
opts.thermostatStrength = 0.25;
opts.thermostatTargetKBT = 0.01;
opts.seeds = 12345;
opts.writeCaseCsv = true;
opts.outputRoot = fullfile('runs','resamp_population_heterogeneity_sweep');
if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
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
        case {'populationheterogeneitystdlist','heterogeneitystdlist','levels'}
            opts.populationHeterogeneityStdList = val;
        case {'populationheterogeneitymin','heterogeneitymin'}
            opts.populationHeterogeneityMin = val;
        case {'populationheterogeneitymax','heterogeneitymax'}
            opts.populationHeterogeneityMax = val;
        case 'nmin'
            opts.NMin = val;
        case 'ntarget'
            opts.NTarget = val;
        case 'nmax'
            opts.NMax = val;
        case 'extractevery'
            opts.extractEvery = val;
        case 'insertevery'
            opts.insertEvery = val;
        case {'extractselectionmode','selectionmode'}
            opts.extractSelectionMode = lower(char(string(val)));
        case {'remapmasssafetyenable','masssafetyenable'}
            opts.remapMassSafetyEnable = logical(val);
        case {'remapmasssafetymode','masssafetymode'}
            opts.remapMassSafetyMode = lower(char(string(val)));
        case {'remapmasssafetyminfactor','masssafetyminfactor'}
            opts.remapMassSafetyMinFactor = val;
        case {'remapmasssafetymaxfactor','masssafetymaxfactor'}
            opts.remapMassSafetyMaxFactor = val;
        case 'thermostatafterremap'
            opts.thermostatAfterRemap = logical(val);
        case 'thermostatstrength'
            opts.thermostatStrength = val;
        case 'thermostattargetkbt'
            opts.thermostatTargetKBT = val;
        case {'seeds','rngseeds'}
            opts.seeds = val;
        case 'writecasecsv'
            opts.writeCaseCsv = logical(val);
        case {'outputroot','outputdir'}
            opts.outputRoot = char(string(val));
        otherwise
            error('Unknown option: %s', key);
    end
end
opts.populationHeterogeneityStdList = opts.populationHeterogeneityStdList(:).';
opts.seeds = opts.seeds(:).';
end

function v = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
