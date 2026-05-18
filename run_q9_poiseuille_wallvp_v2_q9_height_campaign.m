function results = run_q9_poiseuille_wallvp_v2_q9_height_campaign(varargin)
%RUN_Q9_POISEUILLE_WALLVP_V2_Q9_HEIGHT_CAMPAIGN
% Overnight campaign for wall virtual particles v2 including Q9.
%
% Default campaign, with densityFactor=0.8 and visual monitoring ON:
%   Long/reference grid:
%     1) CLASSIC, 32 x 48, 30000 steps
%     2) Q9/general_bc, 32 x 48, 30000 steps
%
%   Height convergence, same wallVP-v2 settings:
%     3) CLASSIC, 48 x 64, 30000 steps
%     4) Q9/general_bc, 48 x 64, 30000 steps
%     5) CLASSIC, 64 x 96, 30000 steps
%     6) Q9/general_bc, 64 x 96, 30000 steps
%
% The 32 x 48 long runs are also the first points of the height series.
% The extra height block therefore starts at 48 x 64 by default.
%
% Typical use:
%   results = run_q9_poiseuille_wallvp_v2_q9_height_campaign();
%
% Quick smoke test:
%   results = run_q9_poiseuille_wallvp_v2_q9_height_campaign( ...
%       'longSteps', 2000, 'heightSteps', 2000, 'visualEvery', 500);
%
% Main summary fields:
%   ratioAccelRecent, ratioAccelLong, nuEff, R2, SNR,
%   finalMeanUx, finalCenterWall, finalLowK, finalKBT.

p = inputParser;
addParameter(p, 'outputRoot', '');
addParameter(p, 'densityFactor', 0.8);
addParameter(p, 'geometryMode', 'shifted_solid_fraction');
addParameter(p, 'longMethods', {'classic','q9'});
addParameter(p, 'heightMethods', {'classic','q9'});
addParameter(p, 'longGrid', [32 48]);
addParameter(p, 'heightGrids', [48 64; 64 96]);
addParameter(p, 'gamma', 20);
addParameter(p, 'longSteps', 30000);
addParameter(p, 'heightSteps', 30000);
addParameter(p, 'sampleEvery', 100);
addParameter(p, 'progressEvery', 1000);
addParameter(p, 'visualEnable', true);
addParameter(p, 'visualEvery', 500);
addParameter(p, 'bodyForceX', 0.005);
addParameter(p, 'wallModeY', 'bounceback');
addParameter(p, 'nuGuess', 0.032);
addParameter(p, 'fitWindowTime', 5);
addParameter(p, 'accelerationWindowTime', 5);
addParameter(p, 'nuReferenceNy', 64);
addParameter(p, 'nuScalingEnable', true);
addParameter(p, 'physicalLy', []);
addParameter(p, 'seed', 11); %#ok<NASGU> % kept for metadata/compatibility if the runner uses it later
addParameter(p, 'saveState', true);
parse(p, varargin{:});
opt = p.Results;

opt.longMethods = normalize_method_list(opt.longMethods);
opt.heightMethods = normalize_method_list(opt.heightMethods);

if isempty(opt.outputRoot)
    opt.outputRoot = fullfile(pwd, sprintf('q9_wallvp_v2_q9_height_%s', datestr(now, 'yyyymmdd_HHMMSS')));
end
if ~exist(opt.outputRoot, 'dir')
    mkdir(opt.outputRoot);
end

logFile = fullfile(opt.outputRoot, 'console_log.txt');
diary(logFile);
c = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n=== wallVP-v2 Q9-inclusive long + height campaign ===\n');
fprintf('outputRoot        : %s\n', opt.outputRoot);
fprintf('densityFactor     : %.6g\n', opt.densityFactor);
fprintf('geometryMode      : %s\n', char(string(opt.geometryMode)));
fprintf('bodyForceX        : %.6g\n', opt.bodyForceX);
fprintf('nuGuess           : %.6g\n', opt.nuGuess);
fprintf('nuReferenceNy     : %d\n', opt.nuReferenceNy);
fprintf('nuScalingEnable   : %d\n', logical(opt.nuScalingEnable));
fprintf('visualEnable      : %d, visualEvery=%d\n', logical(opt.visualEnable), opt.visualEvery);
fprintf('longMethods       : %s\n', strjoin(opt.longMethods, ', '));
fprintf('longGrid          : %d x %d, steps=%d\n', opt.longGrid(1), opt.longGrid(2), opt.longSteps);
fprintf('heightMethods     : %s\n', strjoin(opt.heightMethods, ', '));
fprintf('heightSteps       : %d\n', opt.heightSteps);
fprintf('heightGrids       :');
for ig = 1:size(opt.heightGrids,1)
    fprintf(' %dx%d', opt.heightGrids(ig,1), opt.heightGrids(ig,2));
end
fprintf('\n\n');

cases = build_cases(opt);
rows = repmat(empty_summary_row(), numel(cases), 1);
caseOutputs = struct();

for ic = 1:numel(cases)
    caseDef = cases{ic};
    fprintf('\n######################################################################\n');
    fprintf('Case %d/%d: %s\n', ic, numel(cases), caseDef.label);
    fprintf('phase=%s, method=%s, grid=%dx%d, steps=%d\n', ...
        caseDef.phase, caseDef.method, caseDef.Nx, caseDef.Ny, caseDef.nSteps);
    fprintf('######################################################################\n');

    [row, out] = run_one_case(caseDef, opt);
    rows(ic) = row;

    safeName = matlab.lang.makeValidName(caseDef.label);
    caseOutputs.(safeName) = out;

    summaryTable = struct2table(rows(1:ic)); %#ok<NASGU>
    writetable(summaryTable, fullfile(opt.outputRoot, 'summary.csv'));
    save(fullfile(opt.outputRoot, 'summary.mat'), 'summaryTable', 'opt');
end

summaryTable = struct2table(rows);
writetable(summaryTable, fullfile(opt.outputRoot, 'summary.csv'));

results = struct();
results.outputRoot = opt.outputRoot;
results.logFile = logFile;
results.options = opt;
results.summary = summaryTable;
results.cases = cases;
results.outputs = caseOutputs;

save(fullfile(opt.outputRoot, 'campaign_results.mat'), 'results', '-v7.3');

fprintf('\n=== Campaign complete ===\n');
fprintf('Summary: %s\n', fullfile(opt.outputRoot, 'summary.csv'));
fprintf('MAT    : %s\n', fullfile(opt.outputRoot, 'campaign_results.mat'));
disp(summaryTable);
end

function cases = build_cases(opt)
cases = {};

% Long 32 x 48 runs: both classic and Q9 by default.
for im = 1:numel(opt.longMethods)
    method = lower(char(string(opt.longMethods{im})));
    c = base_case(opt);
    c.phase = 'long32';
    c.method = method;
    c.Nx = opt.longGrid(1);
    c.Ny = opt.longGrid(2);
    c.nSteps = opt.longSteps;
    c.label = sprintf('long32_%s_N%dx%d_df%03d', method, c.Nx, c.Ny, round(1000*opt.densityFactor));
    cases{end+1} = c; %#ok<AGROW>
end

% Height convergence beyond the long grid: both classic and Q9 by default.
for ig = 1:size(opt.heightGrids,1)
    for im = 1:numel(opt.heightMethods)
        method = lower(char(string(opt.heightMethods{im})));
        c = base_case(opt);
        c.phase = 'height';
        c.method = method;
        c.Nx = opt.heightGrids(ig,1);
        c.Ny = opt.heightGrids(ig,2);
        c.nSteps = opt.heightSteps;
        c.label = sprintf('height_%s_N%dx%d_df%03d', method, c.Nx, c.Ny, round(1000*opt.densityFactor));
        cases{end+1} = c; %#ok<AGROW>
    end
end
end

function c = base_case(opt)
c = struct();
c.phase = '';
c.method = 'classic';
c.Nx = NaN;
c.Ny = NaN;
c.gamma = opt.gamma;
c.nSteps = NaN;
c.sampleEvery = opt.sampleEvery;
c.progressEvery = opt.progressEvery;
c.visualEnable = opt.visualEnable;
c.visualEvery = opt.visualEvery;
c.bodyForceX = opt.bodyForceX;
c.wallModeY = opt.wallModeY;
c.geometryMode = opt.geometryMode;
c.densityFactor = opt.densityFactor;
c.nuGuess = opt.nuGuess;
c.fitWindowTime = opt.fitWindowTime;
c.accelerationWindowTime = opt.accelerationWindowTime;
c.nuReferenceNy = opt.nuReferenceNy;
c.nuScalingEnable = opt.nuScalingEnable;
c.physicalLy = opt.physicalLy;
c.saveState = opt.saveState;
c.label = '';
end

function [row, out] = run_one_case(caseDef, opt)
row = empty_summary_row();
row.label = {caseDef.label};
row.phase = {caseDef.phase};
row.method = {caseDef.method};
row.Nx = caseDef.Nx;
row.Ny = caseDef.Ny;
row.gamma = caseDef.gamma;
row.nSteps = caseDef.nSteps;
row.densityFactor = caseDef.densityFactor;
row.bodyForceX = caseDef.bodyForceX;
row.nuGuess = caseDef.nuGuess;
row.outputMat = {''};
row.errorMessage = {''};
out = [];

try
    res = run_q9_poiseuille_wall_virtual_particles_validation( ...
        'methods', {caseDef.method}, ...
        'Nx', caseDef.Nx, ...
        'Ny', caseDef.Ny, ...
        'gamma', caseDef.gamma, ...
        'nSteps', caseDef.nSteps, ...
        'sampleEvery', caseDef.sampleEvery, ...
        'progressEvery', caseDef.progressEvery, ...
        'visualEnable', caseDef.visualEnable, ...
        'visualEvery', caseDef.visualEvery, ...
        'bodyForceX', caseDef.bodyForceX, ...
        'wallModeY', caseDef.wallModeY, ...
        'geometryMode', caseDef.geometryMode, ...
        'densityFactor', caseDef.densityFactor, ...
        'initialPoiseuille', true, ...
        'nuGuess', caseDef.nuGuess, ...
        'fitWindowTime', caseDef.fitWindowTime, ...
        'accelerationWindowTime', caseDef.accelerationWindowTime, ...
        'nuReferenceNy', caseDef.nuReferenceNy, ...
        'nuScalingEnable', caseDef.nuScalingEnable, ...
        'physicalLy', caseDef.physicalLy);

    methodField = lower(char(string(caseDef.method)));
    out = res.(methodField);

    row.ok = true;
    row.actualStep = get_field(out, 'actualLastStep', NaN);
    row.elapsed = get_field(out, 'elapsedWallClock', NaN);
    row.finalMeanUx = last_value(out.diagTable, 'meanUx');
    row.finalCenterWall = last_value(out.diagTable, 'centerMinusWall');
    row.finalLowK = last_value(out.diagTable, 'lowKDensityEnergy');
    row.finalKBT = last_value(out.diagTable, 'kBTCell');
    row.meanCenterWallFitWindow = get_field(out.summary, 'meanCenterMinusWall', NaN);
    row.nuEff = get_field(out.viscosity, 'nuEff', NaN);
    row.nuEffRaw = get_field(out.viscosity, 'nuEffRaw', row.nuEff);
    row.nuEffScaledNyRef = get_field(out.viscosity, 'nuEffScaledNyRef', NaN);
    row.nuEffCellY = get_field(out.viscosity, 'nuEffCellY', NaN);
    row.nuScaleFactorNyRef = get_field(out.viscosity, 'nuScaleFactorNyRef', NaN);
    row.nuReferenceNy = get_field(out.viscosity, 'nuReferenceNy', NaN);
    row.dy = get_field(out.viscosity, 'dy', NaN);
    row.R2 = get_field(out.viscosity, 'R2', NaN);
    row.SNR = get_field(out.viscosity, 'signalToNoise', NaN);
    row.fitT0 = get_field(out.viscosity, 'tMin', NaN);
    row.fitT1 = get_field(out.viscosity, 'tMax', NaN);

    if isfield(out, 'accelerationRecent')
        acc = out.accelerationRecent;
    else
        acc = poiseuille_acceleration_diagnostics(out, caseDef.accelerationWindowTime);
    end
    row.slopeRecent = acc.slopeMeanUx;
    row.ratioAccelRecent = acc.ratioAccel;
    row.accelT0 = acc.t0;
    row.accelT1 = acc.t1;

    if isfield(out, 'accelerationLong')
        accLong = out.accelerationLong;
    else
        accLong = poiseuille_acceleration_diagnostics(out, min(10, max(out.sampleTimes)-min(out.sampleTimes)));
    end
    row.slopeLong = accLong.slopeMeanUx;
    row.ratioAccelLong = accLong.ratioAccel;

    matFile = fullfile(opt.outputRoot, [caseDef.label '.mat']);
    row.outputMat = {matFile};
    if caseDef.saveState
        save(matFile, 'out', 'caseDef', 'row', '-v7.3');
    else
        outLight = out;
        if isfield(outLight, 'state')
            outLight = rmfield(outLight, 'state');
        end
        save(matFile, 'outLight', 'caseDef', 'row', '-v7.3');
    end

catch ME
    row.ok = false;
    row.errorMessage = {ME.message};
    errFile = fullfile(opt.outputRoot, [caseDef.label '_ERROR.txt']);
    fid = fopen(errFile, 'w');
    if fid > 0
        fprintf(fid, '%s\n\n', ME.message);
        for k = 1:numel(ME.stack)
            fprintf(fid, '  in %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
        fclose(fid);
    end
    fprintf(2, '\nERROR in case %s: %s\n', caseDef.label, ME.message);
end
end

function row = empty_summary_row()
row = struct( ...
    'label', {{''}}, 'phase', {{''}}, 'method', {{''}}, 'ok', false, ...
    'Nx', NaN, 'Ny', NaN, 'gamma', NaN, 'nSteps', NaN, ...
    'actualStep', NaN, 'elapsed', NaN, ...
    'densityFactor', NaN, 'bodyForceX', NaN, 'nuGuess', NaN, ...
    'nuEff', NaN, 'nuEffRaw', NaN, 'nuEffScaledNyRef', NaN, ...
    'nuEffCellY', NaN, 'nuScaleFactorNyRef', NaN, 'nuReferenceNy', NaN, 'dy', NaN, ...
    'R2', NaN, 'SNR', NaN, 'fitT0', NaN, 'fitT1', NaN, ...
    'slopeRecent', NaN, 'ratioAccelRecent', NaN, 'accelT0', NaN, 'accelT1', NaN, ...
    'slopeLong', NaN, 'ratioAccelLong', NaN, ...
    'finalMeanUx', NaN, 'finalCenterWall', NaN, 'meanCenterWallFitWindow', NaN, ...
    'finalLowK', NaN, 'finalKBT', NaN, ...
    'outputMat', {{''}}, 'errorMessage', {{''}});
end

function val = last_value(T, name)
val = NaN;
if istable(T) && ismember(name, T.Properties.VariableNames) && height(T) >= 1
    x = T.(name);
    val = x(end);
end
end

function val = get_field(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
end
end

function methods = normalize_method_list(x)
if ischar(x) || isstring(x)
    methods = cellstr(x);
elseif iscell(x)
    methods = cell(size(x));
    for i = 1:numel(x)
        methods{i} = char(string(x{i}));
    end
else
    error('Method lists must be char, string, or cell array.');
end
methods = reshape(methods, 1, []);
for i = 1:numel(methods)
    methods{i} = lower(strrep(strtrim(methods{i}), '-', '_'));
    if ~ismember(methods{i}, {'classic','q6','q9'})
        error('Unknown method "%s". Expected classic, q6, or q9.', methods{i});
    end
end
end
