function results = run_q9_piston_full_model_validation(varargin)
%RUN_Q9_PISTON_FULL_MODEL_VALIDATION Validate full Q6/Q9 piston model.
%
% This campaign applies to the piston the same model structure validated in
% Taylor-Green and Poiseuille:
%   classic SRC/MPCD -> Q6 velocity projection -> Q9 low-k mass-flux correction.
%
% Default methods:
%   classic, q6_projection, q9_full
%
% Optional diagnostic method:
%   q9_mass_only, q9_historical
%
% Example:
%   results = run_q9_piston_full_model_validation( ...
%       'Nx',32,'Ny',32,'gamma',20,'nSteps',10000, ...
%       'compressionTarget',0.05, ...
%       'methods',{'classic','q6_projection','q9_full'}, ...
%       'massFluxProjectionOperator','general_bc', ...
%       'massFluxTargetFilter','elliptic_lowpass', ...
%       'visualEnable',true,'visualEvery',500);

opts = parse_inputs(varargin{:});
opts = set_default(opts, 'outputRoot', fullfile(pwd, ['piston_full_q9_' datestr(now,'yyyymmdd_HHMMSS')]));
opts = set_default(opts, 'methods', {'classic','q6_projection','q9_full'});
opts = set_default(opts, 'Nx', 32);
opts = set_default(opts, 'Ny', 32);
opts = set_default(opts, 'gamma', 20);
opts = set_default(opts, 'nSteps', 5000);
opts = set_default(opts, 'sampleEvery', 100);
opts = set_default(opts, 'progressEvery', 500);
opts = set_default(opts, 'compressionTarget', 0.05);
opts = set_default(opts, 'densityFactor', 0.8);
opts = set_default(opts, 'projectionStrength', 1.0);
opts = set_default(opts, 'massFluxBeta', 5.0e-4);
opts = set_default(opts, 'massFluxProjectionStrength', 1.0);
opts = set_default(opts, 'massFluxProjectionOperator', 'general_bc');
opts = set_default(opts, 'massFluxTargetFilter', 'elliptic_lowpass');
opts = set_default(opts, 'massFluxFinalVelocityProjectionCleanup', false);
opts = set_default(opts, 'massFluxFinalVelocityProjectionStrength', 0.5);
opts = set_default(opts, 'pistonQ9TargetMode', 'reference_gamma');
opts = set_default(opts, 'pistonQ9FreezeEllipticMetric', false);
opts = set_default(opts, 'compressibilityFitStartFraction', 0.2);
opts = set_default(opts, 'compressibilityFitEndFraction', 1.0);
opts = set_default(opts, 'visualEnable', true);
opts = set_default(opts, 'visualEvery', 500);
opts = set_default(opts, 'seed', 11);

if ~exist(opts.outputRoot, 'dir')
    mkdir(opts.outputRoot);
end
logFile = fullfile(opts.outputRoot, 'console_log.txt');
diary(logFile);
diary on;
cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('=== Piston full Q6/Q9 model validation ===\n');
fprintf('outputRoot       : %s\n', opts.outputRoot);
fprintf('methods          : %s\n', strjoin(cellstr(opts.methods), ', '));
fprintf('grid/gamma       : %d x %d / %.6g\n', opts.Nx, opts.Ny, opts.gamma);
fprintf('steps            : %d\n', opts.nSteps);
fprintf('compressionTarget: %.6g\n', opts.compressionTarget);
fprintf('densityFactor    : %.6g\n', opts.densityFactor);
fprintf('Q6 strength      : %.6g\n', opts.projectionStrength);
fprintf('Q9 operator/filter/beta : %s / %s / %.6g\n', ...
    char(string(opts.massFluxProjectionOperator)), char(string(opts.massFluxTargetFilter)), opts.massFluxBeta);
fprintf('Q9 cleanup       : %d / %.6g\n', opts.massFluxFinalVelocityProjectionCleanup, opts.massFluxFinalVelocityProjectionStrength);
fprintf('Q9 target mode   : %s\n', char(string(opts.pistonQ9TargetMode)));
fprintf('freeze elliptic metric: %d\n', opts.pistonQ9FreezeEllipticMetric);
fprintf('Keff fit fraction: %.3g--%.3g\n', opts.compressibilityFitStartFraction, opts.compressibilityFitEndFraction);
fprintf('visual           : %d every %d\n\n', opts.visualEnable, opts.visualEvery);

summaryRows = cell(numel(opts.methods), 1);
results = struct();
results.opts = opts;
results.outputs = struct();

for im = 1:numel(opts.methods)
    method = char(string(opts.methods{im}));
    label = sprintf('piston_%s_N%dx%d_df%03d_fullq9', method, opts.Nx, opts.Ny, round(1000*opts.densityFactor));
    fprintf('\n######################################################################\n');
    fprintf('Case %d/%d: %s\n', im, numel(opts.methods), label);
    fprintf('######################################################################\n\n');

    p = struct();
    p.method = method;
    p.Nx = opts.Nx;
    p.Ny = opts.Ny;
    p.gamma = opts.gamma;
    p.nSteps = opts.nSteps;
    p.sampleEvery = opts.sampleEvery;
    p.progressEvery = opts.progressEvery;
    p.compressionTarget = opts.compressionTarget;
    p.pistonCompressionTarget = opts.compressionTarget;
    p.wallVirtualParticlesDensityFactor = opts.densityFactor;
    p.projectionStrength = opts.projectionStrength;
    p.massFluxProjectionStrength = opts.massFluxProjectionStrength;
    p.massFluxDensityRelaxationBeta = opts.massFluxBeta;
    p.massFluxProjectionOperator = opts.massFluxProjectionOperator;
    p.massFluxTargetFilter = opts.massFluxTargetFilter;
    p.massFluxFinalVelocityProjectionCleanup = opts.massFluxFinalVelocityProjectionCleanup;
    p.massFluxFinalVelocityProjectionStrength = opts.massFluxFinalVelocityProjectionStrength;
    p.pistonQ9TargetMode = opts.pistonQ9TargetMode;
    p.pistonQ9FreezeEllipticMetric = opts.pistonQ9FreezeEllipticMetric;
    p.pistonCompressibilityFitStartFraction = opts.compressibilityFitStartFraction;
    p.pistonCompressibilityFitEndFraction = opts.compressibilityFitEndFraction;
    p.visualEnable = opts.visualEnable;
    p.visualEvery = opts.visualEvery;
    p.visualFigureId = 740 + im - 1;
    p.seed = opts.seed;

    try
        out = run_projection_piston_wallvp_demo(p);
        ok = true;
        errMsg = '';
    catch ME
        out = [];
        ok = false;
        errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
        fprintf(2, 'ERROR in %s:\n%s\n', label, errMsg);
    end

    matFile = fullfile(opts.outputRoot, [label '.mat']);
    save(matFile, 'out', 'p', '-v7.3');
    results.outputs.(matlab.lang.makeValidName(label)) = out;

    if ok && isstruct(out)
        row = make_summary_row(label, method, ok, opts, out, matFile, errMsg);
    else
        row = make_failure_row(label, method, opts, matFile, errMsg);
    end
    summaryRows{im} = row;
end

summary = vertcat(summaryRows{:});
results.summary = summary;
results.comparison = compare_to_classic(results);
writetable(summary, fullfile(opts.outputRoot, 'summary.csv'));
save(fullfile(opts.outputRoot, 'campaign_results.mat'), 'results', '-v7.3');

fprintf('\n=== Piston full Q6/Q9 campaign complete ===\n');
fprintf('Summary: %s\n', fullfile(opts.outputRoot, 'summary.csv'));
disp(summary);
if ~isempty(results.comparison)
    fprintf('\n=== Ratios versus classic ===\n');
    disp(results.comparison);
end
end

function row = make_summary_row(label, method, ok, opts, out, matFile, errMsg)
s = out.summary;
c = out.compressibility;
row = table({label}, {method}, ok, opts.Nx, opts.Ny, opts.gamma, opts.nSteps, out.actualStep, out.elapsedWallClock, ...
    opts.compressionTarget, opts.densityFactor, opts.projectionStrength, opts.massFluxBeta, ...
    {char(string(opts.massFluxProjectionOperator))}, {char(string(opts.massFluxTargetFilter))}, ...
    opts.massFluxFinalVelocityProjectionCleanup, opts.massFluxFinalVelocityProjectionStrength, ...
    {char(string(opts.pistonQ9TargetMode))}, opts.pistonQ9FreezeEllipticMetric, ...
    s.finalYTop, s.finalCompression, ...
    s.finalRhoPhysicalMean, s.meanRhoPhysicalMean, s.finalRhoPhysicalRelError, s.meanRhoPhysicalRelError, ...
    s.finalStdN, s.meanStdN, s.finalLowKDensity, s.meanLowKDensity, ...
    s.finalPkinMean, s.meanPkinMean, s.finalPkinIdealRatio, s.meanPkinIdealRatio, ...
    s.finalPkinTopLayerMean, s.meanPkinTopLayerMean, ...
    s.finalTopWallPressureImpact, s.meanTopWallPressureImpact, ...
    s.finalTopWallPressureVP, s.meanTopWallPressureVP, ...
    s.finalTopWallPressureTotal, s.meanTopWallPressureTotal, ...
    s.finalPistonPowerOnFluid, s.meanPistonPowerOnFluid, ...
    s.finalPistonWorkOnFluidCumulative, ...
    c.KeffPkin, c.R2Pkin, c.KeffPkinTopLayer, c.R2PkinTopLayer, c.KeffTopWallTotal, c.R2TopWallTotal, ...
    s.finalKBTCell, s.meanKBTCell, ...
    s.finalMassFluxBefore, s.meanMassFluxBefore, ...
    s.finalMassFluxAfter, s.meanMassFluxAfter, s.finalMassFluxResidual, s.meanMassFluxResidual, ...
    {matFile}, {errMsg}, ...
    'VariableNames', variable_names());
end

function row = make_failure_row(label, method, opts, matFile, errMsg)
vals = num2cell(nan(1, 51));
row = table({label}, {method}, false, opts.Nx, opts.Ny, opts.gamma, opts.nSteps, vals{:}, {matFile}, {errMsg}, ...
    'VariableNames', variable_names());
end

function names = variable_names()
names = {'label','method','ok','Nx','Ny','gamma','nSteps','actualStep','elapsed', ...
    'compressionTarget','densityFactor','projectionStrength','massFluxBeta', ...
    'massFluxProjectionOperator','massFluxTargetFilter','massFluxCleanup','massFluxCleanupStrength', ...
    'pistonQ9TargetMode','pistonQ9FreezeEllipticMetric', ...
    'finalYTop','finalCompression', ...
    'finalRhoPhysicalMean','meanRhoPhysicalMean','finalRhoPhysicalRelError','meanRhoPhysicalRelError', ...
    'finalStdN','meanStdN','finalLowKDensity','meanLowKDensity', ...
    'finalPkinMean','meanPkinMean','finalPkinIdealRatio','meanPkinIdealRatio', ...
    'finalPkinTopLayerMean','meanPkinTopLayerMean', ...
    'finalTopWallPressureImpact','meanTopWallPressureImpact', ...
    'finalTopWallPressureVP','meanTopWallPressureVP', ...
    'finalTopWallPressureTotal','meanTopWallPressureTotal', ...
    'finalPistonPowerOnFluid','meanPistonPowerOnFluid','finalPistonWorkOnFluidCumulative', ...
    'KeffPkin','R2Pkin','KeffPkinTopLayer','R2PkinTopLayer','KeffTopWallTotal','R2TopWallTotal', ...
    'finalKBTCell','meanKBTCell', ...
    'finalMassFluxBefore','meanMassFluxBefore','finalMassFluxAfter','meanMassFluxAfter','finalMassFluxResidual','meanMassFluxResidual', ...
    'outputMat','errorMessage'};
end

function C = compare_to_classic(results)
C = table();
if ~isfield(results, 'outputs') || ~isstruct(results.outputs)
    return;
end
names = fieldnames(results.outputs);
classicOut = [];
for i = 1:numel(names)
    out = results.outputs.(names{i});
    if isstruct(out) && isfield(out, 'params') && strcmpi(char(string(out.params.method)), 'classic')
        classicOut = out;
        break;
    end
end
if isempty(classicOut)
    return;
end
S0 = classicOut.summary;
rows = {};
for i = 1:numel(names)
    out = results.outputs.(names{i});
    if ~isstruct(out) || ~isfield(out, 'params') || ~isfield(out, 'summary')
        continue;
    end
    method = char(string(out.params.method));
    if strcmpi(method, 'classic')
        continue;
    end
    S = out.summary;
    rows(end+1,:) = {method, ... %#ok<AGROW>
        S.finalLowKDensity / max(S0.finalLowKDensity, eps), ...
        S.meanLowKDensity / max(S0.meanLowKDensity, eps), ...
        S.finalStdN / max(S0.finalStdN, eps), ...
        S.finalPkinMean / max(S0.finalPkinMean, eps), ...
        out.compressibility.KeffPkin / max(classicOut.compressibility.KeffPkin, eps), ...
        S.finalRhoPhysicalRelError, S.finalKBTCell, S.finalMassFluxResidual};
end
if ~isempty(rows)
    C = cell2table(rows, 'VariableNames', {'method','lowKFinalRatio','lowKMeanRatio','stdNFinalRatio', ...
        'PkinFinalRatio','KeffPkinRatio','rhoRelErrFinal','kBTFinal','massFluxResidualFinal'});
end
end

function opts = parse_inputs(varargin)
opts = struct();
if nargin == 1 && isstruct(varargin{1})
    opts = varargin{1};
    return;
end
if mod(nargin, 2) ~= 0
    error('Use name-value pairs or a single struct.');
end
for k = 1:2:nargin
    opts.(char(varargin{k})) = varargin{k+1};
end
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
