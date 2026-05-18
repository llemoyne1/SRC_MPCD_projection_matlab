function results = run_q9_tg32_operator_comparison(varargin)
%RUN_Q9_TG32_OPERATOR_COMPARISON Compare Q9 mass-flux operators on forced TG.
%
% This script runs the same small periodic Taylor--Green case with Q9
% mass-flux projector variants:
%   1. historical FFT periodic projector,
%   2. finite-volume/general-BC elliptic projector with elliptic low-pass,
%   3. optional finite-volume/general-BC projector with exact periodic FFT low-k target mask.
%
% The goal is not to replace the FFT production path immediately, but to
% check whether the future unified elliptic operator gives the same bulk
% behaviour on a purely periodic problem.
%
% Example:
%   results = run_q9_tg32_operator_comparison('nSteps', 5000);
%   results = run_q9_tg32_operator_comparison('nSteps', 10000, 'visualEnable', true);

opts = parse_options(varargin{:});
outDir = fullfile(pwd, sprintf('tg32_q9_operator_comparison_%s', datestr(now,'yyyymmdd_HHMMSS')));
if ~exist(outDir, 'dir'), mkdir(outDir); end

logFile = fullfile(outDir, 'console_log.txt');
diary(logFile);
cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('=== Q9 TG32 mass-flux operator comparison ===\n');
fprintf('outputDir   : %s\n', outDir);
fprintf('grid/gamma  : %d x %d / %g\n', opts.Nx, opts.Ny, opts.gamma);
fprintf('nSteps      : %d\n', opts.nSteps);
fprintf('force       : %.12g\n', opts.taylorGreenForceAmplitude);
fprintf('visual      : %d every %d\n\n', opts.visualEnable, opts.visualEvery);

cases = make_cases(opts);
outputs = struct();
rows = cell(numel(cases), 1);

for ic = 1:numel(cases)
    c = cases(ic);
    fprintf('\n######################################################################\n');
    fprintf('Case %d/%d: %s\n', ic, numel(cases), c.label);
    fprintf('operator=%s, targetFilter=%s\n', c.params.massFluxProjectionOperator, c.params.massFluxTargetFilter);
    fprintf('######################################################################\n\n');
    try
        out = run_projection_forced_taylor_green_demo(c.params);
        outputs.(c.fieldName) = out;
        save(fullfile(outDir, [c.label '.mat']), 'out', '-v7.3');
        rows{ic} = summarize_case(c.label, c.params, out, '');
    catch ME
        warning('run_q9_tg32_operator_comparison:CaseFailed', '%s failed: %s', c.label, ME.message);
        outputs.(c.fieldName) = [];
        rows{ic} = summarize_case(c.label, c.params, [], getReport(ME, 'basic', 'hyperlinks', 'off'));
    end
end

summary = vertcat(rows{:});
writetable(summary, fullfile(outDir, 'summary.csv'));
save(fullfile(outDir, 'comparison_results.mat'), 'opts', 'cases', 'outputs', 'summary', '-v7.3');

fprintf('\n=== Operator comparison complete ===\n');
fprintf('Summary: %s\n', fullfile(outDir, 'summary.csv'));
disp(summary);

results = struct();
results.opts = opts;
results.cases = cases;
results.outputs = outputs;
results.summary = summary;
results.outputDir = outDir;
results.logFile = logFile;
end

function opts = parse_options(varargin)
opts = struct();
opts.Nx = 32;
opts.Ny = 32;
opts.gamma = 20;
opts.seed = 11;
opts.Lx = 1.0;
opts.Ly = 1.0;
opts.dt = 1.0e-3;
opts.kBT = 0.01;
opts.alphaDeg = 90;
opts.nSteps = 5000;
opts.sampleEvery = 50;
opts.progressEvery = 500;
opts.taylorGreenForceAmplitude = 0.12;
opts.taylorGreenInitialAmplitude = 0.10;
opts.massFluxDensityRelaxationBeta = 5.0e-4;
opts.massFluxLowKMaxIndex = 2;
opts.massFluxFinalVelocityProjectionCleanup = true;
opts.massFluxFinalVelocityProjectionStrength = 0.5;
opts.visualEnable = true;
opts.visualEvery = 500;
opts.visualPause = 0.001;
opts.visualMaxParticles = 4000;
opts.meanFieldEnable = true;
opts.meanFieldVisualEnable = true;
opts.taylorGreenViscosityFitFraction = 0.60;
opts.maxWallClockSeconds = Inf;
opts.abortOnUnstable = true;
opts.generalEllipticSolver = 'backslash';
opts.generalEllipticUseSolveCache = true;
opts.generalEllipticFactorization = 'auto';
opts.generalEllipticTargetFilter = 'elliptic_lowpass';
opts.includeGeneralFFTLowKCase = true;

if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(string(varargin{k}));
    if ~isfield(opts, name)
        error('Unknown option: %s', name);
    end
    opts.(name) = varargin{k+1};
end
end

function cases = make_cases(opts)
base = struct();
base.method = 'q9';
base.Lx = opts.Lx;
base.Ly = opts.Ly;
base.Nx = opts.Nx;
base.Ny = opts.Ny;
base.gamma = opts.gamma;
base.seed = opts.seed;
base.dt = opts.dt;
base.kBT = opts.kBT;
base.alphaDeg = opts.alphaDeg;
base.initialPopulationMode = 'exact_per_cell';
base.initialVelocityZeroGlobalMean = true;
base.useRandomGridShift = true;
base.bodyForceX = 0;
base.bodyForceY = 0;
base.taylorGreenInitialAmplitude = opts.taylorGreenInitialAmplitude;
base.taylorGreenAmplitude = opts.taylorGreenInitialAmplitude;
base.taylorGreenThermalNoise = true;
base.taylorGreenForceEnable = true;
base.taylorGreenForceAmplitude = opts.taylorGreenForceAmplitude;
base.taylorGreenForceZeroMeanKick = true;
base.taylorGreenModeX = 1;
base.taylorGreenModeY = 1;
base.nSteps = opts.nSteps;
base.sampleEvery = opts.sampleEvery;
base.progressEvery = opts.progressEvery;
base.visualEnable = logical(opts.visualEnable);
base.visualEvery = opts.visualEvery;
base.visualPause = opts.visualPause;
base.visualMaxParticles = opts.visualMaxParticles;
base.meanFieldEnable = logical(opts.meanFieldEnable);
base.meanFieldVisualEnable = logical(opts.meanFieldVisualEnable);
base.meanFieldStartStep = round(0.2*opts.nSteps);
base.thermostatAfterStep = true;
base.thermostatAfterProjection = true;
base.projectionEnable = true;
base.projectionStrength = 1.0;
base.projectionInterpolationMethod = 'nearest';
base.projectionMomentumCorrectionEnable = true;
base.projectionMomentumCorrectionMode = 'particle_global_exact';
base.massFluxProjectionMode = 'relax_to_uniform_lowk';
base.massFluxProjectionStrength = 1.0;
base.massFluxDensityRelaxationBeta = opts.massFluxDensityRelaxationBeta;
base.massFluxApplyAfterVelocityProjection = true;
base.massFluxLowKMaxIndex = opts.massFluxLowKMaxIndex;
base.lowKMaxIndex = opts.massFluxLowKMaxIndex;
base.massFluxFinalVelocityProjectionCleanup = logical(opts.massFluxFinalVelocityProjectionCleanup);
base.massFluxFinalVelocityProjectionStrength = opts.massFluxFinalVelocityProjectionStrength;
base.taylorGreenViscosityFitFraction = opts.taylorGreenViscosityFitFraction;
base.abortOnUnstable = logical(opts.abortOnUnstable);
base.maxWallClockSeconds = opts.maxWallClockSeconds;

pFFT = base;
pFFT.massFluxProjectionOperator = 'periodic_fft';
pFFT.massFluxTargetFilter = 'lowpass_fft';
pFFT.visualFigureId = 510;
pFFT.meanFieldFigureId = 530;
pFFT.visualTitleSuffix = 'Q9 FFT periodic';

pGEN = base;
pGEN.massFluxProjectionOperator = 'general_bc';
pGEN.massFluxTargetFilter = opts.generalEllipticTargetFilter;
pGEN.massFluxEllipticSolver = opts.generalEllipticSolver;
pGEN.massFluxEllipticAlphaMode = 'constant';
pGEN.massFluxEllipticUseSolveCache = logical(opts.generalEllipticUseSolveCache);
pGEN.massFluxEllipticFactorization = opts.generalEllipticFactorization;
pGEN.massFluxBCLeft = 'periodic';
pGEN.massFluxBCRight = 'periodic';
pGEN.massFluxBCBottom = 'periodic';
pGEN.massFluxBCTop = 'periodic';
pGEN.visualFigureId = 511;
pGEN.meanFieldFigureId = 531;
pGEN.visualTitleSuffix = 'Q9 generalBC elliptic lowpass';

pGENFFT = pGEN;
pGENFFT.massFluxTargetFilter = 'periodic_fft_lowk';
pGENFFT.visualFigureId = 512;
pGENFFT.meanFieldFigureId = 532;
pGENFFT.visualTitleSuffix = 'Q9 generalBC FFT low-k target';

if logical(opts.includeGeneralFFTLowKCase)
    cases = repmat(struct('label','','fieldName','','params',struct()), 3, 1);
else
    cases = repmat(struct('label','','fieldName','','params',struct()), 2, 1);
end
cases(1).label = sprintf('tg_q9_fft_N%dx%d', opts.Nx, opts.Ny);
cases(1).fieldName = 'fft';
cases(1).params = pFFT;
cases(2).label = sprintf('tg_q9_generalbc_elliptic_N%dx%d', opts.Nx, opts.Ny);
cases(2).fieldName = 'general_bc_elliptic';
cases(2).params = pGEN;
if numel(cases) >= 3
    cases(3).label = sprintf('tg_q9_generalbc_fftlowk_N%dx%d', opts.Nx, opts.Ny);
    cases(3).fieldName = 'general_bc_fftlowk';
    cases(3).params = pGENFFT;
end
end

function row = summarize_case(label, params, out, errorMessage)
if isempty(out)
    row = table({label}, {char(string(params.massFluxProjectionOperator))}, false, params.Nx, params.Ny, params.gamma, params.nSteps, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, {errorMessage}, ...
        'VariableNames', {'label','operator','ok','Nx','Ny','gamma','nSteps','actualStep','elapsed','nuEffPreferred','nuEffDerivative','nuEffPlateau','finalLowKDensity','finalKBT','finalAmplitude','finalCoherence','finalHighKFraction','errorMessage'});
    return;
end
s = out.summary;
row = table({label}, {char(string(params.massFluxProjectionOperator))}, true, params.Nx, params.Ny, params.gamma, params.nSteps, out.actualLastStep, out.elapsedWallClock, ...
    s.nuEffFitPreferred, s.nuEffFitDerivativeKnownF, s.nuEffFitPlateau, s.finalLowKDensity, s.finalKBTCell, s.finalAmplitude, s.finalCoherence, s.finalHighKFraction, {errorMessage}, ...
    'VariableNames', {'label','operator','ok','Nx','Ny','gamma','nSteps','actualStep','elapsed','nuEffPreferred','nuEffDerivative','nuEffPlateau','finalLowKDensity','finalKBT','finalAmplitude','finalCoherence','finalHighKFraction','errorMessage'});
end
