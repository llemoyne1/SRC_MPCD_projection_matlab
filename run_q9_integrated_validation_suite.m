function suite = run_q9_integrated_validation_suite(varargin)
%RUN_Q9_INTEGRATED_VALIDATION_SUITE Autonomous Q6/Q9 validation campaign.
%
%   suite = run_q9_integrated_validation_suite(...)
%
% Runs the three current non-regression validations from a single clean code
% snapshot:
%   1) periodic forced Taylor--Green, usually Q9 only;
%   2) Poiseuille channel with wall virtual particles;
%   3) piston smooth-ramp compression with virial EOS and rolling Pwall-Pbulk.
%
% The runner is designed for long unattended runs.  Each case is executed in
% its own subfolder, writes DONE.txt only after all case outputs have been
% saved, and updates suite_state.mat after every case.  If MATLAB or the job is
% interrupted, rerun the same command with 'resume',true: cases with DONE.txt
% are skipped and the campaign continues from the first incomplete case.
%
% Example long run:
%   suite = run_q9_integrated_validation_suite( ...
%       'preset','long15h', ...
%       'visualEnable',false, ...
%       'resume',true);
%
% Example smoke test:
%   suite = run_q9_integrated_validation_suite( ...
%       'preset','smoke', ...
%       'visualEnable',false, ...
%       'resume',false);
%
% The default long15h settings are intentionally conservative.  They can be
% overridden with name-value pairs: tgSteps, poiseuilleSteps,
% pistonInitialHoldSteps, pistonRampSteps, pistonFinalHoldSteps,
% pistonWindowSteps, outputRoot, etc.

opts = parse_suite_inputs(varargin{:});
opts = apply_suite_preset(opts);

if isempty(opts.outputRoot)
    opts.outputRoot = fullfile(pwd, ['q9_integrated_validation_suite_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if ~exist(opts.outputRoot, 'dir')
    mkdir(opts.outputRoot);
end

logFile = fullfile(opts.outputRoot, 'suite_console_log.txt');
try
    diary(logFile);
    diary on;
catch ME
    warning('run_q9_integrated_validation_suite:DiaryFailed', ...
        'Could not start suite diary: %s', ME.message);
end
cleanupObj = onCleanup(@() local_diary_off()); %#ok<NASGU>

tSuite = tic;
fprintf('=== Q9 integrated validation suite ===\n');
fprintf('preset      : %s\n', opts.preset);
fprintf('outputRoot  : %s\n', opts.outputRoot);
fprintf('resume      : %d\n', logical(opts.resume));
fprintf('visualEnable: %d\n\n', logical(opts.visualEnable));

write_manifest(opts);

suite = struct();
suite.opts = opts;
suite.outputRoot = opts.outputRoot;
suite.startedAt = datestr(now, 31);
suite.caseResults = struct();
suite.caseSummaries = table();

caseList = {'TaylorGreen','Poiseuille','PistonRamp'};
for icase = 1:numel(caseList)
    caseName = caseList{icase};
    caseDir = fullfile(opts.outputRoot, caseName);
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    doneFile = fullfile(caseDir, 'DONE.txt');
    errorFile = fullfile(caseDir, 'ERROR.txt');

    if logical(opts.resume) && exist(doneFile, 'file')
        fprintf('\n=== %s already DONE: skipping ===\n', caseName);
        row = load_case_summary_if_available(caseDir, caseName);
        suite.caseSummaries = append_summary_row(suite.caseSummaries, row);
        suite.caseResults.(caseName) = struct('status','skipped_done','caseDir',caseDir);
        save_suite_state(opts.outputRoot, suite, caseName, 'skipped_done');
        continue;
    end

    if exist(errorFile, 'file') && ~logical(opts.resume)
        delete(errorFile);
    end

    fprintf('\n=== Running case %d/%d: %s ===\n', icase, numel(caseList), caseName);
    tCase = tic;
    try
        switch caseName
            case 'TaylorGreen'
                [result, row] = run_case_taylor_green(opts, caseDir);
            case 'Poiseuille'
                [result, row] = run_case_poiseuille(opts, caseDir);
            case 'PistonRamp'
                [result, row] = run_case_piston_ramp(opts, caseDir);
            otherwise
                error('Unknown case: %s', caseName);
        end
        row.elapsedCaseWallClock = toc(tCase);
        row.status = string('done');
        writetable(row, fullfile(caseDir, 'case_summary.csv'));
        write_done_file(doneFile, caseName, row);
        suite.caseResults.(caseName) = result;
        suite.caseSummaries = append_summary_row(suite.caseSummaries, row);
        save_suite_state(opts.outputRoot, suite, caseName, 'done');
        fprintf('=== %s DONE in %.2f s ===\n', caseName, toc(tCase));
    catch ME
        row = make_error_summary_row(caseName, caseDir, ME, toc(tCase));
        suite.caseResults.(caseName) = struct('status','error','caseDir',caseDir,'error',ME);
        suite.caseSummaries = append_summary_row(suite.caseSummaries, row);
        write_error_file(errorFile, caseName, ME);
        save_suite_state(opts.outputRoot, suite, caseName, 'error');
        warning('run_q9_integrated_validation_suite:CaseFailed', ...
            '%s failed: %s', caseName, ME.message);
        if logical(opts.stopOnError)
            rethrow(ME);
        end
    end
end

suite.finishedAt = datestr(now, 31);
suite.elapsedWallClock = toc(tSuite);

summaryFile = fullfile(opts.outputRoot, 'integrated_summary.csv');
writetable(suite.caseSummaries, summaryFile);
write_integrated_report(opts, suite);
try
    save(fullfile(opts.outputRoot, 'integrated_suite_output.mat'), 'suite', '-v7.3');
catch ME
    warning('run_q9_integrated_validation_suite:SaveFailed', ...
        'Could not save integrated_suite_output.mat: %s', ME.message);
end

fprintf('\n=== Integrated suite finished in %.2f s ===\n', suite.elapsedWallClock);
fprintf('Summary: %s\n', summaryFile);
end

%% Case runners ----------------------------------------------------------
function [result, row] = run_case_taylor_green(opts, caseDir)
methods = cellstr(opts.tgMethods);
outs = struct();
rows = table();

for im = 1:numel(methods)
    method = lower(strrep(methods{im}, '-', '_'));
    p = make_tg_params(opts, method, caseDir, im);
    fprintf('\n--- Taylor-Green %s ---\n', upper(method));
    out = run_projection_forced_taylor_green_demo(p);
    outs.(method) = out;
    tsFile = fullfile(caseDir, sprintf('tg_%s_timeseries.csv', method));
    write_tg_timeseries_csv(tsFile, out);
    one = summarize_tg_case(method, out, caseDir);
    rows = [rows; one]; %#ok<AGROW>
    try
        save(fullfile(caseDir, sprintf('tg_%s_output.mat', method)), 'out', 'p', '-v7.3');
    catch ME
        warning('run_q9_integrated_validation_suite:TgSaveFailed', ...
            'Could not save TG %s output: %s', method, ME.message);
    end
end
writetable(rows, fullfile(caseDir, 'tg_summary_by_method.csv'));
result = struct('status','done','caseDir',caseDir,'methods',{methods},'summary',rows);

% Integrated row: prefer Q9, otherwise last method.
pick = preferred_method_index(rows, 'q9');
row = rows(pick, :);
row.caseName = string('TaylorGreen');
row.caseDir = string(caseDir);
row.metricPrimary = row.finalLowKDensity;
row.metricSecondary = row.nuEffFitPreferred;
end

function [result, row] = run_case_poiseuille(opts, caseDir)
methods = cellstr(opts.poiseuilleMethods);
fprintf('\n--- Poiseuille wallVP methods: %s ---\n', strjoin(methods, ', '));

results = run_q9_poiseuille_wall_virtual_particles_validation( ...
    'methods', methods, ...
    'Nx', opts.poiseuilleNx, ...
    'Ny', opts.poiseuilleNy, ...
    'gamma', opts.gamma, ...
    'nSteps', opts.poiseuilleSteps, ...
    'sampleEvery', opts.sampleEvery, ...
    'progressEvery', opts.progressEvery, ...
    'visualEnable', opts.visualEnable, ...
    'visualEvery', opts.visualEvery, ...
    'bodyForceX', opts.poiseuilleBodyForceX, ...
    'wallModeY', opts.poiseuilleWallModeY, ...
    'geometryMode', opts.wallVPGeometryMode, ...
    'densityFactor', opts.wallVPDensityFactor, ...
    'nLayers', opts.wallVPLayers, ...
    'thermal', opts.wallVPThermal, ...
    'stochasticCount', opts.wallVPStochasticCount, ...
    'forceRandomShiftY', opts.wallVPForceRandomShiftY, ...
    'initialPoiseuille', opts.poiseuilleInitialProfileEnable, ...
    'nuGuess', opts.poiseuilleNuGuess, ...
    'fitWindowTime', opts.poiseuilleFitWindowTime, ...
    'accelerationWindowTime', opts.poiseuilleAccelerationWindowTime, ...
    'nuReferenceNy', opts.poiseuilleNuReferenceNy, ...
    'nuScalingEnable', opts.poiseuilleNuScalingEnable);

summary = results.summary;
summary.caseName = repmat(string('Poiseuille'), height(summary), 1);
summary.caseDir = repmat(string(caseDir), height(summary), 1);
writetable(summary, fullfile(caseDir, 'poiseuille_summary_by_method.csv'));
try
    save(fullfile(caseDir, 'poiseuille_output.mat'), 'results', '-v7.3');
catch ME
    warning('run_q9_integrated_validation_suite:PoiseuilleSaveFailed', ...
        'Could not save Poiseuille output: %s', ME.message);
end

result = struct('status','done','caseDir',caseDir,'methods',{methods},'summary',summary);
pick = preferred_method_index(summary, 'q9');
row0 = summary(pick, :);
row = normalize_poiseuille_summary_row(row0, caseDir);
end

function [result, row] = run_case_piston_ramp(opts, caseDir)
fprintf('\n--- Piston smooth ramp Q9_FULL ---\n');

out = run_q9_piston_smooth_wallbulk_validation( ...
    'outputRoot', caseDir, ...
    'method', 'q9_full', ...
    'Nx', opts.pistonNx, ...
    'Ny', opts.pistonNy, ...
    'gamma', opts.gamma, ...
    'compressionFinal', opts.pistonCompressionFinal, ...
    'initialHoldSteps', opts.pistonInitialHoldSteps, ...
    'rampSteps', opts.pistonRampSteps, ...
    'finalHoldSteps', opts.pistonFinalHoldSteps, ...
    'wallBulkResidualWindowSteps', opts.pistonWindowSteps, ...
    'wallBulkResidualLayerCells', opts.pistonLayerCells, ...
    'sampleEvery', opts.sampleEvery, ...
    'progressEvery', opts.progressEvery, ...
    'visualEnable', opts.visualEnable, ...
    'visualEvery', opts.visualEvery, ...
    'plotEnable', opts.plotEnable, ...
    'saveOutput', true, ...
    'densityFactor', opts.wallVPDensityFactor, ...
    'Kvirial', opts.Kvirial, ...
    'virialBeta', opts.virialBeta, ...
    'virialMaxDuFractionThermal', opts.virialMaxDuFractionThermal, ...
    'maxWallClockSeconds', opts.pistonMaxWallClockSeconds);

T = out.diagTable;
strictSummary = recompute_piston_postramp_windows(T, opts);
writetable(strictSummary, fullfile(caseDir, 'piston_postramp_window_summary.csv'));
try
    save(fullfile(caseDir, 'piston_ramp_output_with_strict_summary.mat'), ...
        'out', 'strictSummary', '-v7.3');
catch ME
    warning('run_q9_integrated_validation_suite:PistonSaveFailed', ...
        'Could not save piston output: %s', ME.message);
end

result = struct('status','done','caseDir',caseDir,'summary',out.wallBulk.summaryTable, ...
    'postRampWindowSummary', strictSummary);
row = summarize_piston_case(out, strictSummary, caseDir, opts);
end

%% Parameter builders ----------------------------------------------------
function p = make_tg_params(opts, method, caseDir, methodIndex)
p = struct();
p.method = method;
p.Lx = 1.0; p.Ly = 1.0;
p.Nx = opts.tgNx; p.Ny = opts.tgNy; p.gamma = opts.gamma; p.seed = opts.seed + 100*methodIndex;
p.dt = opts.dt; p.kBT = opts.kBT; p.alphaDeg = 90;
p.initialPopulationMode = 'exact_per_cell';
p.initialVelocityZeroGlobalMean = true;
p.taylorGreenThermalNoise = true;
p.taylorGreenInitialAmplitude = 0.10;
p.taylorGreenAmplitude = p.taylorGreenInitialAmplitude;
p.taylorGreenForceEnable = true;
p.taylorGreenForceAmplitude = opts.tgForceAmplitude;
p.taylorGreenForceZeroMeanKick = true;
p.taylorGreenModeX = 1; p.taylorGreenModeY = 1;
p.bodyForceX = 0; p.bodyForceY = 0;
p.useRandomGridShift = true;
p.projectionInterpolationMethod = 'nearest';
p.thermostatAfterStep = true;
p.thermostatAfterProjection = true;
p.computeFullDiagnosticsEveryStep = false;
p.abortOnUnstable = true;
p.maxStableKBTCell = 0.2;
p.maxStableHighKFraction = 0.9;
p.maxStableDensityRelRMS = 1.2;
p.maxWallClockSeconds = opts.tgMaxWallClockSeconds;
p.nSteps = opts.tgSteps;
p.sampleEvery = opts.sampleEvery;
p.progressEvery = opts.progressEvery;
p.visualEnable = opts.visualEnable;
p.visualEvery = opts.visualEvery;
p.visualSaveFrames = opts.visualSaveFrames;
p.visualFigureId = 410 + methodIndex;
p.visualFigureName = sprintf('Integrated TG %s', upper(method));
p.visualFrameDir = fullfile(caseDir, 'frames', method);
p.visualFramePrefix = sprintf('tg_%s', method);
p.visualPause = 0.001;
p.visualMaxParticles = 8000;
p.visualDensityCLim = 0.6;
p.visualOmegaCLim = NaN;
p.visualQuiverStrideX = max(2, round(p.Nx/18));
p.visualQuiverStrideY = max(2, round(p.Ny/12));
p.visualQuiverScale = 1.2;
p.meanFieldEnable = true;
p.meanFieldVisualEnable = opts.visualEnable;
p.meanFieldStartStep = round(0.2*p.nSteps);
p.meanFieldFigureId = 430 + methodIndex;
p.meanFieldFigureName = sprintf('Integrated TG running mean %s', upper(method));
p.meanFieldSaveFrames = false;
p.meanFieldFrameDir = fullfile(caseDir, 'mean_frames', method);
p.meanFieldFramePrefix = sprintf('tg_mean_%s', method);
p.taylorGreenViscosityFitFraction = 0.60;
p.lowKMaxIndex = 2;
p.massFluxProjectionMode = 'relax_to_uniform_lowk';
p.massFluxProjectionStrength = 1.0;
p.massFluxDensityRelaxationBeta = opts.q9Beta;
p.massFluxApplyAfterVelocityProjection = true;
p.massFluxLowKMaxIndex = 2;
p.massFluxFinalVelocityProjectionCleanup = true;
p.massFluxFinalVelocityProjectionStrength = opts.q9CleanupStrength;
p.projectionStrength = 1.0;
p.projectionEnable = true;
p.projectionMomentumCorrectionEnable = true;
p.projectionMomentumCorrectionMode = 'particle_global_exact';

switch lower(method)
    case 'classic'
        p.projectionEnable = false;
        p.projectionStrength = 0;
        p.massFluxProjectionMode = 'off';
        p.massFluxProjectionStrength = 0;
        p.massFluxDensityRelaxationBeta = 0;
        p.massFluxApplyAfterVelocityProjection = false;
        p.massFluxFinalVelocityProjectionCleanup = false;
        p.massFluxFinalVelocityProjectionStrength = 0;
        p.projectionMomentumCorrectionEnable = false;
    case 'q6'
        p.projectionEnable = true;
        p.projectionStrength = 1;
        p.massFluxProjectionMode = 'off';
        p.massFluxProjectionStrength = 0;
        p.massFluxDensityRelaxationBeta = 0;
        p.massFluxApplyAfterVelocityProjection = false;
        p.massFluxFinalVelocityProjectionCleanup = false;
        p.massFluxFinalVelocityProjectionStrength = 0;
    case 'q9'
        % defaults above
    otherwise
        error('Unknown TG method: %s', method);
end
end

%% Summary helpers -------------------------------------------------------
function row = summarize_tg_case(method, out, caseDir)
S = out.summary;
row = base_summary_row('TaylorGreen', caseDir, method);
row.actualStep = out.actualLastStep;
row.elapsedRunWallClock = out.elapsedWallClock;
row.stoppedEarly = logical(out.stoppedEarly);
row.stopReason = string(out.stopReason);
row.finalKBT = get_field(S, 'finalKBTCell', NaN);
row.meanKBT = get_field(S, 'meanKBTCell', NaN);
row.finalLowKDensity = get_field(S, 'finalLowKDensity', NaN);
row.meanLowKDensity = get_field(S, 'meanLowKDensity', NaN);
row.nuEffFitPreferred = get_field(S, 'nuEffFitPreferred', NaN);
row.finalAmplitude = get_field(S, 'finalAmplitude', NaN);
row.meanAmplitude = get_field(S, 'meanAmplitude', NaN);
row.finalCoherence = get_field(S, 'finalCoherence', NaN);
row.meanCoherence = get_field(S, 'meanCoherence', NaN);
row.notes = string('TG periodic forced vortex');
end

function row = normalize_poiseuille_summary_row(row0, caseDir)
row = base_summary_row('Poiseuille', caseDir, string(row0.method{1}));
row.actualStep = row0.actualStep;
row.elapsedRunWallClock = row0.elapsed;
row.finalKBT = row0.finalKBT;
row.finalLowKDensity = row0.finalLowK;
row.nuEffRaw = row0.nuEffRaw;
row.nuEffScaledNyRef = row0.nuEffScaledNyRef;
row.nuReferenceNy = row0.nuReferenceNy;
row.poiseuilleR2 = row0.R2;
row.poiseuilleSNR = row0.SNR;
row.centerMinusWall = row0.finalCenterWall;
row.ratioAccel = row0.ratioAccel;
row.metricPrimary = row0.nuEffScaledNyRef;
row.metricSecondary = row0.R2;
row.notes = string('Poiseuille wallVP channel');
end

function row = summarize_piston_case(out, strictSummary, caseDir, opts)
row = base_summary_row('PistonRamp', caseDir, 'q9_full');
row.actualStep = get_field(out, 'actualLastStep', get_field(out, 'actualStep', NaN));
row.elapsedRunWallClock = out.elapsedWallClock;
row.finalKBT = tail_value(out.diagTable, 'kBTCell');
row.meanKBT = mean(column_or_nan(out.diagTable, 'kBTCell'), 'omitnan');
row.finalLowKDensity = tail_value(out.diagTable, 'lowKDensity');
row.meanLowKDensity = mean(column_or_nan(out.diagTable, 'lowKDensity'), 'omitnan');
row.finalPkin = tail_value(out.diagTable, 'PkinMean');
row.finalPvir = tail_value(out.diagTable, 'PvirMean');
row.finalPtot = tail_value(out.diagTable, 'PtotMean');
row.finalCompression = tail_value(out.diagTable, 'compression');
row.virialLimitedFraction = tail_value(out.diagTable, 'virialLimitedCellFraction');
row.metricPrimary = NaN;
row.metricSecondary = NaN;
row.notes = string('Piston smooth ramp wall/bulk rolling residual');

% Prefer the largest strict post-ramp window with at least 5 samples.
if ~isempty(strictSummary) && height(strictSummary) > 0
    idx = find(strictSummary.n >= 5 & isfinite(strictSummary.resPtotSymMean));
    if isempty(idx)
        idx = find(isfinite(strictSummary.resPtotSymMean));
    end
    if ~isempty(idx)
        k = idx(end);
        row.pistonWindowSteps = strictSummary.windowSteps(k);
        row.pistonPwallSym = strictSummary.PwallSymMean(k);
        row.pistonPtotLayerSym = strictSummary.PtotLayerSymMean(k);
        row.pistonResPtotSym = strictSummary.resPtotSymMean(k);
        row.pistonResPtotSymStd = strictSummary.resPtotSymStd(k);
        row.pistonResPtotAnti = strictSummary.resPtotAntiMean(k);
        row.metricPrimary = strictSummary.resPtotSymMean(k);
        row.metricSecondary = strictSummary.resPtotSymStd(k);
    end
end
row.pistonRampEndStep = opts.pistonInitialHoldSteps + opts.pistonRampSteps;
end

function row = base_summary_row(caseName, caseDir, method)
row = table();
row.caseName = string(caseName);
row.method = string(method);
row.status = string('pending');
row.caseDir = string(caseDir);
row.actualStep = NaN;
row.elapsedRunWallClock = NaN;
row.elapsedCaseWallClock = NaN;
row.stoppedEarly = false;
row.stopReason = string('');
row.finalKBT = NaN;
row.meanKBT = NaN;
row.finalLowKDensity = NaN;
row.meanLowKDensity = NaN;
row.nuEffFitPreferred = NaN;
row.finalAmplitude = NaN;
row.meanAmplitude = NaN;
row.finalCoherence = NaN;
row.meanCoherence = NaN;
row.nuEffRaw = NaN;
row.nuEffScaledNyRef = NaN;
row.nuReferenceNy = NaN;
row.poiseuilleR2 = NaN;
row.poiseuilleSNR = NaN;
row.centerMinusWall = NaN;
row.ratioAccel = NaN;
row.finalPkin = NaN;
row.finalPvir = NaN;
row.finalPtot = NaN;
row.finalCompression = NaN;
row.virialLimitedFraction = NaN;
row.pistonWindowSteps = NaN;
row.pistonPwallSym = NaN;
row.pistonPtotLayerSym = NaN;
row.pistonResPtotSym = NaN;
row.pistonResPtotSymStd = NaN;
row.pistonResPtotAnti = NaN;
row.pistonRampEndStep = NaN;
row.metricPrimary = NaN;
row.metricSecondary = NaN;
row.notes = string('');
end

%% Piston post-ramp recomputation ---------------------------------------
function S = recompute_piston_postramp_windows(T, opts)
if isempty(T) || height(T) == 0
    S = table();
    return;
end
windowStepsList = unique([opts.pistonWindowSteps, opts.pistonPostRampWindowStepsList(:).']);
windowStepsList = windowStepsList(isfinite(windowStepsList) & windowStepsList > 0);
rampEndStep = opts.pistonInitialHoldSteps + opts.pistonRampSteps;

if ismember('wallImpulseTimeCumulative', T.Properties.VariableNames)
    tcum = T.wallImpulseTimeCumulative;
else
    tcum = T.t;
end
step = T.step;
ItopNorm = column_or_nan(T, 'pressureTopWallTotalCumulativeMean') .* tcum;
IbotNorm = column_or_nan(T, 'pressureBottomWallTotalCumulativeMean') .* tcum;
PkinTop = column_or_nan(T, 'PkinTopLayerMean');
PkinBot = column_or_nan(T, 'PkinBottomLayerMean');
PtotTop = column_or_nan(T, 'PtotTopLayerMean');
PtotBot = column_or_nan(T, 'PtotBottomLayerMean');
if all(~isfinite(PtotTop)) && ismember('PvirTopLayerMean', T.Properties.VariableNames)
    PtotTop = PkinTop + T.PvirTopLayerMean;
end
if all(~isfinite(PtotBot)) && ismember('PvirBottomLayerMean', T.Properties.VariableNames)
    PtotBot = PkinBot + T.PvirBottomLayerMean;
end

S = table();
for w = windowStepsList
    PtopWall = nan(height(T),1);
    PbotWall = nan(height(T),1);
    PkinTopWin = nan(height(T),1); PkinBotWin = nan(height(T),1);
    PtotTopWin = nan(height(T),1); PtotBotWin = nan(height(T),1);
    for k = 1:height(T)
        k0 = find(step <= step(k) - w, 1, 'last');
        if isempty(k0), continue; end
        dtw = tcum(k) - tcum(k0);
        if dtw <= 0, continue; end
        PtopWall(k) = (ItopNorm(k) - ItopNorm(k0)) / dtw;
        PbotWall(k) = (IbotNorm(k) - IbotNorm(k0)) / dtw;
        idx = k0:k;
        tw = T.t(idx);
        if numel(idx) >= 2 && tw(end) > tw(1)
            PkinTopWin(k) = trapz(tw, PkinTop(idx)) / (tw(end)-tw(1));
            PkinBotWin(k) = trapz(tw, PkinBot(idx)) / (tw(end)-tw(1));
            PtotTopWin(k) = trapz(tw, PtotTop(idx)) / (tw(end)-tw(1));
            PtotBotWin(k) = trapz(tw, PtotBot(idx)) / (tw(end)-tw(1));
        else
            PkinTopWin(k) = PkinTop(k);
            PkinBotWin(k) = PkinBot(k);
            PtotTopWin(k) = PtotTop(k);
            PtotBotWin(k) = PtotBot(k);
        end
    end
    PwallSym = 0.5*(PtopWall + PbotWall);
    PkinLayerSym = 0.5*(PkinTopWin + PkinBotWin);
    PtotLayerSym = 0.5*(PtotTopWin + PtotBotWin);
    resPkinSym = PwallSym - PkinLayerSym;
    resPtotSym = PwallSym - PtotLayerSym;
    resPtotTop = PtopWall - PtotTopWin;
    resPtotBot = PbotWall - PtotBotWin;
    resPtotAnti = 0.5*(resPtotTop - resPtotBot);

    idxStrict = column_or_nan(T, 'compression') > 0.999*opts.pistonCompressionFinal & ...
        step - w >= rampEndStep & isfinite(resPtotSym);

    row = table();
    row.windowSteps = w;
    row.n = nnz(idxStrict);
    row.rampEndStep = rampEndStep;
    row.PwallSymMean = mean(PwallSym(idxStrict), 'omitnan');
    row.PwallSymStd = std(PwallSym(idxStrict), 'omitnan');
    row.PkinLayerSymMean = mean(PkinLayerSym(idxStrict), 'omitnan');
    row.PtotLayerSymMean = mean(PtotLayerSym(idxStrict), 'omitnan');
    row.resPkinSymMean = mean(resPkinSym(idxStrict), 'omitnan');
    row.resPkinSymStd = std(resPkinSym(idxStrict), 'omitnan');
    row.resPtotSymMean = mean(resPtotSym(idxStrict), 'omitnan');
    row.resPtotSymStd = std(resPtotSym(idxStrict), 'omitnan');
    row.resPtotAntiMean = mean(resPtotAnti(idxStrict), 'omitnan');
    row.resPtotAntiStd = std(resPtotAnti(idxStrict), 'omitnan');
    S = [S; row]; %#ok<AGROW>
end
end

%% I/O helpers -----------------------------------------------------------
function write_tg_timeseries_csv(filename, out)
S = out.series;
vars = fieldnames(S);
T = table();
for i = 1:numel(vars)
    x = S.(vars{i});
    if isvector(x) && numel(x) == numel(S.step)
        T.(vars{i}) = x(:);
    end
end
writetable(T, filename);
end

function save_suite_state(outputRoot, suite, lastCaseName, lastStatus)
suiteState = struct();
suiteState.updatedAt = datestr(now, 31);
suiteState.lastCaseName = lastCaseName;
suiteState.lastStatus = lastStatus;
suiteState.caseSummaries = suite.caseSummaries;
suiteState.outputRoot = outputRoot;
try
    save(fullfile(outputRoot, 'suite_state.mat'), 'suiteState');
    writetable(suite.caseSummaries, fullfile(outputRoot, 'integrated_summary_partial.csv'));
catch ME
    warning('run_q9_integrated_validation_suite:StateSaveFailed', ...
        'Could not save suite state: %s', ME.message);
end
end

function write_manifest(opts)
file = fullfile(opts.outputRoot, 'manifest.txt');
fid = fopen(file, 'w');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Q9 integrated validation suite\n');
fprintf(fid, 'createdAt: %s\n', datestr(now, 31));
fprintf(fid, 'outputRoot: %s\n', opts.outputRoot);
fprintf(fid, 'preset: %s\n', opts.preset);
fprintf(fid, 'MATLAB version: %s\n', version);
try
    [status, hash] = system('git rev-parse HEAD');
    if status == 0, fprintf(fid, 'git HEAD: %s\n', strtrim(hash)); end
    [status, branch] = system('git branch --show-current');
    if status == 0, fprintf(fid, 'git branch: %s\n', strtrim(branch)); end
    [status, stat] = system('git status --short');
    if status == 0
        fprintf(fid, 'git status --short:\n%s\n', stat);
    end
catch
end
fprintf(fid, '\nOptions:\n');
fn = fieldnames(opts);
for i = 1:numel(fn)
    val = opts.(fn{i});
    fprintf(fid, '%s: %s\n', fn{i}, value_to_string(val));
end
end

function write_done_file(doneFile, caseName, row)
fid = fopen(doneFile, 'w');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s DONE\n', caseName);
fprintf(fid, 'time: %s\n', datestr(now, 31));
if istable(row)
    fprintf(fid, 'summary:\n');
    disp(row);
end
end

function write_error_file(errorFile, caseName, ME)
fid = fopen(errorFile, 'w');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s ERROR\n', caseName);
fprintf(fid, 'time: %s\n', datestr(now, 31));
fprintf(fid, 'message: %s\n', ME.message);
fprintf(fid, 'identifier: %s\n', ME.identifier);
fprintf(fid, 'stack:\n');
for k = 1:numel(ME.stack)
    fprintf(fid, '  %s:%d in %s\n', ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
end
end

function write_integrated_report(opts, suite)
file = fullfile(opts.outputRoot, 'integrated_report.txt');
fid = fopen(file, 'w');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Q9 integrated validation suite report\n');
fprintf(fid, 'started : %s\n', suite.startedAt);
fprintf(fid, 'finished: %s\n', suite.finishedAt);
fprintf(fid, 'elapsed : %.2f s\n', suite.elapsedWallClock);
fprintf(fid, 'preset  : %s\n', opts.preset);
fprintf(fid, 'root    : %s\n\n', opts.outputRoot);
fprintf(fid, 'Case summaries are in integrated_summary.csv.\n');
fprintf(fid, 'Individual case outputs are in TaylorGreen/, Poiseuille/, and PistonRamp/.\n\n');
if ~isempty(suite.caseSummaries)
    for i = 1:height(suite.caseSummaries)
        row = suite.caseSummaries(i,:);
        fprintf(fid, '%s [%s] status=%s metricPrimary=%g metricSecondary=%g\n', ...
            row.caseName, row.method, row.status, row.metricPrimary, row.metricSecondary);
    end
end
end

%% Input parsing ---------------------------------------------------------
function opts = parse_suite_inputs(varargin)
opts = struct();
opts.preset = 'long15h';
opts.outputRoot = '';
opts.resume = true;
opts.stopOnError = false;
opts.visualEnable = false;
opts.visualEvery = 1000;
opts.visualSaveFrames = false;
opts.plotEnable = true;
opts.seed = 11;
opts.dt = 1.0e-3;
opts.kBT = 0.01;
opts.gamma = 20;
opts.sampleEvery = 100;
opts.progressEvery = 1000;
opts.q9Beta = 5.0e-4;
opts.q9CleanupStrength = 0.5;
opts.tgMethods = {'q9'};
opts.tgNx = 64; opts.tgNy = 64; opts.tgSteps = NaN; opts.tgForceAmplitude = 0.08*1.5; opts.tgMaxWallClockSeconds = Inf;
opts.poiseuilleMethods = {'q9'};
opts.poiseuilleNx = 64; opts.poiseuilleNy = 64; opts.poiseuilleSteps = NaN;
opts.poiseuilleBodyForceX = 0.005; opts.poiseuilleWallModeY = 'bounceback';
opts.poiseuilleInitialProfileEnable = true; opts.poiseuilleNuGuess = 0.021;
opts.poiseuilleFitWindowTime = 5; opts.poiseuilleAccelerationWindowTime = 5;
opts.poiseuilleNuReferenceNy = 64; opts.poiseuilleNuScalingEnable = true;
opts.wallVPGeometryMode = 'shifted_solid_fraction'; opts.wallVPDensityFactor = 0.8; opts.wallVPLayers = 1;
opts.wallVPThermal = true; opts.wallVPStochasticCount = true; opts.wallVPForceRandomShiftY = true;
opts.pistonNx = 32; opts.pistonNy = 32; opts.pistonCompressionFinal = 0.01;
opts.pistonInitialHoldSteps = NaN; opts.pistonRampSteps = NaN; opts.pistonFinalHoldSteps = NaN; opts.pistonWindowSteps = NaN;
opts.pistonLayerCells = 3; opts.pistonPostRampWindowStepsList = [1500 2500 3500 5000];
opts.Kvirial = 0.01; opts.virialBeta = 0.02; opts.virialMaxDuFractionThermal = 0.025; opts.pistonMaxWallClockSeconds = Inf;

if nargin == 1 && isstruct(varargin{1})
    user = varargin{1};
    fn = fieldnames(user);
    for i = 1:numel(fn), opts.(fn{i}) = user.(fn{i}); end
    return;
end
if mod(nargin, 2) ~= 0
    error('Use name-value pairs or a single struct.');
end
for i = 1:2:nargin
    name = char(varargin{i});
    opts.(name) = varargin{i+1};
end
end

function opts = apply_suite_preset(opts)
switch lower(char(string(opts.preset)))
    case 'smoke'
        opts.tgSteps = default_if_nan(opts.tgSteps, 300);
        opts.poiseuilleSteps = default_if_nan(opts.poiseuilleSteps, 300);
        opts.pistonInitialHoldSteps = default_if_nan(opts.pistonInitialHoldSteps, 50);
        opts.pistonRampSteps = default_if_nan(opts.pistonRampSteps, 150);
        opts.pistonFinalHoldSteps = default_if_nan(opts.pistonFinalHoldSteps, 100);
        opts.pistonWindowSteps = default_if_nan(opts.pistonWindowSteps, 50);
        opts.sampleEvery = min(opts.sampleEvery, 50);
        opts.progressEvery = min(opts.progressEvery, 100);
    case {'overnight','long15h','long'}
        opts.tgSteps = default_if_nan(opts.tgSteps, 20000);
        opts.poiseuilleSteps = default_if_nan(opts.poiseuilleSteps, 25000);
        opts.pistonInitialHoldSteps = default_if_nan(opts.pistonInitialHoldSteps, 3000);
        opts.pistonRampSteps = default_if_nan(opts.pistonRampSteps, 8000);
        opts.pistonFinalHoldSteps = default_if_nan(opts.pistonFinalHoldSteps, 14000);
        opts.pistonWindowSteps = default_if_nan(opts.pistonWindowSteps, 4000);
    case 'medium'
        opts.tgSteps = default_if_nan(opts.tgSteps, 8000);
        opts.poiseuilleSteps = default_if_nan(opts.poiseuilleSteps, 10000);
        opts.pistonInitialHoldSteps = default_if_nan(opts.pistonInitialHoldSteps, 1500);
        opts.pistonRampSteps = default_if_nan(opts.pistonRampSteps, 4000);
        opts.pistonFinalHoldSteps = default_if_nan(opts.pistonFinalHoldSteps, 6000);
        opts.pistonWindowSteps = default_if_nan(opts.pistonWindowSteps, 2000);
    otherwise
        error('Unknown suite preset: %s', opts.preset);
end
end

%% Generic utilities -----------------------------------------------------
function row = make_error_summary_row(caseName, caseDir, ME, elapsed)
row = base_summary_row(caseName, caseDir, '');
row.status = string('error');
row.elapsedCaseWallClock = elapsed;
row.stopReason = string(ME.message);
row.notes = string(ME.identifier);
end

function row = load_case_summary_if_available(caseDir, caseName)
file = fullfile(caseDir, 'case_summary.csv');
if exist(file, 'file')
    row = readtable(file);
    return;
end
row = base_summary_row(caseName, caseDir, '');
row.status = string('skipped_done');
end

function T = append_summary_row(T, row)
if isempty(T)
    T = row;
else
    % Keep only common variables if a previous partial/error row has a
    % slightly different schema.  The normal base_summary_row schema should
    % make this path rare, but it keeps resume robust.
    missingInT = setdiff(row.Properties.VariableNames, T.Properties.VariableNames);
    for i = 1:numel(missingInT)
        T.(missingInT{i}) = repmat(default_column_value(row.(missingInT{i})), height(T), 1);
    end
    missingInRow = setdiff(T.Properties.VariableNames, row.Properties.VariableNames);
    for i = 1:numel(missingInRow)
        row.(missingInRow{i}) = default_column_value(T.(missingInRow{i}));
    end
    row = row(:, T.Properties.VariableNames);
    T = [T; row];
end
end

function v = default_column_value(example)
if isstring(example)
    v = string('');
elseif iscell(example)
    v = {''};
elseif islogical(example)
    v = false;
else
    v = NaN;
end
end

function idx = preferred_method_index(T, methodName)
idx = height(T);
if ismember('method', T.Properties.VariableNames)
    m = T.method;
    if iscell(m), m = string(m); end
    hit = find(strcmpi(string(m), methodName), 1, 'first');
    if ~isempty(hit), idx = hit; end
end
end

function val = get_field(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
end
end

function x = column_or_nan(T, name)
if istable(T) && ismember(name, T.Properties.VariableNames)
    x = T.(name);
else
    x = nan(height(T),1);
end
end

function val = tail_value(T, name)
x = column_or_nan(T, name);
idx = find(isfinite(x), 1, 'last');
if isempty(idx), val = NaN; else, val = x(idx); end
end

function x = default_if_nan(x, defaultValue)
if isempty(x) || (isnumeric(x) && isscalar(x) && isnan(x))
    x = defaultValue;
end
end

function s = value_to_string(val)
if isnumeric(val) || islogical(val)
    if isscalar(val)
        s = num2str(val);
    else
        s = mat2str(val);
    end
elseif ischar(val)
    s = val;
elseif isstring(val)
    s = char(strjoin(val, ','));
elseif iscell(val)
    try
        s = strjoin(cellfun(@char, val, 'UniformOutput', false), ',');
    catch
        s = '<cell>';
    end
else
    s = class(val);
end
end

function local_diary_off()
try
    diary off;
catch
end
end
