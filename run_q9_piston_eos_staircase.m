%% run_q9_piston_eos_staircase.m
% Quasi-static staircase protocol for Q9 effective-EOS diagnostics.
%
% Instead of using a slow moving piston, this script imposes a sequence of
% target density ratios by changing the active piston height geometrically:
%
%   rho/rho0 = yTop0 / yTop.
%
% At each plateau, particle positions are rescaled affinely in y, the piston
% is then held fixed, and pressure-like quantities are averaged only over the
% final part of the plateau.
%
% Main pressure definitions:
%   Pkin    = mean_cell(rho_cell * kBT_cell)
%   Pwall   = top-wall mechanical pressure from accumulated wall impulse
%   Pexcess = Pwall - Pkin
%
% Default usage:
%   run_q9_piston_eos_staircase
%
% Optional overrides before calling this script:
%   pistonEosRhoRatioList = [1.00 1.02 1.04 1.06];
%   pistonEosUserParams = struct('stepsPerPlateau',4000);
%   run_q9_piston_eos_staircase

clear functions
close all
% clc

%% === Optional user controls ===
if ~exist('pistonEosOutputPrefix', 'var') || isempty(pistonEosOutputPrefix)
    pistonEosOutputPrefix = 'q9_piston_eos_staircase';
end
if ~exist('pistonEosIncludeDecompression', 'var') || isempty(pistonEosIncludeDecompression)
    pistonEosIncludeDecompression = false;
end

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = [char(string(pistonEosOutputPrefix)) '_' tag];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on
cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n=== Q9 piston EOS staircase / plateau protocol ===\n');

%% === Base parameters ===
params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Ly0 = params.Ly;
params.Nx = 32;
params.Ny = 16;
params.gamma = 20;

params.dt = 1.0e-3;
params.kBT = 1.0;
params.alphaDeg = 90;
params.bodyForceX = 0.0;
params.bodyForceY = 0.0;
params.seed = 11;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;

%% === Plateau protocol ===
if exist('pistonEosRhoRatioList', 'var') && ~isempty(pistonEosRhoRatioList)
    rhoRatioUp = pistonEosRhoRatioList(:).';
else
    rhoRatioUp = [1.00 1.02 1.04 1.06 1.08 1.10];
end

if pistonEosIncludeDecompression
    if numel(rhoRatioUp) >= 2
        rhoRatioList = [rhoRatioUp, fliplr(rhoRatioUp(1:end-1))];
    else
        rhoRatioList = rhoRatioUp;
    end
else
    rhoRatioList = rhoRatioUp;
end

params.stepsPerPlateau = 2000;
params.sampleEvery = 50;
params.progressEvery = 500;
params.plateauDiscardFraction = 0.50;
params.rhoRatioList = rhoRatioList;
params.pistonY0 = params.Ly;
params.pistonAffineReposition = false;  % transitions are handled explicitly by this script
params.pistonMotionMode = 'linear';
params.pistonVy = 0.0;
params.pistonStopOnMin = true;
params.maxWallClockSeconds = Inf;

%% === Walls / thermostat ===
params.wallModeY = 'thermalize';
params.pistonTangentialMode = 'specular';
params.wallSigma = sqrt(params.kBT);
params.thermostatAfterProjection = true;
params.thermostatTargetKBT = params.kBT;
params.thermostatStrength = 1.0;
params.thermostatMinParticlesPerCell = 2;
params.thermostatMaxScale = 5.0;

%% === Q9 parameters ===
params.projectionStrength = 1.0;
params.projectedStrength = 1.0;
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = 0.002;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxFinalVelocityProjectionCleanup = true;
params.massFluxFinalVelocityProjectionStrength = 0.5;
params.massFluxTargetFilter = 'lowpass_fft';
params.massFluxLowKMaxIndex = 2;
params.massFluxProjectionRegularization = 1e-12;
params.massFluxMinCellCount = 1.0;
params.projectionInterpolationMethod = 'nearest';
params.projectionTransportDiagnosticsEnable = true;
params.densityTransportDiagnosticsEnable = true;
params.computeDiagnostics = true;
params.computeFullDiagnosticsEveryStep = false;

%% === Diagnostics / plotting ===
params.storeDensityMaps = false;
params.densityBandFraction = 0.20;
params.lowKMaxIndex = 2;
params.smoothWindowSamples = 5;
params.makeFigures = true;
params.useCumulativeWallPressureForFit = true;
params.errorBarUseStderr = true;

% User overrides are applied last.
if exist('pistonEosUserParams', 'var') && isstruct(pistonEosUserParams)
    params = merge_struct(params, pistonEosUserParams);
end

params.rhoRatioList = params.rhoRatioList(:).';
params.nPlateaux = numel(params.rhoRatioList);
params.nSteps = params.nPlateaux * params.stepsPerPlateau;

%% === Initialization ===
rng(params.seed);
[state, initInfo] = projection_initialize_particles(params);
state.piston = piston_state(params.pistonY0, params, 0, "initial", 0);

fprintf('Np                         : %d\n', initInfo.Np);
fprintf('grid                       : %d x %d, gamma=%g\n', params.Nx, params.Ny, params.gamma);
fprintf('rhoRatioList               : %s\n', mat2str(params.rhoRatioList, 5));
fprintf('stepsPerPlateau/sampleEvery: %d / %d\n', params.stepsPerPlateau, params.sampleEvery);
fprintf('total steps                : %d\n', params.nSteps);
fprintf('Q9 beta/lowK/cleanup       : %.6g / %d / %.6g\n', ...
    params.massFluxDensityRelaxationBeta, params.massFluxLowKMaxIndex, ...
    params.massFluxFinalVelocityProjectionStrength);

%% === Staircase loop ===
maxSamples = params.nPlateaux * (floor(params.stepsPerPlateau / params.sampleEvery) + 2) + 1;
rows = repmat(empty_row(), maxSamples, 1);
isamp = 0;
globalStep = 0;
actualLastStep = 0;
wallClock = tic;

% Initial sample at the initial volume, before any relaxation.
isamp = isamp + 1;
rows(isamp) = make_sample_row(0, 0, 0, 0.0, "initial", state, params, [], empty_accumulator());

currentYTop = params.pistonY0;
for ip = 1:params.nPlateaux
    rhoTarget = params.rhoRatioList(ip);
    yTopTarget = params.pistonY0 / rhoTarget;
    branch = branch_name_for_plateau(params.rhoRatioList, ip);

    state = affine_set_active_height(state, currentYTop, yTopTarget, params);
    state.piston = piston_state(yTopTarget, params, globalStep, branch, ip);
    currentYTop = yTopTarget;

    plateauParams = params;
    plateauParams.pistonMotionMode = 'linear';
    plateauParams.pistonY0 = yTopTarget;
    plateauParams.pistonYMin = yTopTarget;
    plateauParams.pistonVy = 0.0;
    plateauParams.pistonStopOnMin = true;

    acc = empty_accumulator();
    for k = 1:params.stepsPerPlateau
        globalStep = globalStep + 1;
        [state, stepDiag] = mpcd_step_projection_piston(state, plateauParams, k);
        state.piston.phase = branch;
        state.piston.phaseIndex = ip;
        state.piston.stepIndex = globalStep;
        actualLastStep = globalStep;
        acc = accumulate_step(acc, stepDiag, plateauParams);

        if mod(k, params.sampleEvery) == 0 || k == params.stepsPerPlateau
            isamp = isamp + 1;
            rows(isamp) = make_sample_row(globalStep, ip, k, globalStep*params.dt, branch, state, plateauParams, stepDiag, acc);
            acc = empty_accumulator();
        end

        if params.progressEvery > 0 && mod(k, params.progressEvery) == 0
            fprintf('  plateau %d/%d, local step %d/%d, rho/rho0=%.5g, yTop=%.6g, Pwall=%.4g, Pkin=%.4g, elapsed=%.1fs\n', ...
                ip, params.nPlateaux, k, params.stepsPerPlateau, rhoTarget, yTopTarget, ...
                finite_or_nan(getf(getf(stepDiag,'wallInfo',struct()), 'pressureTopWall')), rows(isamp).Pkin, toc(wallClock));
        end

        if isfinite(params.maxWallClockSeconds) && toc(wallClock) > params.maxWallClockSeconds
            warning('run_q9_piston_eos_staircase:TimeLimit', ...
                'Stopping at step %d/%d after %.1fs.', globalStep, params.nSteps, toc(wallClock));
            break;
        end
    end

    if actualLastStep < globalStep || (isfinite(params.maxWallClockSeconds) && toc(wallClock) > params.maxWallClockSeconds)
        break;
    end
end

elapsedScript = toc(wallClock);
rows = rows(1:isamp);
T = struct2table(rows);
T = add_smoothed_pressures(T, params.smoothWindowSamples);
[plateauSummary, baseline] = summarize_plateaux(T, params);
T = add_relative_pressures(T, baseline);
plateauSummary = add_plateau_relative_pressures(plateauSummary, baseline);
fitSummary = fit_plateau_eos(plateauSummary);

%% === Save ===
matFile = fullfile(outputDir, 'q9_piston_eos_staircase.mat');
csvFile = fullfile(outputDir, 'q9_piston_eos_staircase_timeseries.csv');
plateauCsvFile = fullfile(outputDir, 'q9_piston_eos_staircase_plateau_summary.csv');
fitCsvFile = fullfile(outputDir, 'q9_piston_eos_staircase_fit_summary.csv');
txtFile = fullfile(outputDir, 'q9_piston_eos_staircase_summary.txt');

save(matFile, 'params', 'initInfo', 'T', 'plateauSummary', 'fitSummary', 'baseline', 'actualLastStep', 'elapsedScript', '-v7.3');
writetable(T, csvFile);
writetable(plateauSummary, plateauCsvFile);
writetable(fitSummary, fitCsvFile);
write_staircase_summary(txtFile, params, initInfo, baseline, plateauSummary, fitSummary, matFile, csvFile, plateauCsvFile, fitCsvFile, txtFile, actualLastStep, elapsedScript);

fprintf('\nTimeseries written to:\n%s\n', csvFile);
fprintf('Plateau summary written to:\n%s\n', plateauCsvFile);
fprintf('Fit summary written to:\n%s\n', fitCsvFile);
fprintf('Text summary written to:\n%s\n', txtFile);
fprintf('MAT written to:\n%s\n', matFile);

if params.makeFigures
    make_staircase_figures(T, plateauSummary, fitSummary, outputDir);
end

diary off

%% ========================================================================
function row = empty_row()
row = struct();
row.step = NaN;
row.t = NaN;
row.plateauIndex = NaN;
row.plateauStep = NaN;
row.branch = "";
row.targetRhoRatio = NaN;
row.yTop = NaN;
row.compression = NaN;
row.volume = NaN;
row.volumeRatio = NaN;
row.rhoMean = NaN;
row.rhoRatio = NaN;
row.Pkin = NaN;
row.Pwall = NaN;
row.PexcessWallKinetic = NaN;
row.PkinSmooth = NaN;
row.PwallSmooth = NaN;
row.PexcessWallKineticSmooth = NaN;
row.PkinRel = NaN;
row.PwallRel = NaN;
row.PexcessRel = NaN;
row.rhoRatioRel = NaN;
row.wallImpulseTopY = NaN;
row.wallHitsTop = NaN;
row.sampleSteps = NaN;
row.sampleDt = NaN;
row.kBTCell = NaN;
row.stdN = NaN;
row.lowKProxy = NaN;
row.rmsDivParticleAfter = NaN;
row.rmsMassFluxResidual = NaN;
row.rmsMassFluxParticleAfter = NaN;
row.meanVy = NaN;
end

function acc = empty_accumulator()
acc.n = 0;
acc.dt = 0.0;
acc.wallImpulseTopY = 0.0;
acc.wallHitsTop = 0.0;
end

function acc = accumulate_step(acc, stepDiag, params)
acc.n = acc.n + 1;
acc.dt = acc.dt + params.dt;
W = getf(stepDiag, 'wallInfo', struct());
acc.wallImpulseTopY = acc.wallImpulseTopY + getf(W, 'impulseOnTopWallY', signed_top_impulse_fallback(W));
acc.wallHitsTop = acc.wallHitsTop + getf(W, 'nTop', 0);
end

function row = make_sample_row(globalStep, plateauIndex, plateauStep, t, branch, state, params, stepDiag, acc)
activeParams = params;
yTop = params.pistonY0;
if isfield(state, 'piston') && isfield(state.piston, 'yTop')
    yTop = state.piston.yTop;
end
activeParams.Ly = yTop;

G = projection_deposit_particles_to_grid(state.x, state.v, activeParams, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
thermal = projection_thermal_diagnostics(state.x, state.v, activeParams, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
Tmap = local_temperature_map(state.x, state.v, activeParams);
Pkin = mean(G.rho(:) .* max(Tmap(:), 0), 'omitnan');
if ~isfinite(Pkin)
    Pkin = size(state.x,1) / max(params.Lx*yTop, eps) * thermal.kBTCellRelative;
end

volume = params.Lx * yTop;
rhoMean = size(state.x,1) / max(volume, eps);
rho0 = params.Nx * params.Ny * params.gamma / max(params.Lx * params.Ly0, eps);

if acc.n > 0 && acc.dt > 0
    Pwall = acc.wallImpulseTopY / max(acc.dt * params.Lx, eps);
else
    Pwall = NaN;
end

row = empty_row();
row.step = globalStep;
row.t = t;
row.plateauIndex = plateauIndex;
row.plateauStep = plateauStep;
row.branch = branch;
row.targetRhoRatio = params.pistonY0 / max(yTop, eps);
row.yTop = yTop;
row.compression = 1.0 - yTop / max(params.pistonY0, eps);
row.volume = volume;
row.volumeRatio = volume / max(params.Lx * params.Ly0, eps);
row.rhoMean = rhoMean;
row.rhoRatio = rhoMean / max(rho0, eps);
row.Pkin = Pkin;
row.Pwall = Pwall;
row.PexcessWallKinetic = Pwall - Pkin;
row.wallImpulseTopY = acc.wallImpulseTopY;
row.wallHitsTop = acc.wallHitsTop;
row.sampleSteps = acc.n;
row.sampleDt = acc.dt;
row.kBTCell = thermal.kBTCellRelative;
row.stdN = std(double(G.N(:)));
row.lowKProxy = lowk_density_proxy(G.N, params.gamma, params.lowKMaxIndex);
row.meanVy = mean(state.v(:,2));
if ~isempty(stepDiag)
    row.rmsDivParticleAfter = getf(stepDiag, 'rmsDivParticleAfter');
    row.rmsMassFluxResidual = getf(stepDiag, 'rmsMassFluxDivResidual');
    row.rmsMassFluxParticleAfter = getf(stepDiag, 'rmsMassFluxDivParticleAfter');
end
end

function T = add_smoothed_pressures(T, w)
w = max(1, round(w));
T.PkinSmooth = movmean(T.Pkin, w, 'omitnan');
T.PwallSmooth = movmean(T.Pwall, w, 'omitnan');
T.PexcessWallKineticSmooth = movmean(T.PexcessWallKinetic, w, 'omitnan');
end

function [S, baseline] = summarize_plateaux(T, params)
plateaux = unique(T.plateauIndex(~isnan(T.plateauIndex) & T.plateauIndex > 0));
rows = repmat(empty_plateau_row(), numel(plateaux), 1);
for m = 1:numel(plateaux)
    ip = plateaux(m);
    idxAll = find(T.plateauIndex == ip);
    idxAll = idxAll(T.plateauStep(idxAll) > 0);
    nDiscard = floor(params.plateauDiscardFraction * numel(idxAll));
    idx = idxAll((nDiscard+1):end);

    r = empty_plateau_row();
    r.plateauIndex = ip;
    r.branch = T.branch(idxAll(end));
    r.targetRhoRatio = mean(T.targetRhoRatio(idx), 'omitnan');
    r.yTop = mean(T.yTop(idx), 'omitnan');
    r.compression = mean(T.compression(idx), 'omitnan');
    r.nSamplesKept = numel(idx);
    r.nStepsRetained = sum(T.sampleSteps(idx), 'omitnan');
    r.retainedTime = sum(T.sampleDt(idx), 'omitnan');
    r.meanRhoRatio = mean(T.rhoRatio(idx), 'omitnan');

    % Kinetic pressure: volume/cell average, sampled on retained windows.
    [r.meanPkin, r.stdPkin, r.stderrPkin] = mean_std_stderr(T.Pkin(idx));

    % Wall pressure has two estimates:
    %   meanPwallSample      : mean of sample-window pressures;
    %   meanPwallCumulative  : total retained impulse / total retained time / Lx.
    % The cumulative estimator is less sensitive to uneven sample-window noise and
    % is used by default for EOS fits and relative pressure diagnostics.
    [r.meanPwallSample, r.stdPwall, r.stderrPwall] = mean_std_stderr(T.Pwall(idx));
    totalImpulseTop = sum(T.wallImpulseTopY(idx), 'omitnan');
    totalTime = sum(T.sampleDt(idx), 'omitnan');
    r.totalWallImpulseTopY = totalImpulseTop;
    r.totalWallHitsTop = sum(T.wallHitsTop(idx), 'omitnan');
    if totalTime > 0
        r.meanPwallCumulative = totalImpulseTop / max(totalTime * params.Lx, eps);
    else
        r.meanPwallCumulative = NaN;
    end

    if get_param(params, 'useCumulativeWallPressureForFit', true)
        r.meanPwall = r.meanPwallCumulative;
    else
        r.meanPwall = r.meanPwallSample;
    end

    % Excess pressure uses the same wall-pressure convention as meanPwall.
    PexcessSample = T.Pwall(idx) - T.Pkin(idx);
    [r.meanPexcessSample, r.stdPexcess, r.stderrPexcess] = mean_std_stderr(PexcessSample);
    r.meanPexcessCumulative = r.meanPwallCumulative - r.meanPkin;
    if get_param(params, 'useCumulativeWallPressureForFit', true)
        r.meanPexcess = r.meanPexcessCumulative;
    else
        r.meanPexcess = r.meanPexcessSample;
    end

    r.meanStdN = mean(T.stdN(idx), 'omitnan');
    r.meanLowKProxy = mean(T.lowKProxy(idx), 'omitnan');
    r.meanDivUAfter = mean(T.rmsDivParticleAfter(idx), 'omitnan');
    r.meanMassFluxResidual = mean(T.rmsMassFluxResidual(idx), 'omitnan');
    r.meanMassFluxAfter = mean(T.rmsMassFluxParticleAfter(idx), 'omitnan');
    r.meanKBTCell = mean(T.kBTCell(idx), 'omitnan');
    r.meanVy = mean(T.meanVy(idx), 'omitnan');
    rows(m) = r;
end
S = struct2table(rows);

idxBase = find(S.plateauIndex == min(S.plateauIndex));
if isempty(idxBase)
    idxBase = 1;
end
baseline = struct();
baseline.plateauIndex = S.plateauIndex(idxBase(1));
baseline.rhoRatio = S.meanRhoRatio(idxBase(1));
baseline.Pkin = S.meanPkin(idxBase(1));
baseline.Pwall = S.meanPwall(idxBase(1));
baseline.PwallSample = S.meanPwallSample(idxBase(1));
baseline.PwallCumulative = S.meanPwallCumulative(idxBase(1));
baseline.Pexcess = S.meanPexcess(idxBase(1));
baseline.PexcessSample = S.meanPexcessSample(idxBase(1));
baseline.PexcessCumulative = S.meanPexcessCumulative(idxBase(1));
end

function r = empty_plateau_row()
r = struct();
r.plateauIndex = NaN;
r.branch = "";
r.targetRhoRatio = NaN;
r.yTop = NaN;
r.compression = NaN;
r.nSamplesKept = NaN;
r.nStepsRetained = NaN;
r.retainedTime = NaN;
r.meanRhoRatio = NaN;
r.meanPkin = NaN;
r.stdPkin = NaN;
r.stderrPkin = NaN;
r.meanPwall = NaN;
r.meanPwallSample = NaN;
r.meanPwallCumulative = NaN;
r.stdPwall = NaN;
r.stderrPwall = NaN;
r.totalWallImpulseTopY = NaN;
r.totalWallHitsTop = NaN;
r.meanPexcess = NaN;
r.meanPexcessSample = NaN;
r.meanPexcessCumulative = NaN;
r.stdPexcess = NaN;
r.stderrPexcess = NaN;
r.PkinRel = NaN;
r.PwallRel = NaN;
r.PwallSampleRel = NaN;
r.PwallCumulativeRel = NaN;
r.PexcessRel = NaN;
r.PexcessSampleRel = NaN;
r.PexcessCumulativeRel = NaN;
r.rhoRatioRel = NaN;
r.meanStdN = NaN;
r.meanLowKProxy = NaN;
r.meanDivUAfter = NaN;
r.meanMassFluxResidual = NaN;
r.meanMassFluxAfter = NaN;
r.meanKBTCell = NaN;
r.meanVy = NaN;
end

function T = add_relative_pressures(T, baseline)
T.rhoRatioRel = T.rhoRatio - baseline.rhoRatio;
T.PkinRel = T.PkinSmooth - baseline.Pkin;
T.PwallRel = T.PwallSmooth - baseline.Pwall;
T.PexcessRel = T.PexcessWallKineticSmooth - baseline.Pexcess;
end

function S = add_plateau_relative_pressures(S, baseline)
S.rhoRatioRel = S.meanRhoRatio - baseline.rhoRatio;
S.PkinRel = S.meanPkin - baseline.Pkin;
S.PwallRel = S.meanPwall - baseline.Pwall;
S.PwallSampleRel = S.meanPwallSample - baseline.PwallSample;
S.PwallCumulativeRel = S.meanPwallCumulative - baseline.PwallCumulative;
S.PexcessRel = S.meanPexcess - baseline.Pexcess;
S.PexcessSampleRel = S.meanPexcessSample - baseline.PexcessSample;
S.PexcessCumulativeRel = S.meanPexcessCumulative - baseline.PexcessCumulative;
end

function [mu, sig, se] = mean_std_stderr(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    mu = NaN; sig = NaN; se = NaN;
    return;
end
mu = mean(x);
if numel(x) >= 2
    sig = std(x);
else
    sig = 0;
end
se = sig / sqrt(max(numel(x), 1));
end

function F = fit_plateau_eos(S)
branches = unique(S.branch);
rows = repmat(empty_fit_row(), 0, 1);
for i = 1:numel(branches)
    br = branches(i);
    idx = S.branch == br & isfinite(S.meanRhoRatio) & isfinite(S.PkinRel) & S.nSamplesKept >= 2;
    if nnz(idx) >= 3
        rows(end+1) = fit_one_branch(S(idx,:), br); %#ok<AGROW>
    end
end
idxAll = isfinite(S.meanRhoRatio) & isfinite(S.PkinRel) & S.nSamplesKept >= 2;
if nnz(idxAll) >= 3
    rows(end+1) = fit_one_branch(S(idxAll,:), "all_plateaux"); %#ok<AGROW>
end
if isempty(rows)
    rows = empty_fit_row();
end
F = struct2table(rows);
end

function r = empty_fit_row()
r = struct();
r.branch = "";
r.nPlateaux = NaN;
r.rhoMin = NaN;
r.rhoMax = NaN;
r.slopePkin = NaN;
r.slopePwall = NaN;
r.slopePexcess = NaN;
r.KeffKinetic = NaN;
r.KeffWall = NaN;
r.KeffExcess = NaN;
r.R2Pkin = NaN;
r.R2Pwall = NaN;
r.R2Pexcess = NaN;
end

function r = fit_one_branch(S, br)
x = S.meanRhoRatio;
x0 = mean(x, 'omitnan');
r = empty_fit_row();
r.branch = br;
r.nPlateaux = height(S);
r.rhoMin = min(x);
r.rhoMax = max(x);
[r.slopePkin, r.R2Pkin] = linear_fit_slope_r2(x, S.PkinRel);
[r.slopePwall, r.R2Pwall] = linear_fit_slope_r2(x, S.PwallRel);
[r.slopePexcess, r.R2Pexcess] = linear_fit_slope_r2(x, S.PexcessRel);
r.KeffKinetic = x0 * r.slopePkin;
r.KeffWall = x0 * r.slopePwall;
r.KeffExcess = x0 * r.slopePexcess;
end

function [slope, r2] = linear_fit_slope_r2(x, y)
mask = isfinite(x) & isfinite(y);
x = x(mask); y = y(mask);
if numel(x) < 2
    slope = NaN; r2 = NaN; return;
end
p = polyfit(x(:), y(:), 1);
yfit = polyval(p, x(:));
slope = p(1);
ssRes = sum((y(:)-yfit(:)).^2);
ssTot = sum((y(:)-mean(y(:))).^2);
if ssTot > 0
    r2 = 1 - ssRes/ssTot;
else
    r2 = NaN;
end
end

function state = affine_set_active_height(state, yOld, yNew, params)
scaleY = yNew / max(yOld, eps);
state.x(:,2) = state.x(:,2) * scaleY;
% Guard against numerical wall overshoots after the affine map.
yeps = max(1e-12, 10*eps(max(yNew,1)));
state.x(:,2) = min(max(state.x(:,2), yeps), yNew - yeps);
state.x(:,1) = mod(state.x(:,1), params.Lx);
end

function p = piston_state(yTop, params, globalStep, branch, phaseIndex)
p = struct();
p.yTop = yTop;
p.yPrev = yTop;
p.Up = 0.0;
p.Ly0 = params.Ly0;
p.activeHeight = yTop;
p.compression = 1.0 - yTop / max(params.pistonY0, eps);
p.stepIndex = globalStep;
p.t = globalStep * params.dt;
p.motionMode = 'staircase';
p.phase = branch;
p.phaseIndex = phaseIndex;
p.dV = 0.0;
p.dVdt = 0.0;
p.activeArea = params.Lx * yTop;
p.rhoPhysicalMean = NaN;
p.gammaMeanGeometric = params.gamma;
end

function br = branch_name_for_plateau(rhoList, ip)
if ip <= 1
    br = "compression";
elseif rhoList(ip) >= rhoList(ip-1)
    br = "compression";
else
    br = "decompression";
end
end

function T = local_temperature_map(x, v, params)
G = projection_deposit_particles_to_grid(x, v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
Nx = params.Nx;
Ny = params.Ny;
dx = params.Lx / Nx;
dy = params.Ly / Ny;
ix = floor(mod(x(:,1), params.Lx) / dx) + 1;
iy = floor(min(max(x(:,2), 0), params.Ly - eps(params.Ly)) / dy) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
cellId = iy + Ny*(ix-1);
UxVec = reshape(G.Ux.', [], 1);
UyVec = reshape(G.Uy.', [], 1);
NVec = reshape(G.N.', [], 1);
ux = UxVec(cellId);
uy = UyVec(cellId);
rel2 = (v(:,1)-ux).^2 + (v(:,2)-uy).^2;
sumRel2 = accumarray(cellId, rel2, [Nx*Ny, 1], @sum, 0);
dof = 2*max(double(NVec)-1, 0);
Tvec = nan(Nx*Ny,1);
mask = dof > 0;
Tvec(mask) = sumRel2(mask)./dof(mask);
T = reshape(Tvec, [Ny, Nx]).';
end

function e = lowk_density_proxy(N, gamma, kmax)
F = fft2(double(N) / max(gamma, eps) - 1.0);
E = abs(F).^2 / numel(F)^2;
mask = false(size(E));
Nx = size(E,1); Ny = size(E,2);
for i = 1:Nx
    ki = min(i-1, Nx-(i-1));
    for j = 1:Ny
        kj = min(j-1, Ny-(j-1));
        if ki <= kmax && kj <= kmax && ~(ki == 0 && kj == 0)
            mask(i,j) = true;
        end
    end
end
e = sum(E(mask));
end

function make_staircase_figures(T, S, F, outputDir)
fig1 = figure('Name', 'Piston EOS staircase: absolute pressure versus density');
hold on;
errorbar(S.meanRhoRatio, S.meanPkin, S.stderrPkin, 'o-');
errorbar(S.meanRhoRatio, S.meanPwall, S.stderrPwall, 'o-');
errorbar(S.meanRhoRatio, S.meanPexcess, S.stderrPexcess, 'o-');
grid on; xlabel('\rho / \rho_0'); ylabel('pressure-like quantity');
legend('Pkin','Pwall cumulative','Pwall-Pkin','Location','best');
title('Plateau-averaged absolute pressures with standard errors');
saveas(fig1, fullfile(outputDir, 'piston_eos_staircase_pressure_vs_density.png'));

fig1b = figure('Name', 'Piston EOS staircase: sampled vs cumulative wall pressure');
plot(S.meanRhoRatio, S.meanPwallSample, 'o-', S.meanRhoRatio, S.meanPwallCumulative, 's-');
grid on; xlabel('\rho / \rho_0'); ylabel('Pwall');
legend('sample mean','cumulative impulse','Location','best');
title('Wall pressure estimators');
saveas(fig1b, fullfile(outputDir, 'piston_eos_staircase_wall_pressure_estimators.png'));

fig2 = figure('Name', 'Piston EOS staircase: relative pressure versus density');
hold on;
errorbar(S.meanRhoRatio, S.PkinRel, S.stderrPkin, 'o-');
errorbar(S.meanRhoRatio, S.PwallRel, S.stderrPwall, 'o-');
errorbar(S.meanRhoRatio, S.PexcessRel, S.stderrPexcess, 'o-');
grid on; xlabel('\rho / \rho_0'); ylabel('\Delta pressure-like quantity');
legend('\Delta Pkin','\Delta Pwall cumulative','\Delta(Pwall-Pkin)','Location','best');
title('Plateau-averaged relative pressures with standard errors');
saveas(fig2, fullfile(outputDir, 'piston_eos_staircase_relative_pressure_vs_density.png'));

fig2b = figure('Name', 'Piston EOS staircase: error bars');
semilogy(S.meanRhoRatio, max(S.stderrPkin, eps), 'o-', S.meanRhoRatio, max(S.stderrPwall, eps), 'o-', S.meanRhoRatio, max(S.stderrPexcess, eps), 'o-');
grid on; xlabel('\rho / \rho_0'); ylabel('standard error');
legend('Pkin','Pwall','Pexcess','Location','best');
title('Plateau pressure standard errors');
saveas(fig2b, fullfile(outputDir, 'piston_eos_staircase_pressure_standard_errors.png'));

fig3 = figure('Name', 'Piston EOS staircase: time series');
tiledlayout(4,1);
nexttile; plot(T.t, T.rhoRatio, '-'); grid on; ylabel('\rho/\rho_0');
nexttile; plot(T.t, T.PkinSmooth, '-', T.t, T.PwallSmooth, '-'); grid on; ylabel('P'); legend('Pkin','Pwall sample','Location','best');
nexttile; plot(T.t, T.PexcessWallKineticSmooth, '-'); grid on; ylabel('Pwall-Pkin');
nexttile; plot(T.t, T.stdN, '-', T.t, T.lowKProxy, '-'); grid on; ylabel('density diag'); xlabel('t'); legend('stdN','low-k proxy','Location','best');
saveas(fig3, fullfile(outputDir, 'piston_eos_staircase_timeseries.png'));

if height(F) > 0
    fig4 = figure('Name', 'Piston EOS staircase: fitted moduli');
    cats = categorical(cellstr(F.branch));
    bar(cats, [F.KeffKinetic, F.KeffWall, F.KeffExcess]);
    grid on; ylabel('K_{eff} fit'); legend('kinetic','wall','excess','Location','best');
    title('Effective compressibility moduli from plateau fits');
    saveas(fig4, fullfile(outputDir, 'piston_eos_staircase_keff_summary.png'));
end
end

function write_staircase_summary(txtFile, params, initInfo, baseline, S, F, matFile, csvFile, plateauCsvFile, fitCsvFile, txtFile2, actualLastStep, elapsedScript)
fid = fopen(txtFile, 'w');
if fid < 0
    warning('Could not open summary text file: %s', txtFile);
    return;
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '=== Q9 piston EOS staircase / plateau protocol ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'Np                                : %d\n', initInfo.Np);
fprintf(fid, 'steps / stepsPerPlateau / sampleEvery : %d / %d / %d\n', params.nSteps, params.stepsPerPlateau, params.sampleEvery);
fprintf(fid, 'actualLastStep                    : %d\n', actualLastStep);
fprintf(fid, 'elapsed seconds                   : %.3f\n', elapsedScript);
fprintf(fid, 'grid, gamma                       : %d x %d, %g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'rhoRatioList                      : %s\n', mat2str(params.rhoRatioList, 6));
fprintf(fid, 'plateauDiscardFraction            : %.3f\n', params.plateauDiscardFraction);
fprintf(fid, 'dt, kBT, alphaDeg                 : %.6g, %.6g, %.6g\n', params.dt, params.kBT, params.alphaDeg);
fprintf(fid, '\n--- Q9 parameters ---\n');
fprintf(fid, 'projectionStrength                : %.6g\n', params.projectionStrength);
fprintf(fid, 'massFluxProjectionMode            : %s\n', char(string(params.massFluxProjectionMode)));
fprintf(fid, 'massFluxDensityRelaxationBeta     : %.6g\n', params.massFluxDensityRelaxationBeta);
fprintf(fid, 'massFluxFinalVelocityProjectionCleanup : %d\n', logical(params.massFluxFinalVelocityProjectionCleanup));
fprintf(fid, 'massFluxFinalVelocityProjectionStrength : %.6g\n', params.massFluxFinalVelocityProjectionStrength);
fprintf(fid, 'massFluxLowKMaxIndex              : %d\n', params.massFluxLowKMaxIndex);
fprintf(fid, '\n--- Pressure definitions ---\n');
fprintf(fid, 'Pkin    = mean_cell(rho_cell * kBT_cell)\n');
fprintf(fid, 'Pwall   = top-wall impulse / (Delta_t_sample * Lx)\n');
fprintf(fid, 'Pwall in plateau summaries uses the cumulative retained impulse by default.\n');
fprintf(fid, 'Pexcess = Pwall - Pkin\n');
fprintf(fid, 'std/stderr are computed from retained sample-window pressures.\n');
fprintf(fid, 'Relative pressures are baseline-subtracted using the first plateau.\n');
fprintf(fid, '\n--- Baseline plateau ---\n');
fprintf(fid, 'plateauIndex                      : %g\n', baseline.plateauIndex);
fprintf(fid, 'rho/rho0                          : %.12g\n', baseline.rhoRatio);
fprintf(fid, 'Pkin / Pwall / Pexcess            : %.12g / %.12g / %.12g\n', baseline.Pkin, baseline.Pwall, baseline.Pexcess);
fprintf(fid, '\n--- Plateau summary ---\n');
for i = 1:height(S)
    fprintf(fid, 'plateau %d [%s], rho/rho0=%.6g, nSamples=%d, PkinRel=%.6g +/- %.3g, PwallRel=%.6g +/- %.3g, PexcessRel=%.6g +/- %.3g, stdN=%.6g, lowK=%.6g\n', ...
        S.plateauIndex(i), char(S.branch(i)), S.meanRhoRatio(i), S.nSamplesKept(i), ...
        S.PkinRel(i), S.stderrPkin(i), S.PwallRel(i), S.stderrPwall(i), S.PexcessRel(i), S.stderrPexcess(i), ...
        S.meanStdN(i), S.meanLowKProxy(i));
    fprintf(fid, '    Pwall sample/cumulative = %.12g / %.12g, hits=%g, retainedTime=%.6g\n', ...
        S.meanPwallSample(i), S.meanPwallCumulative(i), S.totalWallHitsTop(i), S.retainedTime(i));
end
fprintf(fid, '\n--- Fit summary ---\n');
for i = 1:height(F)
    fprintf(fid, 'branch: %s\n', char(F.branch(i)));
    fprintf(fid, '  nPlateaux                       : %d\n', F.nPlateaux(i));
    fprintf(fid, '  Keff kinetic / wall / excess    : %.12g / %.12g / %.12g\n', F.KeffKinetic(i), F.KeffWall(i), F.KeffExcess(i));
    fprintf(fid, '  R2 kinetic / wall / excess      : %.6g / %.6g / %.6g\n', F.R2Pkin(i), F.R2Pwall(i), F.R2Pexcess(i));
end
fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n%s\n%s\n%s\n%s\n', matFile, csvFile, plateauCsvFile, fitCsvFile, txtFile2);
end

function v = signed_top_impulse_fallback(W)
if isstruct(W) && isfield(W, 'dPyTopSigned') && isfinite(W.dPyTopSigned)
    v = -W.dPyTopSigned;
elseif isstruct(W) && isfield(W, 'dPyTop') && isfinite(W.dPyTop)
    v = W.dPyTop;
else
    v = 0.0;
end
end

function y = finite_or_nan(x)
if isempty(x) || ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
    y = NaN;
else
    y = x;
end
end

function value = getf(s, name, defaultValue)
if nargin < 3
    defaultValue = NaN;
end
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end

function out = merge_struct(a, b)
out = a;
fn = fieldnames(b);
for i = 1:numel(fn)
    out.(fn{i}) = b.(fn{i});
end
end
