%% run_q9_piston_eos_beta_sweep.m
% Sweep massFluxDensityRelaxationBeta and estimate the effective EOS modulus.
%
% This script repeatedly calls run_q9_piston_eos_staircase with the resolved
% three-plateau protocol:
%
%   rho/rho0 = [1.00 1.05 1.10]
%   stepsPerPlateau = 6000
%   sampleEvery = 50
%   discard first 50% of each plateau
%
% It aggregates, for each beta, the fitted moduli:
%
%   Kkin, Kwall, Kexcess = Kwall - Kkin-like fit
%
% plus fit quality, pressure standard errors, density diagnostics and output
% folder names.
%
% Default usage:
%   run_q9_piston_eos_beta_sweep
%
% Optional overrides before running:
%   pistonEosBetaList = [0 5e-4 1e-3 2e-3];
%   pistonEosSweepUserParams = struct('stepsPerPlateau',4000);
%   run_q9_piston_eos_beta_sweep

clear functions
close all
% clc

repoRoot = pwd;
addpath(repoRoot);

%% === Sweep controls ===
if exist('pistonEosBetaList', 'var') && ~isempty(pistonEosBetaList)
    betaList = pistonEosBetaList(:).';
else
    % Practical default: includes no Q9 mass-flux stiffness, weak settings,
    % the current piston reference beta=0.002, and one stronger point.
    betaList = [0, 2e-4, 5e-4, 1e-3, 2e-3, 4e-3];
end

rhoRatioList = [1.00 1.05 1.10];
stepsPerPlateau = 6000;
sampleEvery = 50;
discardFraction = 0.50;

betaList = [5e-4, 1e-3, 2e-3];

rhoRatioList = [1.00, 1.05, 1.10];

params.stepsPerPlateau = 6000;
params.sampleEvery = 50;
params.plateauDiscardFraction = 0.50;


if exist('pistonEosSweepUserParams', 'var') && isstruct(pistonEosSweepUserParams)
    if isfield(pistonEosSweepUserParams, 'rhoRatioList')
        rhoRatioList = pistonEosSweepUserParams.rhoRatioList(:).';
    end
    if isfield(pistonEosSweepUserParams, 'stepsPerPlateau')
        stepsPerPlateau = pistonEosSweepUserParams.stepsPerPlateau;
    end
    if isfield(pistonEosSweepUserParams, 'sampleEvery')
        sampleEvery = pistonEosSweepUserParams.sampleEvery;
    end
    if isfield(pistonEosSweepUserParams, 'plateauDiscardFraction')
        discardFraction = pistonEosSweepUserParams.plateauDiscardFraction;
    end
end

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
sweepDir = ['q9_piston_eos_beta_sweep_' tag];
if ~exist(sweepDir, 'dir')
    mkdir(sweepDir);
end

fprintf('\n=== Q9 piston EOS beta sweep ===\n');
fprintf('Output folder       : %s\n', sweepDir);
fprintf('betaList            : %s\n', mat2str(betaList, 6));
fprintf('rhoRatioList        : %s\n', mat2str(rhoRatioList, 6));
fprintf('stepsPerPlateau     : %d\n', stepsPerPlateau);
fprintf('sampleEvery         : %d\n', sampleEvery);
fprintf('discardFraction     : %.3f\n', discardFraction);

oldDir = pwd;
restorePathObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>
cd(sweepDir);
addpath(repoRoot);

sweepRows = repmat(empty_sweep_row(), numel(betaList), 1);
tSweep = tic;

for ib = 1:numel(betaList)
    beta = betaList(ib);
    fprintf('\n--- Sweep beta %d/%d: beta = %.8g ---\n', ib, numel(betaList), beta);
    tOne = tic;

    % Variables consumed by run_q9_piston_eos_staircase.
    pistonEosOutputPrefix = ['q9_piston_eos_beta_' beta_tag(beta)]; %#ok<NASGU>
    pistonEosIncludeDecompression = false; %#ok<NASGU>
    pistonEosRhoRatioList = rhoRatioList; %#ok<NASGU>

    pistonEosUserParams = struct(); %#ok<NASGU>
    pistonEosUserParams.stepsPerPlateau = stepsPerPlateau;
    pistonEosUserParams.sampleEvery = sampleEvery;
    pistonEosUserParams.plateauDiscardFraction = discardFraction;
    pistonEosUserParams.useCumulativeWallPressureForFit = true;
    pistonEosUserParams.makeFigures = false;
    pistonEosUserParams.massFluxDensityRelaxationBeta = beta;
    pistonEosUserParams.massFluxLowKMaxIndex = 2;
    pistonEosUserParams.massFluxFinalVelocityProjectionCleanup = true;
    pistonEosUserParams.massFluxFinalVelocityProjectionStrength = 0.5;
    pistonEosUserParams.massFluxProjectionStrength = 1.0;

    if exist('pistonEosSweepUserParams', 'var') && isstruct(pistonEosSweepUserParams)
        % Apply additional user overrides last, except protocol fields that
        % were already lifted to explicit wrapper variables above.
        tmp = pistonEosSweepUserParams;
        tmp = rmfield_if_present(tmp, {'rhoRatioList','stepsPerPlateau','sampleEvery','plateauDiscardFraction'});
        pistonEosUserParams = merge_struct_local(pistonEosUserParams, tmp); %#ok<NASGU>
    end

    try
        run_q9_piston_eos_staircase
        elapsedOne = toc(tOne);
        sweepRows(ib) = make_sweep_row(beta, ib, outputDir, actualLastStep, elapsedOne, params, plateauSummary, fitSummary, 'ok', '');
        fprintf('beta %.8g done: Kwall=%.6g, Kkin=%.6g, Kextra=%.6g, R2wall=%.4f, elapsed=%.1fs\n', ...
            beta, sweepRows(ib).KeffWall, sweepRows(ib).KeffKinetic, sweepRows(ib).KeffExtraFit, sweepRows(ib).R2Pwall, elapsedOne);
    catch ME
        elapsedOne = toc(tOne);
        warning('Beta %.8g failed: %s', beta, ME.message);
        sweepRows(ib) = make_failed_sweep_row(beta, ib, elapsedOne, ME);
    end

    % Save partial results after every beta so interrupted sweeps are usable.
    sweepSummary = struct2table(sweepRows(1:ib)); %#ok<NASGU>
    writetable(sweepSummary, 'q9_piston_eos_beta_sweep_summary_partial.csv');
    save('q9_piston_eos_beta_sweep_partial.mat', 'sweepSummary', 'betaList', 'rhoRatioList', 'stepsPerPlateau', 'sampleEvery', 'discardFraction');
end

elapsedSweep = toc(tSweep);
sweepSummary = struct2table(sweepRows);

summaryCsv = 'q9_piston_eos_beta_sweep_summary.csv';
summaryTxt = 'q9_piston_eos_beta_sweep_summary.txt';
summaryMat = 'q9_piston_eos_beta_sweep.mat';
writetable(sweepSummary, summaryCsv);
save(summaryMat, 'sweepSummary', 'betaList', 'rhoRatioList', 'stepsPerPlateau', 'sampleEvery', 'discardFraction', 'elapsedSweep');
write_beta_sweep_summary(summaryTxt, sweepSummary, betaList, rhoRatioList, stepsPerPlateau, sampleEvery, discardFraction, elapsedSweep);
make_beta_sweep_figures(sweepSummary, pwd);

fprintf('\n=== Beta sweep complete ===\n');
fprintf('Elapsed seconds: %.1f\n', elapsedSweep);
fprintf('Summary CSV    : %s\n', fullfile(pwd, summaryCsv));
fprintf('Summary TXT    : %s\n', fullfile(pwd, summaryTxt));
fprintf('Summary MAT    : %s\n', fullfile(pwd, summaryMat));

%% ========================================================================
function row = empty_sweep_row()
row = struct();
row.sweepIndex = NaN;
row.beta = NaN;
row.status = "";
row.message = "";
row.outputDir = "";
row.actualLastStep = NaN;
row.elapsedSeconds = NaN;
row.stepsPerPlateau = NaN;
row.sampleEvery = NaN;
row.nPlateaux = NaN;
row.rhoMin = NaN;
row.rhoMax = NaN;
row.massFluxLowKMaxIndex = NaN;
row.cleanupStrength = NaN;
row.massFluxProjectionStrength = NaN;
row.KeffKinetic = NaN;
row.KeffWall = NaN;
row.KeffExcess = NaN;
row.KeffExtraFit = NaN;
row.KeffExtraDifference = NaN;
row.R2Pkin = NaN;
row.R2Pwall = NaN;
row.R2Pexcess = NaN;
row.PkinRelMax = NaN;
row.PwallRelMax = NaN;
row.PexcessRelMax = NaN;
row.stderrPkinMax = NaN;
row.stderrPwallMax = NaN;
row.stderrPexcessMax = NaN;
row.meanStdN = NaN;
row.meanLowKProxy = NaN;
row.meanDivUAfter = NaN;
row.meanMassFluxResidual = NaN;
row.meanMassFluxAfter = NaN;
row.meanKBTCell = NaN;
row.meanWallHitsPerPlateau = NaN;
end

function row = make_failed_sweep_row(beta, ib, elapsedOne, ME)
row = empty_sweep_row();
row.sweepIndex = ib;
row.beta = beta;
row.status = "failed";
row.message = string(ME.message);
row.elapsedSeconds = elapsedOne;
end

function row = make_sweep_row(beta, ib, outputDir, actualLastStep, elapsedOne, params, S, F, status, message)
row = empty_sweep_row();
row.sweepIndex = ib;
row.beta = beta;
row.status = string(status);
row.message = string(message);
row.outputDir = string(outputDir);
row.actualLastStep = actualLastStep;
row.elapsedSeconds = elapsedOne;
row.stepsPerPlateau = getf_local(params, 'stepsPerPlateau', NaN);
row.sampleEvery = getf_local(params, 'sampleEvery', NaN);
row.nPlateaux = height(S);
row.rhoMin = min(S.meanRhoRatio, [], 'omitnan');
row.rhoMax = max(S.meanRhoRatio, [], 'omitnan');
row.massFluxLowKMaxIndex = getf_local(params, 'massFluxLowKMaxIndex', NaN);
row.cleanupStrength = getf_local(params, 'massFluxFinalVelocityProjectionStrength', NaN);
row.massFluxProjectionStrength = getf_local(params, 'massFluxProjectionStrength', NaN);

idxFit = find(string(F.branch) == "compression", 1, 'first');
if isempty(idxFit)
    idxFit = find(string(F.branch) == "all_plateaux", 1, 'first');
end
if ~isempty(idxFit)
    row.KeffKinetic = F.KeffKinetic(idxFit);
    row.KeffWall = F.KeffWall(idxFit);
    row.KeffExcess = F.KeffExcess(idxFit);
    row.KeffExtraFit = F.KeffExcess(idxFit);
    row.KeffExtraDifference = F.KeffWall(idxFit) - F.KeffKinetic(idxFit);
    row.R2Pkin = F.R2Pkin(idxFit);
    row.R2Pwall = F.R2Pwall(idxFit);
    row.R2Pexcess = F.R2Pexcess(idxFit);
end

[~, imax] = max(S.meanRhoRatio);
if ~isempty(imax) && isfinite(imax)
    row.PkinRelMax = S.PkinRel(imax);
    row.PwallRelMax = S.PwallRel(imax);
    row.PexcessRelMax = S.PexcessRel(imax);
    row.stderrPkinMax = S.stderrPkin(imax);
    row.stderrPwallMax = S.stderrPwall(imax);
    row.stderrPexcessMax = S.stderrPexcess(imax);
end

row.meanStdN = mean(S.meanStdN, 'omitnan');
row.meanLowKProxy = mean(S.meanLowKProxy, 'omitnan');
row.meanDivUAfter = mean(S.meanDivUAfter, 'omitnan');
row.meanMassFluxResidual = mean(S.meanMassFluxResidual, 'omitnan');
row.meanMassFluxAfter = mean(S.meanMassFluxAfter, 'omitnan');
row.meanKBTCell = mean(S.meanKBTCell, 'omitnan');
row.meanWallHitsPerPlateau = mean(S.totalWallHitsTop, 'omitnan');
end

function make_beta_sweep_figures(T, outputDir)
mask = string(T.status) == "ok" & isfinite(T.beta);
if ~any(mask)
    return;
end
T = T(mask,:);
[~, ord] = sort(T.beta);
T = T(ord,:);

fig1 = figure('Name', 'Q9 piston EOS beta sweep: moduli');
plot(T.beta, T.KeffKinetic, 'o-', 'DisplayName', 'K kinetic'); hold on;
plot(T.beta, T.KeffWall, 'o-', 'DisplayName', 'K wall');
plot(T.beta, T.KeffExtraDifference, 'o-', 'DisplayName', 'K wall - K kinetic');
grid on; xlabel('\beta = massFluxDensityRelaxationBeta'); ylabel('K_{eff}');
title('Effective compressibility vs Q9 beta'); legend('Location','best');
saveas(fig1, fullfile(outputDir, 'q9_piston_eos_beta_sweep_keff.png'));

fig2 = figure('Name', 'Q9 piston EOS beta sweep: fit quality');
plot(T.beta, T.R2Pwall, 'o-', 'DisplayName', 'R^2 wall'); hold on;
plot(T.beta, T.R2Pkin, 'o-', 'DisplayName', 'R^2 kinetic');
plot(T.beta, T.R2Pexcess, 'o-', 'DisplayName', 'R^2 excess');
grid on; xlabel('\beta = massFluxDensityRelaxationBeta'); ylabel('R^2'); ylim([0 1.05]);
title('EOS fit quality vs Q9 beta'); legend('Location','best');
saveas(fig2, fullfile(outputDir, 'q9_piston_eos_beta_sweep_r2.png'));

fig3 = figure('Name', 'Q9 piston EOS beta sweep: diagnostics');
plot(T.beta, T.meanStdN, 'o-', 'DisplayName', 'mean std(N)'); hold on;
plot(T.beta, T.meanDivUAfter, 'o-', 'DisplayName', 'mean div after');
plot(T.beta, T.meanLowKProxy, 'o-', 'DisplayName', 'mean low-k proxy');
grid on; xlabel('\beta = massFluxDensityRelaxationBeta'); ylabel('diagnostic');
title('Density/divergence diagnostics vs Q9 beta'); legend('Location','best');
saveas(fig3, fullfile(outputDir, 'q9_piston_eos_beta_sweep_diagnostics.png'));

fig4 = figure('Name', 'Q9 piston EOS beta sweep: wall pressure at max compression');
errorbar(T.beta, T.PwallRelMax, T.stderrPwallMax, 'o-', 'DisplayName', '\Delta P wall'); hold on;
errorbar(T.beta, T.PkinRelMax, T.stderrPkinMax, 'o-', 'DisplayName', '\Delta P kinetic');
errorbar(T.beta, T.PexcessRelMax, T.stderrPexcessMax, 'o-', 'DisplayName', '\Delta P excess');
grid on; xlabel('\beta = massFluxDensityRelaxationBeta'); ylabel('\Delta P at max rho');
title('Pressure response at maximum compression'); legend('Location','best');
saveas(fig4, fullfile(outputDir, 'q9_piston_eos_beta_sweep_pressure_at_max_rho.png'));
end

function write_beta_sweep_summary(txtFile, T, betaList, rhoRatioList, stepsPerPlateau, sampleEvery, discardFraction, elapsedSweep)
fid = fopen(txtFile, 'w');
if fid < 0
    warning('Could not open sweep summary file: %s', txtFile);
    return;
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '=== Q9 piston EOS beta sweep ===\n');
fprintf(fid, 'betaList                 : %s\n', mat2str(betaList, 8));
fprintf(fid, 'rhoRatioList             : %s\n', mat2str(rhoRatioList, 8));
fprintf(fid, 'stepsPerPlateau          : %d\n', stepsPerPlateau);
fprintf(fid, 'sampleEvery              : %d\n', sampleEvery);
fprintf(fid, 'plateauDiscardFraction   : %.6g\n', discardFraction);
fprintf(fid, 'elapsed seconds          : %.3f\n', elapsedSweep);
fprintf(fid, '\n--- Per-beta summary ---\n');
for i = 1:height(T)
    fprintf(fid, 'beta %.8g [%s]\n', T.beta(i), string(T.status(i)));
    if string(T.status(i)) == "ok"
        fprintf(fid, '  K kinetic / wall / excess-fit / wall-kinetic : %.10g / %.10g / %.10g / %.10g\n', ...
            T.KeffKinetic(i), T.KeffWall(i), T.KeffExtraFit(i), T.KeffExtraDifference(i));
        fprintf(fid, '  R2 kinetic / wall / excess                  : %.6g / %.6g / %.6g\n', ...
            T.R2Pkin(i), T.R2Pwall(i), T.R2Pexcess(i));
        fprintf(fid, '  PwallRelMax +/- stderr                       : %.10g +/- %.10g\n', ...
            T.PwallRelMax(i), T.stderrPwallMax(i));
        fprintf(fid, '  mean stdN / lowK / divAfter                  : %.10g / %.10g / %.10g\n', ...
            T.meanStdN(i), T.meanLowKProxy(i), T.meanDivUAfter(i));
        fprintf(fid, '  outputDir                                    : %s\n', string(T.outputDir(i)));
    else
        fprintf(fid, '  message                                      : %s\n', string(T.message(i)));
    end
end
end

function s = beta_tag(beta)
if beta == 0
    s = '0';
else
    s = regexprep(sprintf('%.0e', beta), '[+]', '');
    s = strrep(s, '-', 'm');
end
s = regexprep(s, '[^a-zA-Z0-9_]', '_');
end

function v = getf_local(S, name, default)
if isstruct(S) && isfield(S, name)
    v = S.(name);
else
    v = default;
end
end

function S = merge_struct_local(S, U)
if ~isstruct(U)
    return;
end
fn = fieldnames(U);
for k = 1:numel(fn)
    S.(fn{k}) = U.(fn{k});
end
end

function S = rmfield_if_present(S, names)
for k = 1:numel(names)
    if isfield(S, names{k})
        S = rmfield(S, names{k});
    end
end
end
