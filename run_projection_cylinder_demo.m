function out = run_projection_cylinder_demo(params)
%RUN_PROJECTION_CYLINDER_DEMO Cylinder/von-Karman validation for Q6/Q9.
%
%   out = run_projection_cylinder_demo(params)
%
% First Q9-native obstacle case:
%   - x-periodic channel,
%   - bounded y walls,
%   - fixed circular cylinder,
%   - body-force/mean-flow driven wake,
%   - Q6 or Q9 projection after each SRC/MPCD step.
%
% This is a diagnostic benchmark, not yet a calibrated Reynolds-number study.
% It checks whether Q9 reduces density low-k modes without erasing vorticity
% and wake fluctuations.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);
params.solidMask = mpcd_cylinder_mask(params);
params.solidMaskForProjection = params.solidMask;

rng(params.seed);
[state, initInfo] = projection_initialize_particles_cylinder(params);
Np = initInfo.Np;

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
sampleTimes = nan(nSamples, 1);
sampleSteps = nan(nSamples, 1);
meanUxSeries = nan(nSamples, 1);
meanUySeries = nan(nSamples, 1);
meanEnstrophySeries = nan(nSamples, 1);
rmsOmegaSeries = nan(nSamples, 1);
wakeOmegaSeries = nan(nSamples, 1);
wakeUxSeries = nan(nSamples, 1);
wakeUySeries = nan(nSamples, 1);
wakeOmegaUpperSeries = nan(nSamples, 1);
wakeUyUpperSeries = nan(nSamples, 1);
wakeOmegaLowerSeries = nan(nSamples, 1);
wakeUyLowerSeries = nan(nSamples, 1);
wakeOmegaAntiSymSeries = nan(nSamples, 1);
wakeUyAntiSymSeries = nan(nSamples, 1);
kBTCellSeries = nan(nSamples, 1);
UxMeanYProfiles = nan(params.Ny, nSamples);
UyMeanYProfiles = nan(params.Ny, nSamples);
NMeanYProfiles = nan(params.Ny, nSamples);
if params.storeDensityMaps
    NMaps = nan(params.Nx, params.Ny, nSamples);
else
    NMaps = [];
end
if params.storeVorticityMaps
    omegaMaps = nan(params.Nx, params.Ny, nSamples);
else
    omegaMaps = [];
end

% Columns are listed in out.diagColumns.
diagHistory = nan(nSamples, 22);

isamp = 0;
actualLastStep = 0;
stopRequested = false;
wallClockTic = tic;
lastProgressPrint = -inf;
lastStepDiag = [];

for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        sampleDiag = cylinder_sample_diagnostics(state, params, it);

        sampleTimes(isamp) = it * params.dt;
        sampleSteps(isamp) = it;
        meanUxSeries(isamp) = sampleDiag.meanUxFluid;
        meanUySeries(isamp) = sampleDiag.meanUyFluid;
        meanEnstrophySeries(isamp) = sampleDiag.vorticity.meanEnstrophy;
        rmsOmegaSeries(isamp) = sampleDiag.vorticity.rmsOmega;
        wakeOmegaSeries(isamp) = sampleDiag.vorticity.wakeOmega;
        wakeUxSeries(isamp) = sampleDiag.vorticity.wakeUx;
        wakeUySeries(isamp) = sampleDiag.vorticity.wakeUy;
        wakeOmegaUpperSeries(isamp) = sampleDiag.vorticity.wakeOmegaUpper;
        wakeUyUpperSeries(isamp) = sampleDiag.vorticity.wakeUyUpper;
        wakeOmegaLowerSeries(isamp) = sampleDiag.vorticity.wakeOmegaLower;
        wakeUyLowerSeries(isamp) = sampleDiag.vorticity.wakeUyLower;
        wakeOmegaAntiSymSeries(isamp) = sampleDiag.vorticity.wakeOmegaAntiSym;
        wakeUyAntiSymSeries(isamp) = sampleDiag.vorticity.wakeUyAntiSym;
        kBTCellSeries(isamp) = sampleDiag.kBTCell;
        UxMeanYProfiles(:, isamp) = mean(sampleDiag.G.Ux, 1, 'omitnan').';
        UyMeanYProfiles(:, isamp) = mean(sampleDiag.G.Uy, 1, 'omitnan').';
        NMeanYProfiles(:, isamp) = mean(sampleDiag.NForMetrics, 1, 'omitnan').';
        if params.storeDensityMaps
            NMaps(:, :, isamp) = sampleDiag.NForMetrics;
        end
        if params.storeVorticityMaps
            omegaMaps(:, :, isamp) = sampleDiag.vorticity.omega;
        end

        if it == 0
            diagHistory(isamp, :) = [it, it*params.dt, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                sampleDiag.stdNFluid, sampleDiag.stdNFluid, sampleDiag.kBTCell, NaN, NaN, ...
                0, 0, sampleDiag.meanUxFluid, sampleDiag.meanUyFluid, sampleDiag.vorticity.meanEnstrophy, ...
                sampleDiag.vorticity.wakeOmega, NaN];
        elseif ~isempty(lastStepDiag)
            diagHistory(isamp, :) = step_diag_row(it, params, lastStepDiag, sampleDiag);
        end
    end

    if logical(getf(params, 'abortOnUnstable', false)) && it > 0 && exist('sampleDiag', 'var')
        unstable = false;
        reasons = {};
        if abs(sampleDiag.meanUxFluid) > getf(params, 'maxStableAbsMeanUx', Inf)
            unstable = true; reasons{end+1} = 'meanUx'; %#ok<AGROW>
        end
        if sampleDiag.vorticity.meanEnstrophy > getf(params, 'maxStableMeanEnstrophy', Inf)
            unstable = true; reasons{end+1} = 'enstrophy'; %#ok<AGROW>
        end
        if abs(sampleDiag.meanUyFluid) > getf(params, 'maxStableAbsMeanUy', Inf)
            unstable = true; reasons{end+1} = 'meanUy'; %#ok<AGROW>
        end
        if unstable
            warning('run_projection_cylinder_demo:Unstable', ...
                'Stopping early at step %d/%d because stability thresholds were exceeded: %s.', ...
                it, params.nSteps, strjoin(reasons, ', '));
            actualLastStep = it;
            stopRequested = true;
        end
    end
    if stopRequested
        break;
    end

    if it == params.nSteps
        break;
    end

    nextStep = it + 1;
    willSample = mod(nextStep, params.sampleEvery) == 0;
    willProgress = params.progressEvery > 0 && mod(nextStep, params.progressEvery) == 0 && nextStep ~= lastProgressPrint;
    willBeFinal = nextStep == params.nSteps;

    stepParams = params;
    stepParams.computeDiagnostics = params.computeFullDiagnosticsEveryStep || willSample || willProgress || willBeFinal;

    [state, stepDiag] = mpcd_step_projection_cylinder(state, stepParams);
    [state, flowControlInfo] = apply_mean_flow_control(state, params);
    stepDiag.meanFlowControl = flowControlInfo;
    stepDiag.meanVxAfterFlowControl = flowControlInfo.meanVxAfter;
    stepDiag.meanVyAfterFlowControl = flowControlInfo.meanVyAfter;
    actualLastStep = nextStep;
    lastStepDiag = stepDiag;

    if willProgress
        lastProgressPrint = actualLastStep;
        fprintf('  cylinder step %d/%d, t=%.6g, meanUx=%.4g, hits=%d, div=%.3e, mfRes=%.3e, elapsed=%.1fs\n', ...
            actualLastStep, params.nSteps, actualLastStep*params.dt, ...
            finite_or_nan(getf(stepDiag, 'meanVxAfterFlowControl', getf(stepDiag, 'meanVxAfterProjection'))), ...
            round(finite_or_zero(getf(stepDiag, 'cylinderHits'))), ...
            finite_or_nan(getf(stepDiag, 'rmsDivParticleAfter')), ...
            finite_or_nan(getf(stepDiag, 'rmsMassFluxDivResidual')), toc(wallClockTic));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_cylinder_demo:TimeLimit', ...
            'Stopping early at step %d/%d after %.1f s because maxWallClockSeconds=%.1f.', ...
            actualLastStep, params.nSteps, toc(wallClockTic), params.maxWallClockSeconds);
        stopRequested = true;
    end
    if stopRequested
        break;
    end
end

if actualLastStep == 0
    actualLastStep = min(params.nSteps, max(sampleSteps(isfinite(sampleSteps))));
end
elapsedWallClock = toc(wallClockTic);

valid = isfinite(sampleTimes);
sampleTimes = sampleTimes(valid);
sampleSteps = sampleSteps(valid);
meanUxSeries = meanUxSeries(valid);
meanUySeries = meanUySeries(valid);
meanEnstrophySeries = meanEnstrophySeries(valid);
rmsOmegaSeries = rmsOmegaSeries(valid);
wakeOmegaSeries = wakeOmegaSeries(valid);
wakeUxSeries = wakeUxSeries(valid);
wakeUySeries = wakeUySeries(valid);
wakeOmegaUpperSeries = wakeOmegaUpperSeries(valid);
wakeUyUpperSeries = wakeUyUpperSeries(valid);
wakeOmegaLowerSeries = wakeOmegaLowerSeries(valid);
wakeUyLowerSeries = wakeUyLowerSeries(valid);
wakeOmegaAntiSymSeries = wakeOmegaAntiSymSeries(valid);
wakeUyAntiSymSeries = wakeUyAntiSymSeries(valid);
kBTCellSeries = kBTCellSeries(valid);
UxMeanYProfiles = UxMeanYProfiles(:, valid);
UyMeanYProfiles = UyMeanYProfiles(:, valid);
NMeanYProfiles = NMeanYProfiles(:, valid);
if params.storeDensityMaps
    NMaps = NMaps(:, :, valid);
end
if params.storeVorticityMaps
    omegaMaps = omegaMaps(:, :, valid);
end
diagHistory = diagHistory(isfinite(diagHistory(:, 1)), :);

metrics = [];
if params.storeDensityMaps && ~isempty(NMaps)
    metrics = projection_density_homogeneity_metrics(NMaps, params, ...
        'sampleTimes', sampleTimes, ...
        'targetGamma', params.gamma, ...
        'bandFraction', params.densityBandFraction, ...
        'excludeWallCells', params.densityExcludeWallCells, ...
        'lowKMaxIndex', params.lowKMaxIndex);
end

out = struct();
out.params = params;
out.initialization = initInfo;
out.state = state;
out.sampleTimes = sampleTimes;
out.sampleSteps = sampleSteps;
out.meanUxSeries = meanUxSeries;
out.meanUySeries = meanUySeries;
out.meanEnstrophySeries = meanEnstrophySeries;
out.rmsOmegaSeries = rmsOmegaSeries;
out.wakeOmegaSeries = wakeOmegaSeries;
out.wakeUxSeries = wakeUxSeries;
out.wakeUySeries = wakeUySeries;
out.wakeOmegaUpperSeries = wakeOmegaUpperSeries;
out.wakeUyUpperSeries = wakeUyUpperSeries;
out.wakeOmegaLowerSeries = wakeOmegaLowerSeries;
out.wakeUyLowerSeries = wakeUyLowerSeries;
out.wakeOmegaAntiSymSeries = wakeOmegaAntiSymSeries;
out.wakeUyAntiSymSeries = wakeUyAntiSymSeries;
out.kBTCellSeries = kBTCellSeries;
out.UxMeanYProfiles = UxMeanYProfiles;
out.UyMeanYProfiles = UyMeanYProfiles;
out.NMeanYProfiles = NMeanYProfiles;
out.NMaps = NMaps;
out.omegaMaps = omegaMaps;
out.solidMask = params.solidMask;
out.densityMetrics = metrics;
out.diagHistory = diagHistory;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.elapsedWallClock = elapsedWallClock;
out.diagColumns = {'step','t','rmsDivBefore','rmsDivParticleAfter', ...
    'massFluxDivBefore','massFluxDivTarget','massFluxDivProjectedAfter','massFluxDivResidual', ...
    'massFluxDivParticleAfter','dvRms','popStdClassic','popStdProjection','kBTCellAfter', ...
    'densityTransportProjectedRms','thermostatRmsDV','cylinderHits','cylinderDPx', ...
    'meanUxFluid','meanUyFluid','meanEnstrophy','wakeOmega','wakeUy'};
out.shedding = compute_shedding_diagnostics(out);
out.summary = summarize_cylinder_run(out);

fprintf('\n=== run_projection_cylinder_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('initialPopulationMode      : %s\n', initInfo.mode);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('solid cells                : %d / %d\n', nnz(params.solidMask), numel(params.solidMask));
fprintf('steps completed            : %d / %d in %.2f s\n', out.actualLastStep, params.nSteps, out.elapsedWallClock);
fprintf('cylinder R, center         : %.6g at (%.6g, %.6g)\n', params.cylinderRadius, params.cylinderCenterX, params.cylinderCenterY);
fprintf('cylinderWallMode           : %s\n', char(params.cylinderWallMode));
fprintf('projectionStrength         : %.6g\n', params.projectionStrength);
fprintf('massFluxProjectionMode     : %s  strength=%.6g  beta=%.6g\n', ...
    char(params.massFluxProjectionMode), params.massFluxProjectionStrength, params.massFluxDensityRelaxationBeta);
fprintf('cleanup strength           : %.6g\n', params.massFluxFinalVelocityProjectionStrength * double(params.massFluxFinalVelocityProjectionCleanup));
fprintf('mean std(N) fluid          : %.12g\n', out.summary.meanStdN);
fprintf('mean low-k energy          : %.12g\n', out.summary.meanLowKEnergy);
fprintf('mean div(u) particle after : %.12g\n', out.summary.meanRmsDivParticleAfter);
fprintf('mean mass-flux residual    : %.12g\n', out.summary.meanMassFluxDivResidual);
fprintf('mean enstrophy             : %.12g\n', out.summary.meanEnstrophy);
fprintf('wake omega RMS             : %.12g\n', out.summary.wakeOmegaRms);
fprintf('wake offcenter St(Umean)   : %.12g  f=%.12g  SNR=%.12g\n', ...
    out.summary.strouhalWakeUyUpper, out.summary.sheddingFrequencyWakeUyUpper, out.summary.sheddingSNRWakeUyUpper);
fprintf('estimated Re(Umean,nu)     : %.12g\n', out.summary.reEffMeanUx);

if params.makeFigures
    make_figures(out);
end
end


function [state, info] = apply_mean_flow_control(state, params)
mode = lower(strrep(char(string(getf(params, 'meanFlowControlMode', 'off'))), '-', '_'));
info = struct('enabled', false, 'mode', mode, 'targetUx', NaN, 'targetUy', NaN, ...
    'meanVxBefore', mean(state.v(:,1), 'omitnan'), 'meanVyBefore', mean(state.v(:,2), 'omitnan'), ...
    'meanVxAfter', mean(state.v(:,1), 'omitnan'), 'meanVyAfter', mean(state.v(:,2), 'omitnan'), ...
    'deltaUx', 0.0, 'deltaUy', 0.0, 'gain', 0.0);
if strcmp(mode, 'off') || strcmp(mode, 'none')
    return;
end
if ~ismember(mode, {'relax_to_target','target','bulk_velocity_control'})
    error('Unknown meanFlowControlMode: %s', mode);
end
Ux0 = info.meanVxBefore;
Uy0 = info.meanVyBefore;
targetUx = getf(params, 'targetMeanVelocityX', getf(params, 'initialMeanVelocityX', 0.0));
targetUy = getf(params, 'targetMeanVelocityY', 0.0);
tau = max(getf(params, 'meanFlowRelaxationTau', 0.05), eps);
gain = min(max(params.dt / tau, 0.0), 1.0);
maxDelta = getf(params, 'meanFlowCorrectionMax', Inf);

dUx = gain * (targetUx - Ux0);
if isfinite(maxDelta)
    dUx = min(max(dUx, -maxDelta), maxDelta);
end
state.v(:,1) = state.v(:,1) + dUx;

dUy = 0.0;
if logical(getf(params, 'meanFlowControlRemoveMeanY', true))
    dUy = -Uy0;
else
    dUy = gain * (targetUy - Uy0);
end
if isfinite(maxDelta)
    dUy = min(max(dUy, -maxDelta), maxDelta);
end
state.v(:,2) = state.v(:,2) + dUy;

info.enabled = true;
info.targetUx = targetUx;
info.targetUy = targetUy;
info.deltaUx = dUx;
info.deltaUy = dUy;
info.gain = gain;
info.meanVxAfter = mean(state.v(:,1), 'omitnan');
info.meanVyAfter = mean(state.v(:,2), 'omitnan');
end

function row = step_diag_row(it, params, stepDiag, sampleDiag)
row = [it, it*params.dt, ...
    getf(stepDiag,'rmsDivBefore'), getf(stepDiag,'rmsDivParticleAfter'), ...
    getf(stepDiag,'rmsMassFluxDivBefore'), getf(stepDiag,'rmsMassFluxDivTarget'), ...
    getf(stepDiag,'rmsMassFluxDivProjectedAfter'), getf(stepDiag,'rmsMassFluxDivResidual'), ...
    getf(stepDiag,'rmsMassFluxDivParticleAfter'), getf(stepDiag,'dvRms'), ...
    getf(getf(stepDiag,'populationAfterClassic',struct()),'stdN'), ...
    getf(getf(stepDiag,'populationAfterProjection',struct()),'stdN'), ...
    getf(stepDiag,'kBTCellAfterProjection'), getf(stepDiag,'densityTransportProjectedRms'), ...
    getf(stepDiag,'thermostatRmsVelocityChange'), getf(stepDiag,'cylinderHits'), ...
    getf(stepDiag,'cylinderDPx'), sampleDiag.meanUxFluid, sampleDiag.meanUyFluid, ...
    sampleDiag.vorticity.meanEnstrophy, sampleDiag.vorticity.wakeOmega, sampleDiag.vorticity.wakeUy];
end

function sampleDiag = cylinder_sample_diagnostics(state, params, it)
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
thermal = projection_thermal_diagnostics(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
solidMask = params.solidMask;
fluidMask = ~solidMask;
NForMetrics = double(G.N);
NForMetrics(solidMask) = params.gamma;
GforVorticity = G;
GforVorticity.Ux(solidMask) = 0;
GforVorticity.Uy(solidMask) = 0;
vort = projection_vorticity_diagnostics(GforVorticity, params);

sampleDiag = struct();
sampleDiag.it = it;
sampleDiag.t = it * params.dt;
sampleDiag.G = G;
sampleDiag.NForMetrics = NForMetrics;
sampleDiag.vorticity = vort;
sampleDiag.meanUxFluid = mean(G.Ux(fluidMask), 'omitnan');
sampleDiag.meanUyFluid = mean(G.Uy(fluidMask), 'omitnan');
sampleDiag.stdNFluid = std(double(G.N(fluidMask)), 'omitnan');
sampleDiag.outBandFluid = out_band_fraction(double(G.N(fluidMask)), params.gamma, params.densityBandFraction);
sampleDiag.kBTCell = thermal.kBTCellRelative;
end

function summary = summarize_cylinder_run(out)
metrics = out.densityMetrics;
if isempty(metrics)
    meanStdN = NaN;
    meanOutBand = NaN;
    meanLowK = NaN;
    timeAvgRelRms = NaN;
else
    meanStdN = metrics.summary.meanStdN;
    meanOutBand = metrics.summary.meanOutBandFraction;
    meanLowK = metrics.summary.meanLowKEnergy;
    timeAvgRelRms = metrics.summary.timeAvgRelRms;
end

diag = out.diagHistory;
cols = out.diagColumns;
summary = struct();
summary.actualLastStep = out.actualLastStep;
summary.stoppedEarly = out.stoppedEarly;
summary.finalMeanUx = out.meanUxSeries(end);
summary.finalMeanUy = out.meanUySeries(end);
summary.meanMeanUx = mean(out.meanUxSeries, 'omitnan');
summary.meanMeanUy = mean(out.meanUySeries, 'omitnan');
summary.meanStdN = meanStdN;
summary.meanOutBand = meanOutBand;
summary.meanLowKEnergy = meanLowK;
summary.timeAvgRelRms = timeAvgRelRms;
summary.meanEnstrophy = mean(out.meanEnstrophySeries, 'omitnan');
summary.finalEnstrophy = out.meanEnstrophySeries(end);
summary.rmsOmegaMean = mean(out.rmsOmegaSeries, 'omitnan');
summary.wakeOmegaRms = sqrt(mean((out.wakeOmegaSeries - mean(out.wakeOmegaSeries,'omitnan')).^2, 'omitnan'));
summary.wakeUyRms = sqrt(mean((out.wakeUySeries - mean(out.wakeUySeries,'omitnan')).^2, 'omitnan'));
summary.wakeOmegaUpperRms = rms_demeaned(out.wakeOmegaUpperSeries);
summary.wakeUyUpperRms = rms_demeaned(out.wakeUyUpperSeries);
summary.wakeOmegaLowerRms = rms_demeaned(out.wakeOmegaLowerSeries);
summary.wakeUyLowerRms = rms_demeaned(out.wakeUyLowerSeries);
summary.wakeOmegaAntiSymRms = rms_demeaned(out.wakeOmegaAntiSymSeries);
summary.wakeUyAntiSymRms = rms_demeaned(out.wakeUyAntiSymSeries);
summary.wakeOmegaMeanAbs = mean(abs(out.wakeOmegaSeries), 'omitnan');
summary.wakeOmegaSignChanges = count_sign_changes(out.wakeOmegaSeries);
summary.wakeUyUpperSignChanges = count_sign_changes(out.wakeUyUpperSeries);
summary.wakeOmegaAntiSymSignChanges = count_sign_changes(out.wakeOmegaAntiSymSeries);
summary.kBTCellMean = mean(out.kBTCellSeries, 'omitnan');
summary.kBTCellFinal = out.kBTCellSeries(end);
summary.meanCylinderHits = mean_diag_col(diag, cols, 'cylinderHits');
summary.meanCylinderDPx = mean_diag_col(diag, cols, 'cylinderDPx');
summary.meanRmsDivParticleAfter = mean_diag_col(diag, cols, 'rmsDivParticleAfter');
summary.meanMassFluxDivResidual = mean_diag_col(diag, cols, 'massFluxDivResidual');
summary.meanMassFluxDivParticleAfter = mean_diag_col(diag, cols, 'massFluxDivParticleAfter');
summary.meanDensityTransportProjectedRms = mean_diag_col(diag, cols, 'densityTransportProjectedRms');
summary.sheddingFrequencyWakeUyUpper = get_shed_field(out.shedding.wakeUyUpper, 'frequency');
summary.strouhalWakeUyUpper = get_shed_field(out.shedding.wakeUyUpper, 'strouhal');
summary.sheddingSNRWakeUyUpper = get_shed_field(out.shedding.wakeUyUpper, 'spectralSNR');
summary.sheddingFrequencyWakeOmegaAntiSym = get_shed_field(out.shedding.wakeOmegaAntiSym, 'frequency');
summary.strouhalWakeOmegaAntiSym = get_shed_field(out.shedding.wakeOmegaAntiSym, 'strouhal');
summary.sheddingSNROmegaAntiSym = get_shed_field(out.shedding.wakeOmegaAntiSym, 'spectralSNR');
summary.reEffMeanUx = estimate_reynolds(abs(summary.meanMeanUx), out.params);
summary.reEffFinalUx = estimate_reynolds(abs(summary.finalMeanUx), out.params);
summary.targetStrouhal = getf(out.params, 'targetStrouhal', 0.18);
summary.targetRe = getf(out.params, 'targetRe', NaN);
end

function shedding = compute_shedding_diagnostics(out)
Uref = abs(mean(out.meanUxSeries, 'omitnan'));
shedding = struct();
shedding.wakeUyUpper = projection_wake_shedding_diagnostics(out.sampleTimes, out.wakeUyUpperSeries, out.params, Uref);
shedding.wakeOmegaUpper = projection_wake_shedding_diagnostics(out.sampleTimes, out.wakeOmegaUpperSeries, out.params, Uref);
shedding.wakeUyAntiSym = projection_wake_shedding_diagnostics(out.sampleTimes, out.wakeUyAntiSymSeries, out.params, Uref);
shedding.wakeOmegaAntiSym = projection_wake_shedding_diagnostics(out.sampleTimes, out.wakeOmegaAntiSymSeries, out.params, Uref);
end

function Re = estimate_reynolds(U, params)
nu = getf(params, 'nuEffForRe', NaN);
D = 2 * getf(params, 'cylinderRadius', NaN);
if isfinite(U) && isfinite(nu) && nu > 0 && isfinite(D) && D > 0
    Re = U * D / nu;
else
    Re = NaN;
end
end

function x = get_shed_field(shed, name)
if isstruct(shed) && isfield(shed, name) && ~isempty(shed.(name))
    x = shed.(name);
else
    x = NaN;
end
end

function y = rms_demeaned(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    x = x - mean(x, 'omitnan');
    y = sqrt(mean(x.^2, 'omitnan'));
end
end

function v = mean_diag_col(diag, cols, name)
id = find(strcmp(cols, name), 1);
if isempty(id) || isempty(diag)
    v = NaN;
else
    v = mean(diag(:, id), 'omitnan');
end
end

function n = count_sign_changes(x)
x = x(:);
x = x(isfinite(x));
if numel(x) < 2
    n = 0;
    return;
end
x = x - mean(x, 'omitnan');
s = sign(x);
s(s == 0) = [];
if numel(s) < 2
    n = 0;
else
    n = nnz(s(2:end) ~= s(1:end-1));
end
end

function f = out_band_fraction(n, gamma, bandFraction)
lo = gamma * (1 - bandFraction);
hi = gamma * (1 + bandFraction);
f = nnz(n < lo | n > hi) / numel(n);
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 48);
params = set_default(params, 'Ny', 24);
params = set_default(params, 'gamma', 10);
params = set_default(params, 'nSteps', 5000);
params = set_default(params, 'sampleEvery', 50);
params = set_default(params, 'progressEvery', 500);
params = set_default(params, 'maxWallClockSeconds', Inf);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 0.2);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.02);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'wallModeY', 'thermalize');
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'initialMeanVelocityX', 0.08);
params = set_default(params, 'initialMeanVelocityY', 0.0);
params = set_default(params, 'seed', 11);
params = set_default(params, 'cylinderCenterX', 0.35 * params.Lx);
params = set_default(params, 'cylinderCenterY', 0.50 * params.Ly);
params = set_default(params, 'cylinderRadius', 0.10 * params.Ly);
params = set_default(params, 'cylinderWallMode', 'bounceback');
params = set_default(params, 'wakeProbeDistanceDiameters', 2.0);
params = set_default(params, 'wakeProbeYOffsetDiameters', 0.35);
params = set_default(params, 'sheddingTransientFraction', 0.40);
params = set_default(params, 'targetStrouhal', 0.18);
params = set_default(params, 'targetRe', NaN);
params = set_default(params, 'nuEffForRe', NaN);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'projectedStrength', 1.0);
params = set_default(params, 'massFluxProjectionMode', 'off');
params = set_default(params, 'massFluxProjectionStrength', 1.0);
params = set_default(params, 'massFluxDensityRelaxationBeta', 0.0);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', false);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', false);
params = set_default(params, 'massFluxFinalVelocityProjectionStrength', 1.0);
params = set_default(params, 'massFluxTargetFilter', 'lowpass_fft');
params = set_default(params, 'massFluxLowKMaxIndex', 2);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'storeVorticityMaps', true);
params = set_default(params, 'densityBandFraction', 0.20);
params = set_default(params, 'densityExcludeWallCells', 0);
params = set_default(params, 'lowKMaxIndex', 2);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'computeDiagnostics', true);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);

% Mean-flow control is useful for long periodic-cylinder runs.  A constant
% body force can inject energy indefinitely in this closed periodic domain;
% this optional controller keeps the bulk velocity near a prescribed target
% without changing particle positions or thermal fluctuations.
params = set_default(params, 'meanFlowControlMode', 'off');
params = set_default(params, 'targetMeanVelocityX', params.initialMeanVelocityX);
params = set_default(params, 'targetMeanVelocityY', 0.0);
params = set_default(params, 'meanFlowRelaxationTau', 0.05);
params = set_default(params, 'meanFlowCorrectionMax', Inf);
params = set_default(params, 'meanFlowControlRemoveMeanY', true);

% Safety stops for exploratory long runs.  These are disabled by default and
% enabled by the sweep script.
params = set_default(params, 'abortOnUnstable', false);
params = set_default(params, 'maxStableAbsMeanUx', Inf);
params = set_default(params, 'maxStableAbsMeanUy', Inf);
params = set_default(params, 'maxStableMeanEnstrophy', Inf);

% Band-limited shedding search.  Enabled by default because unrestricted FFT
% peak picking tends to select collision-scale peaks in noisy particle runs.
params = set_default(params, 'sheddingUseStrouhalBand', true);
params = set_default(params, 'sheddingStrouhalMin', 0.05);
params = set_default(params, 'sheddingStrouhalMax', 0.50);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
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

function value = finite_or_nan(value)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    value = NaN;
end
end

function value = finite_or_zero(value)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    value = 0;
end
end

function make_figures(out)
figure('Name', 'Cylinder density and vorticity diagnostics');
subplot(2,1,1);
plot(out.sampleTimes, out.meanEnstrophySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('mean enstrophy'); grid on;
subplot(2,1,2);
plot(out.sampleTimes, out.wakeOmegaSeries, 'LineWidth', 1.2);
hold on;
plot(out.sampleTimes, out.wakeOmegaUpperSeries, 'LineWidth', 1.0);
plot(out.sampleTimes, out.wakeOmegaAntiSymSeries, 'LineWidth', 1.0);
xlabel('t'); ylabel('wake \omega'); grid on;
legend({'center','upper','antisym'}, 'Location', 'best');
end
