function out = run_projection_step_channel_demo(params)
%RUN_PROJECTION_STEP_CHANNEL_DEMO Rectangular-step validation for Q6/Q9.
%
%   out = run_projection_step_channel_demo(params)
%
% This benchmark is a lighter structure-preservation test than the cylinder:
% the obstacle is rectangular and grid-aligned, so it produces shear and a
% recirculation zone without a curved under-resolved boundary.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);
params.solidMask = mpcd_step_mask(params);
params.solidMaskForProjection = params.solidMask;

rng(params.seed);
[state, initInfo] = projection_initialize_particles_step_channel(params);
Np = initInfo.Np;

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
sampleTimes = nan(nSamples, 1);
sampleSteps = nan(nSamples, 1);
meanUxSeries = nan(nSamples, 1);
meanUySeries = nan(nSamples, 1);
meanEnstrophySeries = nan(nSamples, 1);
shearOmegaRmsSeries = nan(nSamples, 1);
recirculationLengthSeries = nan(nSamples, 1);
recirculationAreaSeries = nan(nSamples, 1);
recirculationMinUxSeries = nan(nSamples, 1);
recirculationRmsUySeries = nan(nSamples, 1);
probeOmegaSeries = nan(nSamples, 1);
probeUxSeries = nan(nSamples, 1);
probeUySeries = nan(nSamples, 1);
kBTCellSeries = nan(nSamples, 1);
NStdFluidSeries = nan(nSamples, 1);
NOutBandFluidSeries = nan(nSamples, 1);
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
if params.storeVelocityMaps
    UxMaps = nan(params.Nx, params.Ny, nSamples);
    UyMaps = nan(params.Nx, params.Ny, nSamples);
else
    UxMaps = [];
    UyMaps = [];
end

% Columns are listed in out.diagColumns.
diagHistory = nan(nSamples, 20);
isamp = 0;
actualLastStep = 0;
stopRequested = false;
wallClockTic = tic;
lastProgressPrint = -inf;
lastStepDiag = [];

for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        sampleDiag = step_sample_diagnostics(state, params, it);
        sampleTimes(isamp) = it * params.dt;
        sampleSteps(isamp) = it;
        meanUxSeries(isamp) = sampleDiag.meanUxFluid;
        meanUySeries(isamp) = sampleDiag.meanUyFluid;
        meanEnstrophySeries(isamp) = sampleDiag.step.meanEnstrophy;
        shearOmegaRmsSeries(isamp) = sampleDiag.step.shearOmegaRms;
        recirculationLengthSeries(isamp) = sampleDiag.step.recirculationLength;
        recirculationAreaSeries(isamp) = sampleDiag.step.recirculationArea;
        recirculationMinUxSeries(isamp) = sampleDiag.step.recirculationMinUx;
        recirculationRmsUySeries(isamp) = sampleDiag.step.recirculationRmsUy;
        probeOmegaSeries(isamp) = sampleDiag.step.probeOmega;
        probeUxSeries(isamp) = sampleDiag.step.probeUx;
        probeUySeries(isamp) = sampleDiag.step.probeUy;
        kBTCellSeries(isamp) = sampleDiag.kBTCell;
        NStdFluidSeries(isamp) = sampleDiag.stdNFluid;
        NOutBandFluidSeries(isamp) = sampleDiag.outBandFluid;
        if params.storeDensityMaps
            NMaps(:, :, isamp) = sampleDiag.NForMetrics;
        end
        if params.storeVorticityMaps
            omegaMaps(:, :, isamp) = sampleDiag.step.omega;
        end
        if params.storeVelocityMaps
            UxMaps(:, :, isamp) = sampleDiag.G.Ux;
            UyMaps(:, :, isamp) = sampleDiag.G.Uy;
        end
        if logical(params.visualEnable) && mod(it, params.visualEvery) == 0
            projection_step_visualize_frame(state, sampleDiag, params, it);
        end
        if it == 0
            diagHistory(isamp, :) = [it, it*params.dt, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                sampleDiag.stdNFluid, sampleDiag.kBTCell, NaN, NaN, 0, 0, sampleDiag.meanUxFluid, ...
                sampleDiag.meanUyFluid, sampleDiag.step.meanEnstrophy, sampleDiag.step.recirculationLength];
        elseif ~isempty(lastStepDiag)
            diagHistory(isamp, :) = step_diag_row(it, params, lastStepDiag, sampleDiag);
        end
    end

    if logical(params.abortOnUnstable) && it > 0 && exist('sampleDiag', 'var')
        unstable = false;
        reasons = {};
        if abs(sampleDiag.meanUxFluid) > params.maxStableAbsMeanUx
            unstable = true; reasons{end+1} = 'meanUx'; %#ok<AGROW>
        end
        if abs(sampleDiag.meanUyFluid) > params.maxStableAbsMeanUy
            unstable = true; reasons{end+1} = 'meanUy'; %#ok<AGROW>
        end
        if sampleDiag.step.meanEnstrophy > params.maxStableMeanEnstrophy
            unstable = true; reasons{end+1} = 'enstrophy'; %#ok<AGROW>
        end
        if unstable
            warning('run_projection_step_channel_demo:Unstable', ...
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
    [state, stepDiag] = mpcd_step_projection_step_channel(state, stepParams);
    [state, flowControlInfo] = apply_mean_flow_control(state, params);
    stepDiag.meanFlowControl = flowControlInfo;
    stepDiag.meanVxAfterFlowControl = flowControlInfo.meanVxAfter;
    stepDiag.meanVyAfterFlowControl = flowControlInfo.meanVyAfter;
    actualLastStep = nextStep;
    lastStepDiag = stepDiag;

    if willProgress
        lastProgressPrint = actualLastStep;
        fprintf('  step-channel step %d/%d, t=%.6g, meanUx=%.4g, recircL=%.4g, enst=%.4g, div=%.3e, mfRes=%.3e, elapsed=%.1fs\n', ...
            actualLastStep, params.nSteps, actualLastStep*params.dt, ...
            finite_or_nan(getf(stepDiag, 'meanVxAfterFlowControl', getf(stepDiag, 'meanVxAfterProjection'))), ...
            sampleDiag.step.recirculationLength, sampleDiag.step.meanEnstrophy, ...
            finite_or_nan(getf(stepDiag, 'rmsDivParticleAfter')), ...
            finite_or_nan(getf(stepDiag, 'rmsMassFluxDivResidual')), toc(wallClockTic));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_step_channel_demo:TimeLimit', ...
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
shearOmegaRmsSeries = shearOmegaRmsSeries(valid);
recirculationLengthSeries = recirculationLengthSeries(valid);
recirculationAreaSeries = recirculationAreaSeries(valid);
recirculationMinUxSeries = recirculationMinUxSeries(valid);
recirculationRmsUySeries = recirculationRmsUySeries(valid);
probeOmegaSeries = probeOmegaSeries(valid);
probeUxSeries = probeUxSeries(valid);
probeUySeries = probeUySeries(valid);
kBTCellSeries = kBTCellSeries(valid);
NStdFluidSeries = NStdFluidSeries(valid);
NOutBandFluidSeries = NOutBandFluidSeries(valid);
if params.storeDensityMaps
    NMaps = NMaps(:, :, valid);
end
if params.storeVorticityMaps
    omegaMaps = omegaMaps(:, :, valid);
end
if params.storeVelocityMaps
    UxMaps = UxMaps(:, :, valid);
    UyMaps = UyMaps(:, :, valid);
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
out.shearOmegaRmsSeries = shearOmegaRmsSeries;
out.recirculationLengthSeries = recirculationLengthSeries;
out.recirculationAreaSeries = recirculationAreaSeries;
out.recirculationMinUxSeries = recirculationMinUxSeries;
out.recirculationRmsUySeries = recirculationRmsUySeries;
out.probeOmegaSeries = probeOmegaSeries;
out.probeUxSeries = probeUxSeries;
out.probeUySeries = probeUySeries;
out.kBTCellSeries = kBTCellSeries;
out.NStdFluidSeries = NStdFluidSeries;
out.NOutBandFluidSeries = NOutBandFluidSeries;
out.NMaps = NMaps;
out.omegaMaps = omegaMaps;
out.UxMaps = UxMaps;
out.UyMaps = UyMaps;
out.solidMask = params.solidMask;
out.densityMetrics = metrics;
out.diagHistory = diagHistory;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.elapsedWallClock = elapsedWallClock;
out.diagColumns = {'step','t','rmsDivBefore','rmsDivParticleAfter', ...
    'massFluxDivBefore','massFluxDivTarget','massFluxDivProjectedAfter','massFluxDivResidual', ...
    'massFluxDivParticleAfter','dvRms','popStdFluid','kBTCellAfter', ...
    'densityTransportProjectedRms','thermostatRmsDV','stepHits','stepDPx', ...
    'meanUxFluid','meanUyFluid','meanEnstrophy','recirculationLength'};
out.summary = summarize_step_run(out);

fprintf('\n=== run_projection_step_channel_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('initialPopulationMode      : %s\n', initInfo.mode);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('solid cells                : %d / %d\n', nnz(params.solidMask), numel(params.solidMask));
fprintf('steps completed            : %d / %d in %.2f s\n', out.actualLastStep, params.nSteps, out.elapsedWallClock);
fprintf('step x0/x1/height          : %.6g / %.6g / %.6g\n', params.stepX0, params.stepX1, params.stepHeight);
fprintf('stepWallMode               : %s\n', char(params.stepWallMode));
fprintf('projectionStrength         : %.6g\n', params.projectionStrength);
fprintf('massFluxProjectionMode     : %s  strength=%.6g  beta=%.6g\n', ...
    char(params.massFluxProjectionMode), params.massFluxProjectionStrength, params.massFluxDensityRelaxationBeta);
fprintf('cleanup strength           : %.6g\n', params.massFluxFinalVelocityProjectionStrength * double(params.massFluxFinalVelocityProjectionCleanup));
fprintf('mean std(N) fluid          : %.12g\n', out.summary.meanStdN);
fprintf('mean low-k energy          : %.12g\n', out.summary.meanLowKEnergy);
fprintf('mean div(u) particle after : %.12g\n', out.summary.meanRmsDivParticleAfter);
fprintf('mean mass-flux residual    : %.12g\n', out.summary.meanMassFluxDivResidual);
fprintf('mean/final enstrophy       : %.12g / %.12g\n', out.summary.meanEnstrophy, out.summary.finalEnstrophy);
fprintf('mean/final shear omega RMS : %.12g / %.12g\n', out.summary.meanShearOmegaRms, out.summary.finalShearOmegaRms);
fprintf('mean/final recirc length   : %.12g / %.12g\n', out.summary.meanRecirculationLength, out.summary.finalRecirculationLength);

if params.makeFigures
    make_figures(out);
end
end

function row = step_diag_row(it, params, stepDiag, sampleDiag)
row = [it, it*params.dt, ...
    getf(stepDiag,'rmsDivBefore'), getf(stepDiag,'rmsDivParticleAfter'), ...
    getf(stepDiag,'rmsMassFluxDivBefore'), getf(stepDiag,'rmsMassFluxDivTarget'), ...
    getf(stepDiag,'rmsMassFluxDivProjectedAfter'), getf(stepDiag,'rmsMassFluxDivResidual'), ...
    getf(stepDiag,'rmsMassFluxDivParticleAfter'), getf(stepDiag,'dvRms'), ...
    sampleDiag.stdNFluid, sampleDiag.kBTCell, ...
    getf(stepDiag,'densityTransportProjectedRms'), ...
    getf(stepDiag,'thermostatRmsVelocityChange'), ...
    getf(stepDiag,'stepHits'), getf(stepDiag,'stepDPx'), ...
    sampleDiag.meanUxFluid, sampleDiag.meanUyFluid, sampleDiag.step.meanEnstrophy, ...
    sampleDiag.step.recirculationLength];
end

function sampleDiag = step_sample_diagnostics(state, params, it)
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
solidMask = params.solidMask;
G.N(solidMask) = params.gamma;
G.rho(solidMask) = params.gamma / max(G.dx * G.dy, eps);
G.Px(solidMask) = 0;
G.Py(solidMask) = 0;
G.Ux(solidMask) = 0;
G.Uy(solidMask) = 0;
G.valid(solidMask) = true;
thermal = projection_thermal_diagnostics(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
pop = projection_population_diagnostics(state.x, params, ...
    'periodicX', true, 'periodicY', false, 'targetGamma', params.gamma, 'bandFraction', params.densityBandFraction);
step = projection_step_diagnostics(G, params);
NForMetrics = double(G.N);
NForMetrics(solidMask) = params.gamma;
fluidN = double(pop.N(~solidMask));
band = params.gamma * params.densityBandFraction;
sampleDiag = struct();
sampleDiag.it = it;
sampleDiag.t = it * params.dt;
sampleDiag.G = G;
sampleDiag.NForMetrics = NForMetrics;
sampleDiag.thermal = thermal;
sampleDiag.pop = pop;
sampleDiag.step = step;
sampleDiag.meanUxFluid = mean(G.Ux(~solidMask), 'omitnan');
sampleDiag.meanUyFluid = mean(G.Uy(~solidMask), 'omitnan');
sampleDiag.stdNFluid = std(fluidN, 'omitnan');
sampleDiag.outBandFluid = mean(abs(fluidN - params.gamma) > band, 'omitnan');
sampleDiag.kBTCell = thermal.kBTCellRelative;
end

function summary = summarize_step_run(out)
metrics = out.densityMetrics;
if isempty(metrics)
    meanStdN = mean(out.NStdFluidSeries, 'omitnan');
    meanOutBand = mean(out.NOutBandFluidSeries, 'omitnan');
    meanLowK = NaN;
    timeAvgRelRms = NaN;
else
    meanStdN = mean(out.NStdFluidSeries, 'omitnan');
    meanOutBand = mean(out.NOutBandFluidSeries, 'omitnan');
    meanLowK = metrics.summary.meanLowKEnergy;
    timeAvgRelRms = metrics.summary.timeAvgRelRms;
end

diag = out.diagHistory;
cols = out.diagColumns;
summary = struct();
summary.actualLastStep = out.actualLastStep;
summary.stoppedEarly = out.stoppedEarly;
summary.meanStdN = meanStdN;
summary.meanOutBand = meanOutBand;
summary.meanLowKEnergy = meanLowK;
summary.timeAvgRelRms = timeAvgRelRms;
summary.meanMeanUx = mean(out.meanUxSeries, 'omitnan');
summary.finalMeanUx = out.meanUxSeries(end);
summary.meanMeanUy = mean(out.meanUySeries, 'omitnan');
summary.finalMeanUy = out.meanUySeries(end);
summary.meanEnstrophy = mean(out.meanEnstrophySeries, 'omitnan');
summary.finalEnstrophy = out.meanEnstrophySeries(end);
summary.meanShearOmegaRms = mean(out.shearOmegaRmsSeries, 'omitnan');
summary.finalShearOmegaRms = out.shearOmegaRmsSeries(end);
summary.meanRecirculationLength = mean(out.recirculationLengthSeries, 'omitnan');
summary.finalRecirculationLength = out.recirculationLengthSeries(end);
summary.meanRecirculationArea = mean(out.recirculationAreaSeries, 'omitnan');
summary.finalRecirculationArea = out.recirculationAreaSeries(end);
summary.meanRecirculationMinUx = mean(out.recirculationMinUxSeries, 'omitnan');
summary.finalRecirculationMinUx = out.recirculationMinUxSeries(end);
summary.meanRecirculationRmsUy = mean(out.recirculationRmsUySeries, 'omitnan');
summary.finalRecirculationRmsUy = out.recirculationRmsUySeries(end);
summary.probeOmegaRms = sqrt(mean(out.probeOmegaSeries.^2, 'omitnan'));
summary.probeUyRms = sqrt(mean(out.probeUySeries.^2, 'omitnan'));
summary.meanKBTCell = mean(out.kBTCellSeries, 'omitnan');
summary.finalKBTCell = out.kBTCellSeries(end);
summary.meanRmsDivParticleAfter = mean_diag_col(diag, cols, 'rmsDivParticleAfter');
summary.meanMassFluxDivResidual = mean_diag_col(diag, cols, 'massFluxDivResidual');
summary.meanMassFluxDivParticleAfter = mean_diag_col(diag, cols, 'massFluxDivParticleAfter');
summary.meanDensityTransportProjectedRms = mean_diag_col(diag, cols, 'densityTransportProjectedRms');
summary.meanStepHits = mean_diag_col(diag, cols, 'stepHits');
end

function v = mean_diag_col(diag, cols, name)
id = find(strcmp(cols, name), 1);
if isempty(id) || isempty(diag)
    v = NaN;
else
    v = mean(diag(:, id), 'omitnan');
end
end

function [state, info] = apply_mean_flow_control(state, params)
mode = lower(strrep(char(string(getf(params, 'meanFlowControlMode', 'off'))), '-', '_'));
info = struct();
info.mode = mode;
info.meanVxBefore = mean(state.v(:,1), 'omitnan');
info.meanVyBefore = mean(state.v(:,2), 'omitnan');
info.dVxApplied = 0;
info.dVyApplied = 0;
if strcmp(mode, 'off')
    info.meanVxAfter = info.meanVxBefore;
    info.meanVyAfter = info.meanVyBefore;
    return;
end
targetUx = getf(params, 'targetMeanVelocityX', info.meanVxBefore);
tau = max(getf(params, 'meanFlowRelaxationTau', 0.5), eps);
maxCorr = getf(params, 'meanFlowCorrectionMax', Inf);
dvx = (targetUx - info.meanVxBefore) / tau;
dvx = max(min(dvx, maxCorr), -maxCorr);
state.v(:,1) = state.v(:,1) + dvx;
info.dVxApplied = dvx;
if logical(getf(params, 'meanFlowControlRemoveMeanY', true))
    dvy = -info.meanVyBefore / tau;
    dvy = max(min(dvy, maxCorr), -maxCorr);
    state.v(:,2) = state.v(:,2) + dvy;
    info.dVyApplied = dvy;
end
info.meanVxAfter = mean(state.v(:,1), 'omitnan');
info.meanVyAfter = mean(state.v(:,2), 'omitnan');
end

function make_figures(out)
t = out.sampleTimes;
figure('Name','Step-channel density/vorticity diagnostics','Color','w');
subplot(2,2,1);
plot(t, out.meanEnstrophySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('mean enstrophy'); grid on;
subplot(2,2,2);
plot(t, out.shearOmegaRmsSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('shear \omega RMS'); grid on;
subplot(2,2,3);
plot(t, out.recirculationLengthSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('recirculation length'); grid on;
subplot(2,2,4);
plot(t, out.probeUySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('probe Uy'); grid on;

if ~isempty(out.omegaMaps)
    figure('Name','Final step-channel vorticity','Color','w');
    imagesc(out.omegaMaps(:,:,end).'); axis image xy; colorbar;
    hold on;
    mask = out.solidMask.';
    contour(mask, [0.5 0.5], 'k', 'LineWidth', 1.0);
    title('Final omega with step mask'); xlabel('ix'); ylabel('iy');
end
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 48);
params = set_default(params, 'Ny', 24);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'nSteps', 3000);
params = set_default(params, 'sampleEvery', 50);
params = set_default(params, 'progressEvery', 1000);
params = set_default(params, 'maxWallClockSeconds', Inf);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 0.02);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.0);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'wallModeY', 'bounceback');
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'stepX0', 0.15 * params.Lx);
params = set_default(params, 'stepX1', 0.35 * params.Lx);
params = set_default(params, 'stepHeight', 0.25 * params.Ly);
params = set_default(params, 'stepWallMode', 'bounceback');
params = set_default(params, 'stepRecirculationMaxDistance', 0.50 * params.Lx);
params = set_default(params, 'stepRecirculationYMax', min(2.0 * params.stepHeight, 0.55 * params.Ly));
params = set_default(params, 'stepShearLayerMaxDistance', 0.40 * params.Lx);
params = set_default(params, 'stepShearLayerHalfHeight', max(2 * params.Ly / params.Ny, 0.08 * params.Ly));
params = set_default(params, 'stepProbeDistanceHeights', 3.0);
params = set_default(params, 'initialPopulationMode', 'exact_per_fluid_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'initialMeanVelocityX', 0.04);
params = set_default(params, 'initialMeanVelocityY', 0.0);
params = set_default(params, 'meanFlowControlMode', 'relax_to_target');
params = set_default(params, 'targetMeanVelocityX', params.initialMeanVelocityX);
params = set_default(params, 'meanFlowRelaxationTau', 0.5);
params = set_default(params, 'meanFlowCorrectionMax', 5.0e-5);
params = set_default(params, 'meanFlowControlRemoveMeanY', true);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'massFluxProjectionMode', 'off');
params = set_default(params, 'massFluxProjectionStrength', 1.0);
params = set_default(params, 'massFluxDensityRelaxationBeta', 0.0);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', true);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', false);
params = set_default(params, 'massFluxFinalVelocityProjectionStrength', 1.0);
params = set_default(params, 'massFluxTargetFilter', 'lowpass_fft');
params = set_default(params, 'massFluxLowKMaxIndex', 2);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'thermostatMinParticlesPerCell', 3);
params = set_default(params, 'densityBandFraction', 0.20);
params = set_default(params, 'densityExcludeWallCells', 0);
params = set_default(params, 'lowKMaxIndex', 2);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'storeVorticityMaps', true);
params = set_default(params, 'storeVelocityMaps', false);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'visualEnable', false);
params = set_default(params, 'visualEvery', params.sampleEvery);
params = set_default(params, 'visualFigureId', 101);
params = set_default(params, 'visualFigureName', 'Step-channel live visualization');
params = set_default(params, 'visualTitleSuffix', '');
params = set_default(params, 'visualMaxParticles', 3000);
params = set_default(params, 'visualParticleSize', 5);
params = set_default(params, 'visualDensityCLim', 0.5);
params = set_default(params, 'visualOmegaCLim', NaN);
params = set_default(params, 'visualOmegaClipSigma', 3.0);
params = set_default(params, 'visualSpeedCLim', NaN);
params = set_default(params, 'visualQuiverStrideX', 4);
params = set_default(params, 'visualQuiverStrideY', 3);
params = set_default(params, 'visualQuiverScale', 1.5);
params = set_default(params, 'visualPause', 0.0);
params = set_default(params, 'visualSaveFrames', false);
params = set_default(params, 'visualFrameDir', 'step_visual_frames');
params = set_default(params, 'visualFramePrefix', 'step_frame');
params = set_default(params, 'visualFrameResolution', 120);
params = set_default(params, 'seed', 11);
params = set_default(params, 'abortOnUnstable', true);
params = set_default(params, 'maxStableAbsMeanUx', 0.5);
params = set_default(params, 'maxStableAbsMeanUy', 0.5);
params = set_default(params, 'maxStableMeanEnstrophy', 5.0e3);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function v = getf(s, name, defaultValue)
if nargin < 3
    defaultValue = NaN;
end
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end

function y = finite_or_nan(x)
if isempty(x) || ~isnumeric(x) || ~isfinite(x)
    y = NaN;
else
    y = x;
end
end
