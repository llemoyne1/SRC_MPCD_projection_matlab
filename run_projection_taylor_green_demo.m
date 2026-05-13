function out = run_projection_taylor_green_demo(params)
%RUN_PROJECTION_TAYLOR_GREEN_DEMO Taylor-Green structure test for Q6/Q9.
%
%   out = run_projection_taylor_green_demo(params)
%
% Fully periodic benchmark designed to test whether Q9 preserves coherent
% vortical structures while reducing low-k density modes.  No walls and no
% obstacle are present, so the test isolates projection effects from boundary
% condition errors.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);
rng(params.seed);
[state, initInfo] = projection_initialize_particles_taylor_green(params);
Np = initInfo.Np;

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
sampleTimes = nan(nSamples, 1);
sampleSteps = nan(nSamples, 1);
meanUxSeries = nan(nSamples, 1);
meanUySeries = nan(nSamples, 1);
kBTCellSeries = nan(nSamples, 1);
tgAmplitudeSeries = nan(nSamples, 1);
tgModeEnergySeries = nan(nSamples, 1);
tgCoherenceSeries = nan(nSamples, 1);
tgResidualEnergySeries = nan(nSamples, 1);
tgTotalHydroEnergySeries = nan(nSamples, 1);
enstrophySeries = nan(nSamples, 1);
omegaRmsSeries = nan(nSamples, 1);
highKVelocityFractionSeries = nan(nSamples, 1);
lowKVelocityFractionSeries = nan(nSamples, 1);
tgModeVelocityFractionSeries = nan(nSamples, 1);
rmsDivGridSeries = nan(nSamples, 1);
NStdSeries = nan(nSamples, 1);
NOutBandSeries = nan(nSamples, 1);
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
diagHistory = nan(nSamples, 18);
isamp = 0;
actualLastStep = 0;
stopRequested = false;
wallClockTic = tic;
lastProgressPrint = -inf;
lastStepDiag = [];

for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        sampleDiag = tg_sample_diagnostics(state, params, it);
        sampleTimes(isamp) = sampleDiag.t;
        sampleSteps(isamp) = it;
        meanUxSeries(isamp) = sampleDiag.meanUx;
        meanUySeries(isamp) = sampleDiag.meanUy;
        kBTCellSeries(isamp) = sampleDiag.kBTCell;
        tgAmplitudeSeries(isamp) = sampleDiag.tg.modeAmplitude;
        tgModeEnergySeries(isamp) = sampleDiag.tg.modeEnergy;
        tgCoherenceSeries(isamp) = sampleDiag.tg.modeCoherence;
        tgResidualEnergySeries(isamp) = sampleDiag.tg.residualEnergy;
        tgTotalHydroEnergySeries(isamp) = sampleDiag.tg.totalHydroEnergy;
        enstrophySeries(isamp) = sampleDiag.tg.enstrophy;
        omegaRmsSeries(isamp) = sampleDiag.tg.omegaRms;
        highKVelocityFractionSeries(isamp) = sampleDiag.tg.highKEnergyFractionVelocity;
        lowKVelocityFractionSeries(isamp) = sampleDiag.tg.lowKEnergyFractionVelocity;
        tgModeVelocityFractionSeries(isamp) = sampleDiag.tg.tgModeEnergyFractionVelocity;
        rmsDivGridSeries(isamp) = sampleDiag.tg.rmsDiv;
        NStdSeries(isamp) = sampleDiag.pop.stdN;
        NOutBandSeries(isamp) = sampleDiag.pop.outBandFraction;
        if params.storeDensityMaps
            NMaps(:, :, isamp) = sampleDiag.G.N;
        end
        if params.storeVorticityMaps
            omegaMaps(:, :, isamp) = sampleDiag.tg.omega;
        end
        if it == 0
            diagHistory(isamp, :) = [it, it*params.dt, NaN, NaN, NaN, NaN, NaN, NaN, ...
                NaN, NaN, sampleDiag.pop.stdN, sampleDiag.kBTCell, sampleDiag.tg.modeAmplitude, ...
                sampleDiag.tg.modeCoherence, sampleDiag.tg.enstrophy, sampleDiag.tg.highKEnergyFractionVelocity, ...
                sampleDiag.meanUx, sampleDiag.meanUy];
        elseif ~isempty(lastStepDiag)
            diagHistory(isamp, :) = step_diag_row(it, params, lastStepDiag, sampleDiag);
        end
    end

    if logical(params.abortOnUnstable) && it > 0 && exist('sampleDiag', 'var')
        unstable = false;
        reasons = {};
        if abs(sampleDiag.meanUx) > params.maxStableAbsMeanUx
            unstable = true; reasons{end+1} = 'meanUx'; %#ok<AGROW>
        end
        if abs(sampleDiag.meanUy) > params.maxStableAbsMeanUy
            unstable = true; reasons{end+1} = 'meanUy'; %#ok<AGROW>
        end
        if sampleDiag.tg.enstrophy > params.maxStableMeanEnstrophy
            unstable = true; reasons{end+1} = 'enstrophy'; %#ok<AGROW>
        end
        if sampleDiag.tg.highKEnergyFractionVelocity > params.maxStableHighKFraction
            unstable = true; reasons{end+1} = 'highK'; %#ok<AGROW>
        end
        if unstable
            warning('run_projection_taylor_green_demo:Unstable', ...
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
    [state, stepDiag] = mpcd_step_projection_periodic_q9(state, stepParams);
    actualLastStep = nextStep;
    lastStepDiag = stepDiag;

    if willProgress
        lastProgressPrint = actualLastStep;
        fprintf('  TG step %d/%d, t=%.6g, A=%.5g, coh=%.4g, enst=%.4g, lowKdens=%.3e, div=%.3e, mfRes=%.3e, elapsed=%.1fs\n', ...
            actualLastStep, params.nSteps, actualLastStep*params.dt, ...
            sampleDiag.tg.modeAmplitude, sampleDiag.tg.modeCoherence, sampleDiag.tg.enstrophy, ...
            current_lowk_density(NMaps, isamp), finite_or_nan(getf(stepDiag, 'rmsDivParticleAfter')), ...
            finite_or_nan(getf(stepDiag, 'rmsMassFluxDivResidual')), toc(wallClockTic));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_taylor_green_demo:TimeLimit', ...
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
kBTCellSeries = kBTCellSeries(valid);
tgAmplitudeSeries = tgAmplitudeSeries(valid);
tgModeEnergySeries = tgModeEnergySeries(valid);
tgCoherenceSeries = tgCoherenceSeries(valid);
tgResidualEnergySeries = tgResidualEnergySeries(valid);
tgTotalHydroEnergySeries = tgTotalHydroEnergySeries(valid);
enstrophySeries = enstrophySeries(valid);
omegaRmsSeries = omegaRmsSeries(valid);
highKVelocityFractionSeries = highKVelocityFractionSeries(valid);
lowKVelocityFractionSeries = lowKVelocityFractionSeries(valid);
tgModeVelocityFractionSeries = tgModeVelocityFractionSeries(valid);
rmsDivGridSeries = rmsDivGridSeries(valid);
NStdSeries = NStdSeries(valid);
NOutBandSeries = NOutBandSeries(valid);
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
        'excludeWallCells', 0, ...
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
out.kBTCellSeries = kBTCellSeries;
out.tgAmplitudeSeries = tgAmplitudeSeries;
out.tgModeEnergySeries = tgModeEnergySeries;
out.tgCoherenceSeries = tgCoherenceSeries;
out.tgResidualEnergySeries = tgResidualEnergySeries;
out.tgTotalHydroEnergySeries = tgTotalHydroEnergySeries;
out.enstrophySeries = enstrophySeries;
out.omegaRmsSeries = omegaRmsSeries;
out.highKVelocityFractionSeries = highKVelocityFractionSeries;
out.lowKVelocityFractionSeries = lowKVelocityFractionSeries;
out.tgModeVelocityFractionSeries = tgModeVelocityFractionSeries;
out.rmsDivGridSeries = rmsDivGridSeries;
out.NStdSeries = NStdSeries;
out.NOutBandSeries = NOutBandSeries;
out.NMaps = NMaps;
out.omegaMaps = omegaMaps;
out.densityMetrics = metrics;
out.diagHistory = diagHistory;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.elapsedWallClock = elapsedWallClock;
out.diagColumns = {'step','t','rmsDivBefore','rmsDivParticleAfter', ...
    'massFluxDivBefore','massFluxDivTarget','massFluxDivProjectedAfter','massFluxDivResidual', ...
    'massFluxDivParticleAfter','dvRms','popStdAfter','kBTCellAfterProjection', ...
    'tgAmplitude','tgCoherence','enstrophy','highKVelocityFraction','meanUx','meanUy'};
out.summary = summarize_tg_run(out);

fprintf('\n=== run_projection_taylor_green_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('initialPopulationMode      : %s\n', params.initialPopulationMode);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('steps completed            : %d / %d in %.2f s\n', actualLastStep, params.nSteps, elapsedWallClock);
fprintf('TG amplitude/mode          : %.12g / (%d,%d)\n', params.taylorGreenAmplitude, params.taylorGreenModeX, params.taylorGreenModeY);
fprintf('projectionStrength         : %.12g\n', params.projectionStrength);
fprintf('massFluxProjectionMode     : %s  strength=%.12g  beta=%.12g\n', ...
    params.massFluxProjectionMode, params.massFluxProjectionStrength, params.massFluxDensityRelaxationBeta);
fprintf('cleanup strength           : %.12g\n', params.massFluxFinalVelocityProjectionStrength * double(logical(params.massFluxFinalVelocityProjectionCleanup)));
fprintf('mean std(N)                : %.12g\n', out.summary.meanStdN);
fprintf('mean low-k density energy  : %.12g\n', out.summary.meanLowKEnergy);
fprintf('mean div(u) particle after : %.12g\n', out.summary.meanRmsDivParticleAfter);
fprintf('mean mass-flux residual    : %.12g\n', out.summary.meanMassFluxDivResidual);
fprintf('mean/final TG amplitude    : %.12g / %.12g\n', out.summary.meanTGAmplitude, out.summary.finalTGAmplitude);
fprintf('mean/final TG coherence    : %.12g / %.12g\n', out.summary.meanTGCoherence, out.summary.finalTGCoherence);
fprintf('mean/final enstrophy       : %.12g / %.12g\n', out.summary.meanEnstrophy, out.summary.finalEnstrophy);
fprintf('mean high-k vel fraction   : %.12g\n', out.summary.meanHighKVelocityFraction);

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
    getf(getf(stepDiag,'populationAfterProjection',struct()),'stdN'), ...
    getf(stepDiag,'kBTCellAfterProjection'), sampleDiag.tg.modeAmplitude, ...
    sampleDiag.tg.modeCoherence, sampleDiag.tg.enstrophy, ...
    sampleDiag.tg.highKEnergyFractionVelocity, sampleDiag.meanUx, sampleDiag.meanUy];
end

function sampleDiag = tg_sample_diagnostics(state, params, it)
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
thermal = projection_thermal_diagnostics(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
pop = projection_population_diagnostics(state.x, params, ...
    'periodicX', true, 'periodicY', true, 'targetGamma', params.gamma, 'bandFraction', params.densityBandFraction);
tg = projection_taylor_green_diagnostics(G, params);
sampleDiag = struct();
sampleDiag.it = it;
sampleDiag.t = it * params.dt;
sampleDiag.G = G;
sampleDiag.thermal = thermal;
sampleDiag.pop = pop;
sampleDiag.tg = tg;
sampleDiag.meanUx = mean(G.Ux(:), 'omitnan');
sampleDiag.meanUy = mean(G.Uy(:), 'omitnan');
sampleDiag.kBTCell = thermal.kBTCellRelative;
end

function summary = summarize_tg_run(out)
metrics = out.densityMetrics;
if isempty(metrics)
    meanStdN = mean(out.NStdSeries, 'omitnan');
    meanOutBand = mean(out.NOutBandSeries, 'omitnan');
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
summary.meanStdN = meanStdN;
summary.meanOutBand = meanOutBand;
summary.meanLowKEnergy = meanLowK;
summary.timeAvgRelRms = timeAvgRelRms;
summary.meanMeanUx = mean(out.meanUxSeries, 'omitnan');
summary.meanMeanUy = mean(out.meanUySeries, 'omitnan');
summary.finalMeanUx = out.meanUxSeries(end);
summary.finalMeanUy = out.meanUySeries(end);
summary.meanKBTCell = mean(out.kBTCellSeries, 'omitnan');
summary.finalKBTCell = out.kBTCellSeries(end);
summary.initialTGAmplitude = out.tgAmplitudeSeries(1);
summary.meanTGAmplitude = mean(out.tgAmplitudeSeries, 'omitnan');
summary.finalTGAmplitude = out.tgAmplitudeSeries(end);
summary.finalAbsTGAmplitude = abs(out.tgAmplitudeSeries(end));
summary.amplitudeRetention = abs(out.tgAmplitudeSeries(end)) / max(abs(out.tgAmplitudeSeries(1)), eps);
summary.meanTGModeEnergy = mean(out.tgModeEnergySeries, 'omitnan');
summary.finalTGModeEnergy = out.tgModeEnergySeries(end);
summary.modeEnergyRetention = out.tgModeEnergySeries(end) / max(out.tgModeEnergySeries(1), eps);
summary.meanTGCoherence = mean(out.tgCoherenceSeries, 'omitnan');
summary.finalTGCoherence = out.tgCoherenceSeries(end);
summary.meanResidualEnergy = mean(out.tgResidualEnergySeries, 'omitnan');
summary.finalResidualEnergy = out.tgResidualEnergySeries(end);
summary.meanTotalHydroEnergy = mean(out.tgTotalHydroEnergySeries, 'omitnan');
summary.finalTotalHydroEnergy = out.tgTotalHydroEnergySeries(end);
summary.meanEnstrophy = mean(out.enstrophySeries, 'omitnan');
summary.finalEnstrophy = out.enstrophySeries(end);
summary.meanOmegaRms = mean(out.omegaRmsSeries, 'omitnan');
summary.finalOmegaRms = out.omegaRmsSeries(end);
summary.meanHighKVelocityFraction = mean(out.highKVelocityFractionSeries, 'omitnan');
summary.finalHighKVelocityFraction = out.highKVelocityFractionSeries(end);
summary.meanLowKVelocityFraction = mean(out.lowKVelocityFractionSeries, 'omitnan');
summary.finalLowKVelocityFraction = out.lowKVelocityFractionSeries(end);
summary.meanTGModeVelocityFraction = mean(out.tgModeVelocityFractionSeries, 'omitnan');
summary.finalTGModeVelocityFraction = out.tgModeVelocityFractionSeries(end);
summary.meanRmsDivGrid = mean(out.rmsDivGridSeries, 'omitnan');
summary.finalRmsDivGrid = out.rmsDivGridSeries(end);
summary.meanRmsDivParticleAfter = mean_diag_col(diag, cols, 'rmsDivParticleAfter');
summary.meanMassFluxDivResidual = mean_diag_col(diag, cols, 'massFluxDivResidual');
summary.meanMassFluxDivParticleAfter = mean_diag_col(diag, cols, 'massFluxDivParticleAfter');
summary.meanDensityTransportProjectedRms = mean_diag_col(diag, cols, 'densityTransportProjectedRms');
summary.meanDvRms = mean_diag_col(diag, cols, 'dvRms');
end

function v = mean_diag_col(diag, cols, name)
id = find(strcmp(cols, name), 1);
if isempty(id) || isempty(diag)
    v = NaN;
else
    v = mean(diag(:, id), 'omitnan');
end
end

function make_figures(out)
t = out.sampleTimes;
figure('Name','Taylor-Green diagnostics','Color','w');
subplot(2,2,1);
plot(t, out.tgAmplitudeSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('TG amplitude'); grid on;
subplot(2,2,2);
plot(t, out.tgCoherenceSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('mode coherence'); grid on;
subplot(2,2,3);
plot(t, out.enstrophySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('enstrophy'); grid on;
subplot(2,2,4);
plot(t, out.highKVelocityFractionSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('high-k velocity fraction'); grid on;

if ~isempty(out.omegaMaps)
    figure('Name','Final Taylor-Green vorticity','Color','w');
    imagesc(out.omegaMaps(:,:,end).'); axis image xy; colorbar;
    title('Final omega'); xlabel('ix'); ylabel('iy');
end
end

function y = current_lowk_density(NMaps, isamp)
if isempty(NMaps) || isamp < 1 || size(NMaps,3) < isamp
    y = NaN;
else
    N = NMaps(:,:,isamp);
    rel = N / max(mean(N(:),'omitnan'), eps) - 1;
    H = fft2(rel);
    y = sum(abs(H(:)).^2) / numel(H)^2;
end
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 1.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 32);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'nSteps', 10000);
params = set_default(params, 'sampleEvery', 50);
params = set_default(params, 'progressEvery', 1000);
params = set_default(params, 'maxWallClockSeconds', Inf);
params = set_default(params, 'dt', 2.0e-3);
params = set_default(params, 'kBT', 0.05);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.0);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'useRandomGridShift', true);
params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'taylorGreenAmplitude', 0.20);
params = set_default(params, 'taylorGreenModeX', 1);
params = set_default(params, 'taylorGreenModeY', 1);
params = set_default(params, 'taylorGreenThermalNoise', true);
params = set_default(params, 'taylorGreenHighKCut', max(4, min(params.Nx, params.Ny)/4));
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
params = set_default(params, 'storeVorticityMaps', false);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'seed', 11);
params = set_default(params, 'abortOnUnstable', true);
params = set_default(params, 'maxStableAbsMeanUx', 1.0);
params = set_default(params, 'maxStableAbsMeanUy', 1.0);
params = set_default(params, 'maxStableMeanEnstrophy', 1.0e4);
params = set_default(params, 'maxStableHighKFraction', 0.80);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function v = getf(s, name, defaultValue)
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
