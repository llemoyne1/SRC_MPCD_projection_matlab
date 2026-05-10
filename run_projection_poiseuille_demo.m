function out = run_projection_poiseuille_demo(params)
%RUN_PROJECTION_POISEUILLE_DEMO Poiseuille prototype with pressure projection.
%
%   out = run_projection_poiseuille_demo(params)
%
% Standalone Q3 demo: x-periodic channel, y walls, body force, classic SRD,
% algebraic pressure projection, and population diagnostics before/after the
% projection substep.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);

rng(params.seed);
[state, initInfo] = projection_initialize_particles(params);
Np = initInfo.Np;

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
sampleTimes = nan(nSamples, 1);
sampleSteps = nan(nSamples, 1);
UxProfiles = nan(params.Ny, nSamples);
UyProfiles = nan(params.Ny, nSamples);
NProfiles = nan(params.Ny, nSamples);
if params.storeDensityMaps
    NMaps = nan(params.Nx, params.Ny, nSamples);
else
    NMaps = [];
end
diagHistory = nan(nSamples, 51);
% Main columns are documented in out.diagColumns below. The first 24 columns
% are the original Q3/Q3b diagnostics; later columns add thermal, thermostat
% and continuous density-transport diagnostics.

isamp = 0;
stopRequested = false;
actualLastStep = 0;
wallClockTic = tic;
lastProgressPrint = -inf;
for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
            'periodicX', true, 'periodicY', false, 'minCount', 1);
        UxProfiles(:, isamp) = mean(G.Ux, 1).';
        UyProfiles(:, isamp) = mean(G.Uy, 1).';
        NProfiles(:, isamp) = mean(G.N, 1).';
        if params.storeDensityMaps
            NMaps(:, :, isamp) = double(G.N);
        end
        sampleTimes(isamp) = it * params.dt;
        sampleSteps(isamp) = it;
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

    [state, stepDiag] = mpcd_step_projection_poiseuille(state, stepParams);
    actualLastStep = nextStep;

    if willProgress
        lastProgressPrint = actualLastStep;
        if isfield(stepDiag, 'kBTCellAfterProjection') && isfinite(stepDiag.kBTCellAfterProjection)
            kBTForPrint = stepDiag.kBTCellAfterProjection;
        else
            kBTForPrint = NaN;
        end
        fprintf('  step %d/%d, t=%.6g, elapsed=%.1fs, div red=%.3e, kBTcell=%.3g, fullDiag=%d\n', ...
            actualLastStep, params.nSteps, actualLastStep*params.dt, toc(wallClockTic), ...
            stepDiag.divReductionProjected, kBTForPrint, stepParams.computeDiagnostics);
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_poiseuille_demo:TimeLimit', ...
            'Stopping early at step %d/%d after %.1f s because maxWallClockSeconds=%.1f.', ...
            actualLastStep, params.nSteps, toc(wallClockTic), params.maxWallClockSeconds);
        stopRequested = true;
    end

    if willSample
        % Fill the row that has just been created on next loop pass. The
        % diagnostic arrays are aligned to sample step after this step.
        row = floor(nextStep / params.sampleEvery) + 1;
        if row <= nSamples
            diagHistory(row, :) = [nextStep, nextStep*params.dt, ...
                stepDiag.rmsDivBefore, stepDiag.rmsDivProjectedAfter, stepDiag.rmsDivParticleAfter, ...
                stepDiag.kBTBeforeProjection, stepDiag.kBTAfterProjection, ...
                stepDiag.populationBeforeStep.stdN, stepDiag.populationAfterClassic.stdN, stepDiag.populationAfterProjection.stdN, ...
                stepDiag.populationBeforeStep.nEmptyCells, stepDiag.populationAfterClassic.nEmptyCells, stepDiag.populationAfterProjection.nEmptyCells, ...
                stepDiag.populationProjectionDeltaMaxAbs, stepDiag.dvRms, ...
                stepDiag.populationTransportClassicDeltaRms, stepDiag.populationTransportProjectedDeltaRms, ...
                stepDiag.populationTransportProjectedMinusClassicRms, ...
                stepDiag.populationTransportClassicDeltaMaxAbs, stepDiag.populationTransportProjectedDeltaMaxAbs, ...
                stepDiag.populationTransportProjectedMinusClassicMaxAbs, ...
                stepDiag.populationTransportClassicStdAfter, stepDiag.populationTransportProjectedStdAfter, ...
                stepDiag.populationTransportProjectedVsClassicStdDelta, ...
                stepDiag.kBTGlobalBeforeProjection, stepDiag.kBTGlobalAfterProjectionRaw, stepDiag.kBTGlobalAfterProjection, ...
                stepDiag.kBTCellBeforeProjection, stepDiag.kBTCellAfterProjectionRaw, stepDiag.kBTCellAfterProjection, ...
                stepDiag.hydroKEBeforeProjection, stepDiag.hydroKEAfterProjectionRaw, stepDiag.hydroKEAfterProjection, ...
                stepDiag.thermalKEBeforeProjection, stepDiag.thermalKEAfterProjectionRaw, stepDiag.thermalKEAfterProjection, ...
                stepDiag.totalKEBeforeProjection, stepDiag.totalKEAfterProjectionRaw, stepDiag.totalKEAfterProjection, ...
                stepDiag.thermostatRmsVelocityChange, stepDiag.thermostatMeanScale, stepDiag.thermostatNCells, ...
                stepDiag.densityTransportClassicRms, stepDiag.densityTransportProjectedRms, ...
                stepDiag.densityTransportProjectedMinusClassicRms, ...
                stepDiag.densityTransportClassicMaxAbs, stepDiag.densityTransportProjectedMaxAbs, ...
                stepDiag.densityTransportProjectedMinusClassicMaxAbs, ...
                stepDiag.densityTransportClassicPredictedStdAfter, stepDiag.densityTransportProjectedPredictedStdAfter, ...
                stepDiag.densityTransportProjectedVsClassicStdDelta];
        end
    end

    if stopRequested
        break;
    end
end

if actualLastStep == 0
    actualLastStep = min(params.nSteps, max(sampleSteps(isfinite(sampleSteps))));
end
elapsedWallClock = toc(wallClockTic);

% Trim, in case sample schedule changed.
valid = isfinite(sampleTimes);
sampleTimes = sampleTimes(valid);
sampleSteps = sampleSteps(valid);
UxProfiles = UxProfiles(:, valid);
UyProfiles = UyProfiles(:, valid);
NProfiles = NProfiles(:, valid);
if params.storeDensityMaps
    NMaps = NMaps(:, :, valid);
end
diagHistory = diagHistory(isfinite(diagHistory(:, 1)), :);

yCenters = ((0:params.Ny-1).' + 0.5) * params.Ly / params.Ny;

out = struct();
out.params = params;
out.initialization = initInfo;
out.state = state;
out.sampleTimes = sampleTimes;
out.sampleSteps = sampleSteps;
out.yCenters = yCenters;
out.UxProfiles = UxProfiles;
out.UyProfiles = UyProfiles;
out.NProfiles = NProfiles;
out.NMaps = NMaps;
out.diagHistory = diagHistory;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.elapsedWallClock = elapsedWallClock;
out.diagColumns = {'step','t','rmsDivBefore','rmsDivProjected','rmsDivParticle', ...
    'kBTBefore','kBTAfter','popStdBefore','popStdClassic','popStdProjection', ...
    'emptyBefore','emptyClassic','emptyProjection','popDeltaProjectionMaxAbs','dvRms', ...
    'popForecastClassicDeltaRms','popForecastProjectedDeltaRms', ...
    'popForecastProjectedMinusClassicRms', ...
    'popForecastClassicDeltaMax','popForecastProjectedDeltaMax', ...
    'popForecastProjectedMinusClassicMax', ...
    'popForecastClassicStdAfter','popForecastProjectedStdAfter', ...
    'popForecastProjectedVsClassicStdDelta', ...
    'kBTGlobalBefore','kBTGlobalAfterRaw','kBTGlobalAfter', ...
    'kBTCellBefore','kBTCellAfterRaw','kBTCellAfter', ...
    'hydroKEBefore','hydroKEAfterRaw','hydroKEAfter', ...
    'thermalKEBefore','thermalKEAfterRaw','thermalKEAfter', ...
    'totalKEBefore','totalKEAfterRaw','totalKEAfter', ...
    'thermostatRmsDV','thermostatMeanScale','thermostatNCells', ...
    'densityTransportClassicRms','densityTransportProjectedRms', ...
    'densityTransportProjectedMinusClassicRms', ...
    'densityTransportClassicMax','densityTransportProjectedMax', ...
    'densityTransportProjectedMinusClassicMax', ...
    'densityTransportClassicStdAfter','densityTransportProjectedStdAfter', ...
    'densityTransportProjectedVsClassicStdDelta'};
out.summary = summarize_run(out);
out.viscosity = analyze_projection_poiseuille_viscosity(out, ...
    'excludeWallCells', params.excludeWallCellsFit, ...
    'fitStartFraction', params.fitStartFraction);

fprintf('\n=== run_projection_poiseuille_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('initialPopulationMode      : %s\n', initInfo.mode);
fprintf('initial std(N), outband20  : %.6g / %.6g\n', initInfo.initialStdN, initInfo.initialOutBand20);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('steps completed            : %d / %d in %.2f s\n', out.actualLastStep, params.nSteps, out.elapsedWallClock);
if out.stoppedEarly
    fprintf('stopped early              : 1\n');
end
fprintf('wallModeY                  : %s\n', params.wallModeY);
fprintf('projectionStrength         : %.6g\n', params.projectionStrength);
fprintf('last rms div before        : %.12e\n', out.summary.lastRmsDivBefore);
fprintf('last rms div projected     : %.12e\n', out.summary.lastRmsDivProjected);
fprintf('last rms div particle grid : %.12e\n', out.summary.lastRmsDivParticle);
fprintf('last kBT after projection  : %.12e\n', out.summary.lastKBTAfter);
fprintf('last kBT cell raw/final    : %.12e / %.12e\n', out.summary.lastKBTCellAfterRaw, out.summary.lastKBTCellAfter);
fprintf('thermostat enabled         : %d\n', params.thermostatAfterProjection);
fprintf('thermostat mean scale      : %.6g\n', out.summary.lastThermostatMeanScale);
fprintf('last pop std classic/proj  : %.6g / %.6g\n', out.summary.lastPopStdClassic, out.summary.lastPopStdProjection);
fprintf('max pop delta projection   : %.6g\n', out.summary.maxPopDeltaProjection);
fprintf('forecast pop rms classic/proj/diff : %.6g / %.6g / %.6g\n', ...
    out.summary.lastPopForecastClassicDeltaRms, ...
    out.summary.lastPopForecastProjectedDeltaRms, ...
    out.summary.lastPopForecastProjectedMinusClassicRms);
fprintf('forecast pop std classic/proj      : %.6g / %.6g\n', ...
    out.summary.lastPopForecastClassicStdAfter, ...
    out.summary.lastPopForecastProjectedStdAfter);
fprintf('continuous rho rms classic/proj/diff: %.6g / %.6g / %.6g\n', ...
    out.summary.lastDensityTransportClassicRms, ...
    out.summary.lastDensityTransportProjectedRms, ...
    out.summary.lastDensityTransportProjectedMinusClassicRms);
fprintf('nu_eff fit                 : %.12e\n', out.viscosity.nuEff);
fprintf('fit R2                     : %.6f\n', out.viscosity.R2);

if params.makeFigures
    make_figures(out);
end
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 16);
params = set_default(params, 'gamma', 10);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 1.0);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.02);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'nSteps', 500);
params = set_default(params, 'sampleEvery', 10);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'projectionTransportDiagnosticsEnable', true);
params = set_default(params, 'densityTransportDiagnosticsEnable', true);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'storeDensityMaps', false);
params = set_default(params, 'initialPopulationMode', 'random_uniform');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'thermostatAfterProjection', false);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'thermostatMinParticlesPerCell', 2);
params = set_default(params, 'thermostatMaxScale', 10.0);
params = set_default(params, 'fitStartFraction', 0.5);
params = set_default(params, 'wallModeY', 'bounceback');
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'seed', 2);
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'progressEvery', 0);
params = set_default(params, 'maxWallClockSeconds', Inf);
params = set_default(params, 'excludeWallCellsFit', 2);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function summary = summarize_run(out)
H = out.diagHistory;
summary = struct();
if isempty(H)
    summary.lastRmsDivBefore = NaN;
    summary.lastRmsDivProjected = NaN;
    summary.lastRmsDivParticle = NaN;
    summary.lastKBTAfter = NaN;
    summary.lastPopStdClassic = NaN;
    summary.lastPopStdProjection = NaN;
    summary.maxPopDeltaProjection = NaN;
    summary.lastKBTCellAfterRaw = NaN;
    summary.lastKBTCellAfter = NaN;
    summary.lastThermostatMeanScale = NaN;
    summary.lastDensityTransportClassicRms = NaN;
    summary.lastDensityTransportProjectedRms = NaN;
    summary.lastDensityTransportProjectedMinusClassicRms = NaN;
    return;
end
summary.lastRmsDivBefore = H(end, 3);
summary.lastRmsDivProjected = H(end, 4);
summary.lastRmsDivParticle = H(end, 5);
summary.lastKBTAfter = H(end, 7);
summary.lastPopStdClassic = H(end, 9);
summary.lastPopStdProjection = H(end, 10);
summary.maxPopDeltaProjection = max(H(:, 14));
summary.lastPopForecastClassicDeltaRms = H(end, 16);
summary.lastPopForecastProjectedDeltaRms = H(end, 17);
summary.lastPopForecastProjectedMinusClassicRms = H(end, 18);
summary.lastPopForecastProjectedMinusClassicMax = H(end, 21);
summary.lastPopForecastClassicStdAfter = H(end, 22);
summary.lastPopForecastProjectedStdAfter = H(end, 23);
summary.lastPopForecastProjectedVsClassicStdDelta = H(end, 24);
summary.meanPopForecastClassicDeltaRms = mean(H(:, 16), 'omitnan');
summary.meanPopForecastProjectedDeltaRms = mean(H(:, 17), 'omitnan');
summary.meanPopForecastProjectedMinusClassicRms = mean(H(:, 18), 'omitnan');
summary.maxPopForecastProjectedMinusClassicMax = max(H(:, 21));
summary.lastKBTGlobalBefore = H(end, 25);
summary.lastKBTGlobalAfterRaw = H(end, 26);
summary.lastKBTGlobalAfter = H(end, 27);
summary.lastKBTCellBefore = H(end, 28);
summary.lastKBTCellAfterRaw = H(end, 29);
summary.lastKBTCellAfter = H(end, 30);
summary.lastHydroKEBefore = H(end, 31);
summary.lastHydroKEAfterRaw = H(end, 32);
summary.lastHydroKEAfter = H(end, 33);
summary.lastThermalKEBefore = H(end, 34);
summary.lastThermalKEAfterRaw = H(end, 35);
summary.lastThermalKEAfter = H(end, 36);
summary.lastTotalKEBefore = H(end, 37);
summary.lastTotalKEAfterRaw = H(end, 38);
summary.lastTotalKEAfter = H(end, 39);
summary.lastThermostatRmsDV = H(end, 40);
summary.lastThermostatMeanScale = H(end, 41);
summary.lastThermostatNCells = H(end, 42);
summary.lastDensityTransportClassicRms = H(end, 43);
summary.lastDensityTransportProjectedRms = H(end, 44);
summary.lastDensityTransportProjectedMinusClassicRms = H(end, 45);
summary.lastDensityTransportClassicMax = H(end, 46);
summary.lastDensityTransportProjectedMax = H(end, 47);
summary.lastDensityTransportProjectedMinusClassicMax = H(end, 48);
summary.lastDensityTransportClassicStdAfter = H(end, 49);
summary.lastDensityTransportProjectedStdAfter = H(end, 50);
summary.lastDensityTransportProjectedVsClassicStdDelta = H(end, 51);
summary.meanKBTCellAfter = mean(H(:, 30), 'omitnan');
summary.meanDensityTransportClassicRms = mean(H(:, 43), 'omitnan');
summary.meanDensityTransportProjectedRms = mean(H(:, 44), 'omitnan');
summary.meanDensityTransportProjectedMinusClassicRms = mean(H(:, 45), 'omitnan');
summary.meanPopStdClassic = mean(H(:, 9), 'omitnan');
summary.meanPopStdProjection = mean(H(:, 10), 'omitnan');
summary.meanDivReductionParticle = mean(H(:, 5) ./ max(H(:, 3), eps), 'omitnan');
end

function make_figures(out)
figure('Name', 'Projection Poiseuille profile');
hold on;
plot(out.yCenters, out.viscosity.UxMean, 'o-', 'DisplayName', 'time-avg Ux');
plot(out.yCenters, out.viscosity.UxFit, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('fit nu=%.4g', out.viscosity.nuEff));
xlabel('y'); ylabel('Ux'); grid on; legend('Location', 'best');
title('Poiseuille profile with pressure projection');

if ~isempty(out.diagHistory)
    H = out.diagHistory;
    figure('Name', 'Projection Poiseuille diagnostics');
    tiledlayout(5, 1);
    nexttile;
    semilogy(H(:, 2), H(:, 3), '-', H(:, 2), H(:, 5), '-');
    ylabel('rms div'); grid on; legend('before', 'particle after');
    nexttile;
    plot(H(:, 2), H(:, 7), '-', H(:, 2), H(:, 30), '--'); ylabel('kBT'); grid on; legend('global after', 'cell-relative after');
    nexttile;
    plot(H(:, 2), H(:, 9), '-', H(:, 2), H(:, 10), '--');
    ylabel('population std'); grid on; legend('after classic', 'after projection');
    nexttile;
    plot(H(:, 2), H(:, 16), '-', H(:, 2), H(:, 17), '--', H(:, 2), H(:, 18), ':');
    ylabel('forecast pop rms'); grid on; legend('classic-now', 'projected-now', 'projected-classic');
    nexttile;
    plot(H(:, 2), H(:, 43), '-', H(:, 2), H(:, 44), '--', H(:, 2), H(:, 45), ':');
    xlabel('t'); ylabel('cont. rho rms'); grid on; legend('classic', 'projected', 'diff');
end
end
