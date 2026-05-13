function out = run_projection_piston_demo(params)
%RUN_PROJECTION_PISTON_DEMO Moving-piston validation for Q6/Q9 projection.
%
%   out = run_projection_piston_demo(params)
%
% Minimal Q9-native piston case.  The top wall moves from pistonY0 toward
% pistonYMin.  The active computational height is yTop(t), with fixed Nx-by-Ny
% resolution.  This first piston benchmark tests whether Q9 controls coherent
% density modes during a slow compression without reactivating the historical
% redistribution/liquid-closure operator.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);

rng(params.seed);
[state, initInfo] = projection_initialize_particles(params);
state.piston = struct('yTop', params.pistonY0, 'yPrev', params.pistonY0, ...
    'Up', params.pistonVy, 'Ly0', params.Ly, 'compression', 0.0, ...
    'activeHeight', params.pistonY0);
Np = initInfo.Np;

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
sampleTimes = nan(nSamples, 1);
sampleSteps = nan(nSamples, 1);
yTopSeries = nan(nSamples, 1);
compressionSeries = nan(nSamples, 1);
rhoPhysicalSeries = nan(nSamples, 1);
PkinMeanSeries = nan(nSamples, 1);
kBTCellSeries = nan(nSamples, 1);
UxProfiles = nan(params.Ny, nSamples);
UyProfiles = nan(params.Ny, nSamples);
NProfiles = nan(params.Ny, nSamples);
if params.storeDensityMaps
    NMaps = nan(params.Nx, params.Ny, nSamples);
else
    NMaps = [];
end

% Columns are listed in out.diagColumns.
diagHistory = nan(nSamples, 21);

isamp = 0;
actualLastStep = 0;
stopRequested = false;
wallClockTic = tic;
lastProgressPrint = -inf;
lastStepDiag = [];

for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        sampleParams = active_params_from_state(params, state);
        sampleDiag = piston_sample_diagnostics(state, sampleParams, it);

        sampleTimes(isamp) = it * params.dt;
        sampleSteps(isamp) = it;
        yTopSeries(isamp) = sampleDiag.yTop;
        compressionSeries(isamp) = sampleDiag.compression;
        rhoPhysicalSeries(isamp) = sampleDiag.rhoPhysicalMean;
        PkinMeanSeries(isamp) = sampleDiag.PkinMean;
        kBTCellSeries(isamp) = sampleDiag.kBTCell;
        UxProfiles(:, isamp) = mean(sampleDiag.G.Ux, 1).';
        UyProfiles(:, isamp) = mean(sampleDiag.G.Uy, 1).';
        NProfiles(:, isamp) = mean(sampleDiag.G.N, 1).';
        if params.storeDensityMaps
            NMaps(:, :, isamp) = double(sampleDiag.G.N);
        end

        if it == 0
            diagHistory(isamp, :) = [it, it*params.dt, sampleDiag.yTop, sampleDiag.compression, ...
                NaN, NaN, NaN, NaN, NaN, NaN, ...
                sampleDiag.stdN, sampleDiag.stdN, sampleDiag.kBTCell, NaN, NaN, ...
                NaN, NaN, sampleDiag.meanVy, sampleDiag.rhoPhysicalMean, sampleDiag.PkinMean, NaN];
        elseif ~isempty(lastStepDiag)
            diagHistory(isamp, :) = step_diag_row(it, params, lastStepDiag, sampleDiag);
        end
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

    [state, stepDiag] = mpcd_step_projection_piston(state, stepParams, nextStep);
    actualLastStep = nextStep;
    lastStepDiag = stepDiag;

    if willProgress
        lastProgressPrint = actualLastStep;
        fprintf('  piston step %d/%d, t=%.6g, yTop=%.6g, comp=%.4g, div red=%.3e, massFlux res=%.3e, elapsed=%.1fs\n', ...
            actualLastStep, params.nSteps, actualLastStep*params.dt, ...
            stepDiag.piston.yTop, stepDiag.piston.compression, ...
            finite_or_nan(stepDiag.divReductionProjected), ...
            finite_or_nan(stepDiag.rmsMassFluxDivResidual), toc(wallClockTic));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_piston_demo:TimeLimit', ...
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
yTopSeries = yTopSeries(valid);
compressionSeries = compressionSeries(valid);
rhoPhysicalSeries = rhoPhysicalSeries(valid);
PkinMeanSeries = PkinMeanSeries(valid);
kBTCellSeries = kBTCellSeries(valid);
UxProfiles = UxProfiles(:, valid);
UyProfiles = UyProfiles(:, valid);
NProfiles = NProfiles(:, valid);
if params.storeDensityMaps
    NMaps = NMaps(:, :, valid);
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
out.yTopSeries = yTopSeries;
out.compressionSeries = compressionSeries;
out.rhoPhysicalSeries = rhoPhysicalSeries;
out.PkinMeanSeries = PkinMeanSeries;
out.kBTCellSeries = kBTCellSeries;
out.UxProfiles = UxProfiles;
out.UyProfiles = UyProfiles;
out.NProfiles = NProfiles;
out.NMaps = NMaps;
out.densityMetrics = metrics;
out.diagHistory = diagHistory;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.elapsedWallClock = elapsedWallClock;
out.diagColumns = {'step','t','yTop','compression','rmsDivBefore','rmsDivParticleAfter', ...
    'massFluxDivBefore','massFluxDivTarget','massFluxDivProjectedAfter','massFluxDivResidual', ...
    'popStdClassic','popStdProjection','kBTCellAfter','densityTransportProjectedRms', ...
    'dvRms','nTopHits','dPyTop','meanVyAfter','rhoPhysicalMean','PkinMean', ...
    'massFluxDivParticleAfter'};
out.summary = summarize_piston_run(out);

fprintf('\n=== run_projection_piston_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('initialPopulationMode      : %s\n', initInfo.mode);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('steps completed            : %d / %d in %.2f s\n', out.actualLastStep, params.nSteps, out.elapsedWallClock);
fprintf('piston yTop final          : %.12g\n', out.summary.finalYTop);
fprintf('compression final          : %.12g\n', out.summary.finalCompression);
fprintf('projectionStrength         : %.6g\n', params.projectionStrength);
fprintf('massFluxProjectionMode     : %s  strength=%.6g  beta=%.6g\n', ...
    char(params.massFluxProjectionMode), params.massFluxProjectionStrength, params.massFluxDensityRelaxationBeta);
fprintf('mean std(N)                : %.12g\n', out.summary.meanStdN);
fprintf('mean low-k energy          : %.12g\n', out.summary.meanLowKEnergy);
fprintf('mean rho transport proj    : %.12g\n', out.summary.meanDensityTransportProjectedRms);
fprintf('mean div(u) particle after : %.12g\n', out.summary.meanRmsDivParticleAfter);
fprintf('mean mass-flux residual    : %.12g\n', out.summary.meanMassFluxDivResidual);
fprintf('mean mass-flux particle aft: %.12g\n', out.summary.meanMassFluxDivParticleAfter);
fprintf('mean Pkin                  : %.12g\n', out.summary.meanPkinMean);

if params.makeFigures
    make_figures(out);
end
end

function row = step_diag_row(it, params, stepDiag, sampleDiag)
row = [it, it*params.dt, sampleDiag.yTop, sampleDiag.compression, ...
    getf(stepDiag,'rmsDivBefore'), getf(stepDiag,'rmsDivParticleAfter'), ...
    getf(stepDiag,'rmsMassFluxDivBefore'), getf(stepDiag,'rmsMassFluxDivTarget'), ...
    getf(stepDiag,'rmsMassFluxDivProjectedAfter'), getf(stepDiag,'rmsMassFluxDivResidual'), ...
    getf(getf(stepDiag,'populationAfterClassic',struct()),'stdN'), ...
    getf(getf(stepDiag,'populationAfterProjection',struct()),'stdN'), ...
    getf(stepDiag,'kBTCellAfterProjection'), getf(stepDiag,'densityTransportProjectedRms'), ...
    getf(stepDiag,'dvRms'), getf(getf(stepDiag,'wallInfo',struct()),'nTop'), ...
    getf(getf(stepDiag,'wallInfo',struct()),'dPyTop'), getf(stepDiag,'meanVyAfterProjection'), ...
    sampleDiag.rhoPhysicalMean, sampleDiag.PkinMean, getf(stepDiag,'rmsMassFluxDivParticleAfter')];
end

function sampleDiag = piston_sample_diagnostics(state, params, it)
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
thermal = projection_thermal_diagnostics(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
if isfield(state, 'piston') && isfield(state.piston, 'yTop')
    yTop = state.piston.yTop;
    Ly0 = state.piston.Ly0;
else
    yTop = params.Ly;
    Ly0 = get_param(params, 'Ly0', params.Ly);
end
rhoPhysicalMean = size(state.x, 1) / max(params.Lx * yTop, eps);
Tmap = local_temperature_map(state.x, state.v, params);
PkinMean = mean(G.rho(:) .* max(Tmap(:), 0), 'omitnan');
if ~isfinite(PkinMean)
    PkinMean = rhoPhysicalMean * thermal.kBTCellRelative;
end
sampleDiag = struct();
sampleDiag.it = it;
sampleDiag.yTop = yTop;
sampleDiag.compression = 1.0 - yTop / Ly0;
sampleDiag.rhoPhysicalMean = rhoPhysicalMean;
sampleDiag.PkinMean = PkinMean;
sampleDiag.kBTCell = thermal.kBTCellRelative;
sampleDiag.meanVy = mean(state.v(:,2));
sampleDiag.stdN = std(double(G.N(:)));
sampleDiag.G = G;
sampleDiag.thermal = thermal;
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

function params = active_params_from_state(params0, state)
params = params0;
params.Ly0 = get_param(params0, 'Ly0', get_param(params0, 'LyReference', params0.Ly));
if isfield(state, 'piston') && isfield(state.piston, 'yTop') && ~isempty(state.piston.yTop)
    params.Ly = state.piston.yTop;
else
    params.Ly = params.pistonY0;
end
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Ly0', params.Ly);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 16);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 1.0);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.0);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'nSteps', 3000);
params = set_default(params, 'sampleEvery', 50);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'massFluxProjectionMode', 'off');
params = set_default(params, 'massFluxProjectionStrength', params.projectionStrength);
params = set_default(params, 'massFluxDensityRelaxationBeta', 0.0);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', false);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', false);
params = set_default(params, 'massFluxFinalVelocityProjectionStrength', 1.0);
params = set_default(params, 'massFluxTargetFilter', 'none');
params = set_default(params, 'massFluxLowKMaxIndex', 2);
params = set_default(params, 'massFluxProjectionRegularization', 1e-12);
params = set_default(params, 'massFluxMinCellCount', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'projectionTransportDiagnosticsEnable', true);
params = set_default(params, 'densityTransportDiagnosticsEnable', true);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'densityBandFraction', 0.20);
params = set_default(params, 'densityExcludeWallCells', 0);
params = set_default(params, 'lowKMaxIndex', 2);
params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'thermostatMinParticlesPerCell', 2);
params = set_default(params, 'thermostatMaxScale', 5.0);
params = set_default(params, 'wallModeY', 'thermalize');
params = set_default(params, 'wallSigma', sqrt(params.kBT));
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'pistonTangentialMode', 'specular');
params = set_default(params, 'pistonY0', params.Ly);
params = set_default(params, 'pistonCompressionTarget', 0.10);
params = set_default(params, 'pistonYMin', params.Ly * (1.0 - params.pistonCompressionTarget));
if ~isfield(params, 'pistonVy') || isempty(params.pistonVy)
    params.pistonVy = -(params.pistonY0 - params.pistonYMin) / max(params.nSteps * params.dt, eps);
end
params = set_default(params, 'pistonStopOnMin', true);
params = set_default(params, 'pistonAffineReposition', true);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'seed', 11);
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'progressEvery', 500);
params = set_default(params, 'maxWallClockSeconds', Inf);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function summary = summarize_piston_run(out)
H = out.diagHistory;
summary = struct();
summary.finalYTop = last_or_nan(out.yTopSeries);
summary.finalCompression = last_or_nan(out.compressionSeries);
summary.finalRhoPhysicalMean = last_or_nan(out.rhoPhysicalSeries);
summary.meanRhoPhysicalMean = mean(out.rhoPhysicalSeries, 'omitnan');
summary.meanPkinMean = mean(out.PkinMeanSeries, 'omitnan');
summary.finalPkinMean = last_or_nan(out.PkinMeanSeries);
summary.meanKBTCell = mean(out.kBTCellSeries, 'omitnan');
if ~isempty(H)
    summary.meanRmsDivBefore = mean(H(:,5), 'omitnan');
    summary.meanRmsDivParticleAfter = mean(H(:,6), 'omitnan');
    summary.meanDivReductionParticle = summary.meanRmsDivParticleAfter / max(summary.meanRmsDivBefore, eps);
    summary.meanDensityTransportProjectedRms = mean(H(:,14), 'omitnan');
    summary.meanMassFluxDivBefore = mean(H(:,7), 'omitnan');
    summary.meanMassFluxDivTarget = mean(H(:,8), 'omitnan');
    summary.meanMassFluxDivProjectedAfter = mean(H(:,9), 'omitnan');
    summary.meanMassFluxDivResidual = mean(H(:,10), 'omitnan');
    summary.meanMassFluxDivParticleAfter = mean(H(:,21), 'omitnan');
    summary.meanPopStdClassic = mean(H(:,11), 'omitnan');
    summary.meanPopStdProjection = mean(H(:,12), 'omitnan');
    summary.meanTopHits = mean(H(:,16), 'omitnan');
else
    summary.meanRmsDivBefore = NaN;
    summary.meanRmsDivParticleAfter = NaN;
    summary.meanDivReductionParticle = NaN;
    summary.meanDensityTransportProjectedRms = NaN;
    summary.meanMassFluxDivBefore = NaN;
    summary.meanMassFluxDivTarget = NaN;
    summary.meanMassFluxDivProjectedAfter = NaN;
    summary.meanMassFluxDivResidual = NaN;
    summary.meanMassFluxDivParticleAfter = NaN;
    summary.meanPopStdClassic = NaN;
    summary.meanPopStdProjection = NaN;
    summary.meanTopHits = NaN;
end
if isfield(out, 'densityMetrics') && ~isempty(out.densityMetrics)
    summary.meanStdN = out.densityMetrics.summary.meanStdN;
    summary.finalStdN = out.densityMetrics.summary.finalStdN;
    summary.meanOutBandFraction = out.densityMetrics.summary.meanOutBandFraction;
    summary.meanLowKEnergy = out.densityMetrics.summary.meanLowKEnergy;
    summary.timeAvgRelRms = out.densityMetrics.summary.timeAvgRelRms;
else
    summary.meanStdN = NaN;
    summary.finalStdN = NaN;
    summary.meanOutBandFraction = NaN;
    summary.meanLowKEnergy = NaN;
    summary.timeAvgRelRms = NaN;
end
end

function make_figures(out)
figure('Name', 'Piston kinematics');
tiledlayout(3,1);
nexttile; plot(out.sampleTimes, out.yTopSeries, '-o'); grid on; ylabel('yTop');
nexttile; plot(out.sampleTimes, out.compressionSeries, '-o'); grid on; ylabel('compression');
nexttile; plot(out.sampleTimes, out.rhoPhysicalSeries, '-o'); grid on; ylabel('rho mean'); xlabel('t');

if ~isempty(out.diagHistory)
    H = out.diagHistory;
    figure('Name', 'Piston Q6/Q9 diagnostics');
    tiledlayout(4,1);
    nexttile; plot(H(:,2), H(:,11), '-', H(:,2), H(:,12), '--'); grid on; ylabel('std N'); legend('classic','projection');
    nexttile; semilogy(H(:,2), H(:,5), '-', H(:,2), H(:,6), '--'); grid on; ylabel('div u'); legend('before','after particle');
    nexttile; semilogy(H(:,2), H(:,7), '-', H(:,2), H(:,10), '--'); grid on; ylabel('div Nu'); legend('before','residual');
    nexttile; plot(H(:,2), H(:,20), '-'); grid on; ylabel('Pkin mean'); xlabel('t');
end
end

function y = finite_or_nan(x)
if isempty(x) || ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
    y = NaN;
else
    y = x;
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
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

function v = last_or_nan(x)
if isempty(x)
    v = NaN;
else
    idx = find(isfinite(x), 1, 'last');
    if isempty(idx)
        v = NaN;
    else
        v = x(idx);
    end
end
end
