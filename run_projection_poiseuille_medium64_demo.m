function out = run_projection_poiseuille_medium64_demo(params)
%RUN_PROJECTION_POISEUILLE_MEDIUM64_DEMO Poiseuille validation at TG medium64 parameters.
%
%   out = run_projection_poiseuille_medium64_demo(params)
%
% Minimal channel benchmark aligned with the validated forced Taylor--Green
% medium64 setup: Nx=Ny=64, gamma=20, dt=1e-3, kBT=0.01, exact per-cell
% initialization, corrected cell thermostat, and optional Q6/Q9 projections.
%
% Supported methods:
%   params.method = 'classic', 'q6', or 'q9'.
%
% Live monitoring extension:
%   Set params.visualEnable=true to display a 2x3 monitoring dashboard with
%   particles, density, velocity field, mean profile + fit, and time-series
%   diagnostics (kBT, low-k density, div(u), center-wall velocity, nu_eff).
%   A "Stop run" button is added to the figure for graceful early stopping.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);
params = apply_method_params(params);

rng(params.seed);
if logical(get_field(params, 'visualEnable', false))
    reset_visual_stop_request(params);
end
[state, initInfo] = projection_initialize_particles_poiseuille(params);
Np = initInfo.Np;

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
yCenters = ((0:params.Ny-1).' + 0.5) * params.Ly / params.Ny;

UxProfiles = nan(params.Ny, nSamples);
UyProfiles = nan(params.Ny, nSamples);
NProfiles = nan(params.Ny, nSamples);
diagRows = repmat(empty_sample_diag(), nSamples, 1);
fitRows = repmat(empty_fit_diag(), nSamples, 1);

isamp = 0;
actualLastStep = 0;
lastProgressPrint = -inf;
wallClockTic = tic;
stopRequested = false;
stopReason = '';
lastStepDiag = struct();
lastSampleDiag = [];
lastSampleGrid = [];

for it = 0:params.nSteps
    sampledThisIter = false;
    if mod(it, params.sampleEvery) == 0
        sampledThisIter = true;
        isamp = isamp + 1;
        [sdiag, G] = sample_poiseuille_diagnostics(state, params, it, lastStepDiag);
        diagRows(isamp) = sdiag;
        UxProfiles(:, isamp) = mean(G.Ux, 1, 'omitnan').';
        UyProfiles(:, isamp) = mean(G.Uy, 1, 'omitnan').';
        NProfiles(:, isamp) = mean(G.N, 1, 'omitnan').';
        partialOut = build_partial_out(params, state, yCenters, UxProfiles, UyProfiles, NProfiles, diagRows, isamp);
        fitRows(isamp) = compute_live_fit(partialOut, params, it);
        lastSampleDiag = sdiag;
        lastSampleGrid = G;
    end

    if logical(params.visualEnable) && sampledThisIter && (mod(it, params.visualEvery) == 0 || it == params.nSteps)
        partialOut = build_partial_out(params, state, yCenters, UxProfiles, UyProfiles, NProfiles, diagRows, isamp);
        partialOut.fitTable = struct2table(fitRows(1:isamp));
        partialOut.viscosity = fit_struct_to_summary(fitRows(isamp), partialOut, params);
        keepRunning = projection_poiseuille_visualize_frame(state, lastSampleGrid, lastSampleDiag, partialOut, params, it);
        if ~keepRunning
            stopRequested = true;
            stopReason = 'visualStopRequested';
        end
    end

    if stopRequested
        actualLastStep = it;
        break;
    end

    if it == params.nSteps
        actualLastStep = it;
        break;
    end

    nextStep = it + 1;
    stepParams = params;
    % Keep the projection block in its fast/minimal diagnostic mode. The
    % runner samples physical diagnostics itself at sampleEvery.
    stepParams.computeDiagnostics = false;

    switch lower(strrep(char(string(params.method)), '-', '_'))
        case {'classic','src','srd'}
            [state, lastStepDiag] = mpcd_step_classic_poiseuille(state, stepParams);
        case {'q6','projection','velocity_projection'}
            [state, lastStepDiag] = mpcd_step_projection_poiseuille_q9(state, stepParams);
        case {'q9','mass_flux','mass_flux_projection'}
            [state, lastStepDiag] = mpcd_step_projection_poiseuille_q9(state, stepParams);
        otherwise
            error('Unknown method: %s', char(string(params.method)));
    end
    actualLastStep = nextStep;

    willProgress = params.progressEvery > 0 && mod(nextStep, params.progressEvery) == 0 && nextStep ~= lastProgressPrint;
    if willProgress
        lastProgressPrint = nextStep;
        tmp = sample_poiseuille_diagnostics(state, params, nextStep, lastStepDiag);
        fprintf('  Poiseuille %s step %d/%d, t=%.4g, meanUx=%.5g, center-wall=%.5g, kBT=%.4g, lowK=%.3e, elapsed=%.1fs\n', ...
            upper(char(string(params.method))), nextStep, params.nSteps, nextStep*params.dt, ...
            tmp.meanUx, tmp.centerMinusWall, tmp.kBTCell, tmp.lowKDensityEnergy, toc(wallClockTic));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClockTic) > params.maxWallClockSeconds
        warning('run_projection_poiseuille_medium64_demo:TimeLimit', ...
            'Stopping early at step %d/%d after %.1f s because maxWallClockSeconds=%.1f.', ...
            actualLastStep, params.nSteps, toc(wallClockTic), params.maxWallClockSeconds);
        stopRequested = true;
        stopReason = 'maxWallClockSeconds';
    end
end

elapsedWallClock = toc(wallClockTic);
if actualLastStep == 0
    actualLastStep = min(params.nSteps, max([diagRows.step]));
end

valid = isfinite([diagRows.step]).';
diagRows = diagRows(valid);
fitRows = fitRows(valid);
UxProfiles = UxProfiles(:, valid);
UyProfiles = UyProfiles(:, valid);
NProfiles = NProfiles(:, valid);

diagTable = struct2table(diagRows);
fitTable = struct2table(fitRows);
sampleTimes = diagTable.t;
sampleSteps = diagTable.step;

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
out.diagTable = diagTable;
out.fitTable = fitTable;
out.actualLastStep = actualLastStep;
out.stoppedEarly = actualLastStep < params.nSteps;
out.stopReason = stopReason;
out.elapsedWallClock = elapsedWallClock;
fitArgs = poiseuille_fit_args(params);
out.viscosity = analyze_projection_poiseuille_viscosity(out, fitArgs{:});
out.summary = summarize_poiseuille_run(out);
% Two explicit acceleration diagnostics are stored/printed to avoid
% confusing the fit-window/global trend with the recent stationarity trend.
out.accelerationRecent = poiseuille_acceleration_diagnostics(out, params.accelerationWindowTime);
out.accelerationLong = poiseuille_acceleration_diagnostics(out, min(10, max(out.sampleTimes)-min(out.sampleTimes)));

fprintf('\n=== run_projection_poiseuille_medium64_demo (%s) ===\n', upper(char(string(params.method))));
fprintf('Np/grid/gamma                : %d / %d x %d / %.6g\n', Np, params.Nx, params.Ny, params.gamma);
fprintf('dt, kBT, bodyForceX          : %.6g / %.6g / %.6g\n', params.dt, params.kBT, params.bodyForceX);
fprintf('wallModeY                    : %s\n', char(string(params.wallModeY)));
fprintf('steps completed              : %d / %d in %.2f s\n', out.actualLastStep, params.nSteps, out.elapsedWallClock);
fprintf('mean/final center-wall Ux    : %.12g / %.12g\n', out.summary.meanCenterMinusWall, out.summary.finalCenterMinusWall);
fprintf('mean/final kBT cell          : %.12g / %.12g\n', out.summary.meanKBTCell, out.summary.finalKBTCell);
fprintf('mean/final low-k density     : %.12g / %.12g\n', out.summary.meanLowKDensity, out.summary.finalLowKDensity);
fprintf('nu_eff raw Poiseuille fit     : %.12g  (R2=%.6g, SNR=%.6g)\n', ...
    out.viscosity.nuEffRaw, out.viscosity.R2, out.viscosity.signalToNoise);
if isfield(out.viscosity, 'nuEffScaledNyRef') && isfinite(out.viscosity.nuEffScaledNyRef)
    fprintf('nu_eff scaled @ NyRef=%d      : %.12g  (scale=%.6g, dy=%.6g, dyRef=%.6g)\n', ...
        out.viscosity.nuReferenceNy, out.viscosity.nuEffScaledNyRef, ...
        out.viscosity.nuScaleFactorNyRef, out.viscosity.dy, out.viscosity.dyRef);
end
if isfield(out.viscosity, 'nuEffCellY') && isfinite(out.viscosity.nuEffCellY)
    fprintf('nu_eff cell-y units           : %.12g  (= raw/dy^2)\n', out.viscosity.nuEffCellY);
end
accRecent = out.accelerationRecent;
if isfinite(accRecent.ratioAccel)
    fprintf('d<ux>/dt / bodyForceX recent  : %.12g  (d<ux>/dt=%.12g, window %.4g--%.4g, windowTime=%.4g)\n', ...
        accRecent.ratioAccel, accRecent.slopeMeanUx, accRecent.t0, accRecent.t1, params.accelerationWindowTime);
end
accLong = out.accelerationLong;
if isfinite(accLong.ratioAccel)
    fprintf('d<ux>/dt / bodyForceX long    : %.12g  (d<ux>/dt=%.12g, window %.4g--%.4g)\n', ...
        accLong.ratioAccel, accLong.slopeMeanUx, accLong.t0, accLong.t1);
end
end

function params = set_default_params(params)
% Taylor--Green medium64 matched SRC parameters.
params = set_default(params,'method','q9');
params = set_default(params,'Lx',1.0); params = set_default(params,'Ly',1.0);
params = set_default(params,'Nx',64); params = set_default(params,'Ny',64);
params = set_default(params,'gamma',20); params = set_default(params,'seed',11);
params = set_default(params,'dt',1.0e-3); params = set_default(params,'kBT',0.01);
params = set_default(params,'alphaDeg',90);
params = set_default(params,'bodyForceX',0.01); params = set_default(params,'bodyForceY',0.0);
params = set_default(params,'wallModeY','bounceback');
params = set_default(params,'Ubottom',0.0); params = set_default(params,'Utop',0.0);
params = set_default(params,'wallSigma',sqrt(params.kBT));
params = set_default(params,'useRandomGridShiftX',true);
params = set_default(params,'useRandomGridShiftY',false);

% Wall virtual particles for SRC/MPCD no-slip wall coupling. Disabled by
% default for backward compatibility; validation scripts enable it.  The
% shifted_solid_fraction mode is the physically preferred mode and is
% compatible with future stencil/OpenMP implementations.
params = set_default(params,'wallVirtualParticlesEnable',false);
params = set_default(params,'wallVirtualParticlesGeometryMode','shifted_solid_fraction');
params = set_default(params,'wallVirtualParticlesForceRandomShiftY',true);
params = set_default(params,'wallVirtualParticlesCells',1);
params = set_default(params,'wallVirtualParticlesDensityFactor',1.0);
params = set_default(params,'wallVirtualParticlesPerCell',[]);
params = set_default(params,'wallVirtualParticlesThermal',true);
params = set_default(params,'wallVirtualParticlesKBT',params.kBT);
params = set_default(params,'wallVirtualParticlesStochasticCount',true);
params = set_default(params,'wallVirtualParticlesIncludeBottom',true);
params = set_default(params,'wallVirtualParticlesIncludeTop',true);

params = set_default(params,'initialPopulationMode','exact_per_cell');
params = set_default(params,'initialVelocityZeroGlobalMean',true);
params = set_default(params,'initialMeanVelocityX',0.0);
params = set_default(params,'initialMeanVelocityY',0.0);
params = set_default(params,'initialPoiseuilleProfileEnable',false);
params = set_default(params,'initialPoiseuilleNuGuess',0.03);
params = set_default(params,'initialPoiseuilleScale',1.0);
params = set_default(params,'initialPoiseuilleSubtractMean',false);

% Q6/Q9 settings copied from the current forced TG validation.
params = set_default(params,'projectionEnable',true);
params = set_default(params,'projectionStrength',1.0);
params = set_default(params,'projectionInterpolationMethod','nearest');
params = set_default(params,'massFluxProjectionMode','relax_to_uniform_lowk');
params = set_default(params,'massFluxProjectionStrength',1.0);
params = set_default(params,'massFluxDensityRelaxationBeta',5e-4);
params = set_default(params,'massFluxApplyAfterVelocityProjection',true);
params = set_default(params,'massFluxProjectionOperator','general_bc');
params = set_default(params,'massFluxTargetFilter','elliptic_lowpass');
params = set_default(params,'massFluxLowKMaxIndex',2);
params = set_default(params,'massFluxProjectionRegularization',1e-12);
params = set_default(params,'massFluxMinCellCount',1.0);
params = set_default(params,'massFluxFinalVelocityProjectionCleanup',true);
params = set_default(params,'massFluxFinalVelocityProjectionStrength',0.5);
params = set_default(params,'projectionMomentumCorrectionEnable',true);
params = set_default(params,'projectionMomentumCorrectionMode','particle_global_exact');

% Corrected common thermostat.
params = set_default(params,'thermostatAfterStep',true);
params = set_default(params,'thermostatAfterProjection',true);
params = set_default(params,'thermostatTargetKBT',params.kBT);
params = set_default(params,'thermostatStrength',1.0);
params = set_default(params,'thermostatMinParticlesPerCell',2);
params = set_default(params,'thermostatMaxScale',10.0);

params = set_default(params,'nSteps',5000);
params = set_default(params,'sampleEvery',100);
params = set_default(params,'progressEvery',1000);
params = set_default(params,'computeDiagnostics',false);
params = set_default(params,'fitStartFraction',0.5);
params = set_default(params,'fitStartTime',[]);
params = set_default(params,'fitWindowTime',[]);
params = set_default(params,'fitWindowFraction',[]);
params = set_default(params,'fitLastNProfiles',[]);
params = set_default(params,'excludeWallCellsFit',3);
params = set_default(params,'accelerationWindowTime',5.0);

% Poiseuille viscosity scaling. nuEff remains the raw fit in y-coordinate
% units; nuEffScaledNyRef = nuEffRaw*(Ny/NyRef)^2 for fixed Ly.
params = set_default(params,'poiseuilleNuScalingEnable',true);
params = set_default(params,'poiseuilleNuReferenceNy',64);
params = set_default(params,'poiseuillePhysicalLy',params.Ly);
params = set_default(params,'poiseuilleNuScalingMode','Ny_over_reference');

params = set_default(params,'lowKMaxIndex',2);
params = set_default(params,'maxWallClockSeconds',Inf);

% Live monitoring / visualization.
params = set_default(params,'visualEnable',true);
params = set_default(params,'visualEvery',max(params.sampleEvery, 500));
params = set_default(params,'visualPause',0.001);
params = set_default(params,'visualFigureId',510);
params = set_default(params,'visualFigureName','Poiseuille live monitoring');
params = set_default(params,'visualTitleSuffix','');
params = set_default(params,'visualMaxParticles',6000);
params = set_default(params,'visualParticleSize',5);
params = set_default(params,'visualDensityCLim',0.5);
params = set_default(params,'visualUxCLim',NaN);
params = set_default(params,'visualLowKFloor',1e-12);
params = set_default(params,'visualQuiverStrideX',max(2, round(params.Nx/18)));
params = set_default(params,'visualQuiverStrideY',max(2, round(params.Ny/12)));
params = set_default(params,'visualQuiverScale',1.0);
params = set_default(params,'visualSaveFrames',false);
params = set_default(params,'visualFrameDir','poiseuille_live_frames');
params = set_default(params,'visualFramePrefix','poiseuille');
params = set_default(params,'visualFrameResolution',120);
params = set_default(params,'visualStatusTrendWindow',5);
params = set_default(params,'visualSteadySlopeThreshold',5e-4);
params = set_default(params,'visualSteadyLowKGrowthThreshold',2.0);
params = set_default(params,'visualDangerKBT',2.0*params.kBT);
params = set_default(params,'visualDangerLowK',1e-2);
params = set_default(params,'visualDangerRmsDivU',1.0);
params = set_default(params,'visualDangerOutBand',0.8);
params = set_default(params,'visualReferenceNu',NaN);
params = set_default(params,'visualUseScaledNu',true);
params = set_default(params,'visualAccelerationWindowTime',5.0);
end

function params = apply_method_params(params)
method = lower(strrep(char(string(params.method)), '-', '_'));
switch method
    case {'classic','src','srd'}
        params.method = 'classic';
        params.projectionEnable = false;
        params.massFluxProjectionMode = 'off';
        params.massFluxFinalVelocityProjectionCleanup = false;
        params.thermostatAfterStep = true;
        params.thermostatAfterProjection = false;
    case {'q6','projection','velocity_projection'}
        params.method = 'q6';
        params.projectionEnable = true;
        params.projectionStrength = 1.0;
        params.massFluxProjectionMode = 'off';
        params.massFluxFinalVelocityProjectionCleanup = false;
        params.thermostatAfterStep = true;
        params.thermostatAfterProjection = true;
    case {'q9','mass_flux','mass_flux_projection'}
        params.method = 'q9';
        params.projectionEnable = true;
        params.projectionStrength = 1.0;
        params.massFluxProjectionMode = 'relax_to_uniform_lowk';
        params.massFluxProjectionStrength = 1.0;
        params.massFluxApplyAfterVelocityProjection = true;
        if ~isfield(params,'massFluxProjectionOperator') || isempty(params.massFluxProjectionOperator)
            params.massFluxProjectionOperator = 'general_bc';
        end
        if ~isfield(params,'massFluxTargetFilter') || isempty(params.massFluxTargetFilter)
            params.massFluxTargetFilter = 'elliptic_lowpass';
        end
        if ~isfield(params,'massFluxFinalVelocityProjectionCleanup') || isempty(params.massFluxFinalVelocityProjectionCleanup)
            params.massFluxFinalVelocityProjectionCleanup = true;
        end
        params.thermostatAfterStep = true;
        params.thermostatAfterProjection = true;
    otherwise
        error('Unknown method: %s', char(string(params.method)));
end
end

function [diag, G] = sample_poiseuille_diagnostics(state, params, it, lastStepDiag)
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
pop = projection_population_diagnostics(state.x, params, 'periodicX', true, 'periodicY', false);
therm = projection_thermal_diagnostics(state.x, state.v, params, 'periodicX', true, 'periodicY', false, 'minCount', 1);
projDiag = projection_project_grid_periodic_x_neumann_y(G.Ux, G.Uy, params);

UxY = mean(G.Ux, 1, 'omitnan').';
yCenters = ((0:params.Ny-1).' + 0.5) * params.Ly / params.Ny;
centerVelocity = interp1(yCenters, UxY, params.Ly/2, 'linear', 'extrap');
wallMeanVelocity = mean([UxY(1), UxY(end)], 'omitnan');
centerMinusWall = centerVelocity - wallMeanVelocity;

rhoRel = double(G.N) / params.gamma - 1;
[lowKE, totalE] = low_k_density_energy(rhoRel, params.lowKMaxIndex);

diag = empty_sample_diag();
diag.step = it;
diag.t = it * params.dt;
diag.meanUx = mean(state.v(:,1), 'omitnan');
diag.meanUy = mean(state.v(:,2), 'omitnan');
diag.centerVelocity = centerVelocity;
diag.wallMeanVelocity = wallMeanVelocity;
diag.centerMinusWall = centerMinusWall;
diag.kBTCell = therm.kBTCellRelative;
diag.kBTLocalMean = therm.localKBTMean;
diag.kBTGlobal = therm.kBTGlobal;
diag.hydroKE = therm.hydroKineticEnergy;
diag.thermalKE = therm.thermalKineticEnergy;
diag.popStd = pop.stdN;
diag.popOutBand = pop.outBandFraction;
diag.nEmptyCells = pop.nEmptyCells;
diag.densityRelRMS = sqrt(mean(rhoRel(:).^2, 'omitnan'));
diag.lowKDensityEnergy = lowKE;
diag.totalDensityEnergy = totalE;
diag.lowKDensityFraction = lowKE / max(totalE, eps);
diag.rmsDivU = projDiag.rmsDivBefore;
diag.maxAbsDivU = projDiag.maxAbsDivBefore;
diag.momentumRawDV = get_field(lastStepDiag, 'projectionMomentumTotalRawMeanDVNorm', NaN);
diag.momentumResidualDV = get_field(lastStepDiag, 'projectionMomentumTotalResidualMeanDVNorm', NaN);
diag.wallHitsBottom = get_nested_field(lastStepDiag, {'wallInfo','nBot'}, NaN);
diag.wallHitsTop = get_nested_field(lastStepDiag, {'wallInfo','nTop'}, NaN);
end

function diag = empty_sample_diag()
diag = struct('step', NaN, 't', NaN, ...
    'meanUx', NaN, 'meanUy', NaN, ...
    'centerVelocity', NaN, 'wallMeanVelocity', NaN, 'centerMinusWall', NaN, ...
    'kBTCell', NaN, 'kBTLocalMean', NaN, 'kBTGlobal', NaN, ...
    'hydroKE', NaN, 'thermalKE', NaN, ...
    'popStd', NaN, 'popOutBand', NaN, 'nEmptyCells', NaN, ...
    'densityRelRMS', NaN, 'lowKDensityEnergy', NaN, 'totalDensityEnergy', NaN, 'lowKDensityFraction', NaN, ...
    'rmsDivU', NaN, 'maxAbsDivU', NaN, ...
    'momentumRawDV', NaN, 'momentumResidualDV', NaN, ...
    'wallHitsBottom', NaN, 'wallHitsTop', NaN);
end

function fitDiag = empty_fit_diag()
fitDiag = struct('step', NaN, 't', NaN, ...
    'tMin', NaN, 'tMax', NaN, 'nProfiles', NaN, ...
    'nuEff', NaN, 'nuEffRaw', NaN, 'nuEffScaledNyRef', NaN, ...
    'nuEffCellY', NaN, 'nuScaleFactorNyRef', NaN, 'nuReferenceNy', NaN, 'dy', NaN, ...
    'R2', NaN, 'signalToNoise', NaN, ...
    'centerMinusWall', NaN, 'physicalCandidate', NaN);
end

function partialOut = build_partial_out(params, state, yCenters, UxProfiles, UyProfiles, NProfiles, diagRows, isamp)
valid = 1:isamp;
partialOut = struct();
partialOut.params = params;
partialOut.state = state;
partialOut.yCenters = yCenters;
partialOut.UxProfiles = UxProfiles(:, valid);
partialOut.UyProfiles = UyProfiles(:, valid);
partialOut.NProfiles = NProfiles(:, valid);
partialOut.diagTable = struct2table(diagRows(valid));
partialOut.sampleTimes = partialOut.diagTable.t;
partialOut.sampleSteps = partialOut.diagTable.step;
end

function fitDiag = compute_live_fit(partialOut, params, it)
fitDiag = empty_fit_diag();
fitDiag.step = it;
fitDiag.t = it * params.dt;
if isempty(partialOut.diagTable) || height(partialOut.diagTable) < 3
    return;
end
try
    fitArgs = poiseuille_fit_args(params);
    visc = analyze_projection_poiseuille_viscosity(partialOut, fitArgs{:});
    fitDiag.tMin = visc.tMin;
    fitDiag.tMax = visc.tMax;
    fitDiag.nProfiles = visc.nProfiles;
    fitDiag.nuEff = visc.nuEff;
    fitDiag.nuEffRaw = visc.nuEffRaw;
    fitDiag.nuEffScaledNyRef = visc.nuEffScaledNyRef;
    fitDiag.nuEffCellY = visc.nuEffCellY;
    fitDiag.nuScaleFactorNyRef = visc.nuScaleFactorNyRef;
    fitDiag.nuReferenceNy = visc.nuReferenceNy;
    fitDiag.dy = visc.dy;
    fitDiag.R2 = visc.R2;
    fitDiag.signalToNoise = visc.signalToNoise;
    fitDiag.centerMinusWall = visc.centerMinusWall;
    fitDiag.physicalCandidate = double(visc.physicalCandidate);
catch
    % Keep NaNs until enough samples accumulate.
end
end

function visc = fit_struct_to_summary(fitDiag, partialOut, params)
% Construct a viscosity-like summary for visualization from the latest fit
% row. If the row is invalid, try to recompute on the fly.
visc = struct('tMin',NaN,'tMax',NaN,'nProfiles',NaN,'nuEff',NaN,'R2',NaN, ...
    'nuEffRaw',NaN,'nuEffScaledNyRef',NaN,'nuEffCellY',NaN, ...
    'nuScaleFactorNyRef',NaN,'nuReferenceNy',NaN,'dy',NaN, ...
    'signalToNoise',NaN,'centerMinusWall',NaN,'UxMean',nan(size(partialOut.yCenters)), ...
    'UxFit',nan(size(partialOut.yCenters)),'y',partialOut.yCenters,'fitMask',true(size(partialOut.yCenters)));
if isfinite(fitDiag.nuEff)
    try
        fitArgs = poiseuille_fit_args(params);
        tmp = analyze_projection_poiseuille_viscosity(partialOut, fitArgs{:});
        visc = tmp;
    catch
        visc.tMin = fitDiag.tMin;
        visc.tMax = fitDiag.tMax;
        visc.nProfiles = fitDiag.nProfiles;
        visc.nuEff = fitDiag.nuEff;
        visc.nuEffRaw = get_field(fitDiag, 'nuEffRaw', fitDiag.nuEff);
        visc.nuEffScaledNyRef = get_field(fitDiag, 'nuEffScaledNyRef', NaN);
        visc.nuEffCellY = get_field(fitDiag, 'nuEffCellY', NaN);
        visc.nuScaleFactorNyRef = get_field(fitDiag, 'nuScaleFactorNyRef', NaN);
        visc.nuReferenceNy = get_field(fitDiag, 'nuReferenceNy', NaN);
        visc.dy = get_field(fitDiag, 'dy', NaN);
        visc.R2 = fitDiag.R2;
        visc.signalToNoise = fitDiag.signalToNoise;
        visc.centerMinusWall = fitDiag.centerMinusWall;
    end
else
    try
        fitArgs = poiseuille_fit_args(params);
        visc = analyze_projection_poiseuille_viscosity(partialOut, fitArgs{:});
    catch
    end
end
end

function args = poiseuille_fit_args(params)
args = {'excludeWallCells', params.excludeWallCellsFit, ...
    'fitStartFraction', params.fitStartFraction};
if isfield(params, 'fitStartTime') && ~isempty(params.fitStartTime)
    args = [args, {'fitStartTime', params.fitStartTime}];
end
if isfield(params, 'fitWindowTime') && ~isempty(params.fitWindowTime)
    args = [args, {'fitWindowTime', params.fitWindowTime}];
end
if isfield(params, 'fitWindowFraction') && ~isempty(params.fitWindowFraction)
    args = [args, {'fitWindowFraction', params.fitWindowFraction}];
end
if isfield(params, 'fitLastNProfiles') && ~isempty(params.fitLastNProfiles)
    args = [args, {'fitLastNProfiles', params.fitLastNProfiles}];
end
if isfield(params, 'poiseuilleNuScalingEnable') && ~isempty(params.poiseuilleNuScalingEnable)
    args = [args, {'nuScalingEnable', params.poiseuilleNuScalingEnable}];
end
if isfield(params, 'poiseuilleNuReferenceNy') && ~isempty(params.poiseuilleNuReferenceNy)
    args = [args, {'nuReferenceNy', params.poiseuilleNuReferenceNy}];
end
if isfield(params, 'poiseuillePhysicalLy') && ~isempty(params.poiseuillePhysicalLy)
    args = [args, {'physicalLy', params.poiseuillePhysicalLy}];
end
if isfield(params, 'poiseuilleNuScalingMode') && ~isempty(params.poiseuilleNuScalingMode)
    args = [args, {'nuScalingMode', params.poiseuilleNuScalingMode}];
end
end

function [lowKE, totalE] = low_k_density_energy(rhoRel, kmax)
A = double(rhoRel);
A(~isfinite(A)) = 0;
F = fft2(A);
[Nx, Ny] = size(A);
kx = (0:Nx-1).';
ky = 0:Ny-1;
kx = min(kx, Nx-kx);
ky = min(ky, Ny-ky);
mask = (kx <= kmax) & (ky <= kmax);
mask(1,1) = false;
E = abs(F).^2 / numel(F)^2;
lowKE = sum(E(mask));
totalE = sum(E(:), 'omitnan');
end

function summary = summarize_poiseuille_run(out)
T = out.diagTable;
summary = struct();
if isempty(T)
    return;
end
idx = T.t >= out.viscosity.tMin;
if nnz(idx) < 1
    idx = true(height(T), 1);
end
summary.meanUx = mean(T.meanUx(idx), 'omitnan');
summary.finalMeanUx = T.meanUx(end);
summary.meanCenterVelocity = mean(T.centerVelocity(idx), 'omitnan');
summary.finalCenterVelocity = T.centerVelocity(end);
summary.meanWallVelocity = mean(T.wallMeanVelocity(idx), 'omitnan');
summary.finalWallVelocity = T.wallMeanVelocity(end);
summary.meanCenterMinusWall = mean(T.centerMinusWall(idx), 'omitnan');
summary.finalCenterMinusWall = T.centerMinusWall(end);
summary.meanKBTCell = mean(T.kBTCell(idx), 'omitnan');
summary.finalKBTCell = T.kBTCell(end);
summary.meanLowKDensity = mean(T.lowKDensityEnergy(idx), 'omitnan');
summary.finalLowKDensity = T.lowKDensityEnergy(end);
summary.meanDensityRelRMS = mean(T.densityRelRMS(idx), 'omitnan');
summary.finalDensityRelRMS = T.densityRelRMS(end);
summary.meanPopStd = mean(T.popStd(idx), 'omitnan');
summary.finalPopStd = T.popStd(end);
summary.meanPopOutBand = mean(T.popOutBand(idx), 'omitnan');
summary.finalPopOutBand = T.popOutBand(end);
summary.meanRmsDivU = mean(T.rmsDivU(idx), 'omitnan');
summary.finalRmsDivU = T.rmsDivU(end);
summary.meanMomentumRawDV = mean(T.momentumRawDV(idx), 'omitnan');
summary.meanMomentumResidualDV = mean(T.momentumResidualDV(idx), 'omitnan');
end


function reset_visual_stop_request(params)
%RESET_VISUAL_STOP_REQUEST Clear stale stop flag when reusing a MATLAB figure.
try
    figId = get_field(params, 'visualFigureId', 510);
    fig = figure(figId);
    setappdata(fig, 'stopRequested', false);
catch
    % Non-fatal: visualization will still try to initialize the flag at t=0.
end
end

function v = get_field(S, name, defaultValue)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    v = S.(name);
else
    v = defaultValue;
end
end

function v = get_nested_field(S, names, defaultValue)
v = defaultValue;
if ~isstruct(S)
    return;
end
tmp = S;
for k = 1:numel(names)
    if isstruct(tmp) && isfield(tmp, names{k})
        tmp = tmp.(names{k});
    else
        return;
    end
end
v = tmp;
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end
