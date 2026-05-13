%% run_q9_piston_eos_cycle.m
% Q9-native slow piston cycle for an effective equation-of-state diagnostic.
%
% The run performs a slow compression-hold-decompression-hold cycle and
% measures three pressure-like quantities:
%
%   Pkin             : kinetic/cell pressure, rho*kBT_cell, from grid moments.
%   Pwall            : mechanical piston pressure from wall impulse.
%   PexcessWallKinetic: residual effective pressure, Pwall - Pkin.
%
% A fourth diagnostic is also stored for interpretation:
%
%   PprojectionWorkRaw = -Delta E_projection_raw / Delta V,
%
% which is only meaningful during ramps where Delta V ~= 0.  It is not used
% as the primary EOS pressure because it is protocol-dependent and can be
% contaminated by thermostat and redistribution of kinetic energy.

clear functions
close all
% clc

%% === Output folder ===
tag = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = ['q9_piston_eos_cycle_' tag];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

diary(fullfile(outputDir, 'console_log.txt'));
diary on
cleanupObj = onCleanup(@() diary('off'));

fprintf('\n=== Q9 piston EOS cycle ===\n');

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

%% === Slow compression / decompression cycle ===
params.pistonMotionMode = 'cycle';
params.pistonY0 = params.Ly;
params.pistonCompressionTarget = 0.05;
params.pistonYMin = params.pistonY0 * (1.0 - params.pistonCompressionTarget);

% Default cycle: slow enough to be close to quasi-static, but still usable in MATLAB.
params.pistonCycleCompressSteps = 20000;
params.pistonCycleHoldSteps = 5000;
params.pistonCycleDecompressSteps = 20000;
params.pistonCycleFinalHoldSteps = 5000;
params.nSteps = params.pistonCycleCompressSteps + params.pistonCycleHoldSteps + ...
    params.pistonCycleDecompressSteps + params.pistonCycleFinalHoldSteps;

params.sampleEvery = 100;
params.progressEvery = 1000;
params.maxWallClockSeconds = Inf;
params.pistonStopOnMin = true;
params.pistonAffineReposition = true;

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
params.smoothWindowSamples = 21;
params.makeFigures = true;

%% === Initialization ===
rng(params.seed);
[state, initInfo] = projection_initialize_particles(params);
state.piston = struct('yTop', params.pistonY0, 'yPrev', params.pistonY0, ...
    'Up', 0.0, 'Ly0', params.Ly, 'compression', 0.0, ...
    'activeHeight', params.pistonY0, 'phase', 'initial', 'phaseIndex', 0);

fprintf('Np                         : %d\n', initInfo.Np);
fprintf('grid                       : %d x %d, gamma=%g\n', params.Nx, params.Ny, params.gamma);
fprintf('cycle steps C/H/D/H        : %d / %d / %d / %d\n', ...
    params.pistonCycleCompressSteps, params.pistonCycleHoldSteps, ...
    params.pistonCycleDecompressSteps, params.pistonCycleFinalHoldSteps);
fprintf('nSteps/sampleEvery         : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf('target compression         : %.6g\n', params.pistonCompressionTarget);
fprintf('Q9 beta/lowK/cleanup       : %.6g / %d / %.6g\n', ...
    params.massFluxDensityRelaxationBeta, params.massFluxLowKMaxIndex, ...
    params.massFluxFinalVelocityProjectionStrength);

%% === Time loop ===
nSamples = floor(params.nSteps / params.sampleEvery) + 1;
rows = repmat(empty_row(), nSamples, 1);
isamp = 0;
acc = empty_accumulator();
actualLastStep = 0;
wallClock = tic;

% Initial sample.
isamp = isamp + 1;
rows(isamp) = make_sample_row(0, 0.0, state, params, [], acc);
acc = empty_accumulator();

for it = 1:params.nSteps
    stepParams = params;
    stepParams.computeDiagnostics = true;
    [state, stepDiag] = mpcd_step_projection_piston(state, stepParams, it);
    actualLastStep = it;

    acc = accumulate_step(acc, stepDiag, params);

    if mod(it, params.sampleEvery) == 0 || it == params.nSteps
        isamp = isamp + 1;
        rows(isamp) = make_sample_row(it, it*params.dt, state, params, stepDiag, acc);
        acc = empty_accumulator();
    end

    if params.progressEvery > 0 && mod(it, params.progressEvery) == 0
        p = state.piston;
        fprintf('  eos piston step %d/%d, t=%.4g, phase=%s, comp=%.4g, yTop=%.6g, Pwall=%.4g, Pkin=%.4g, elapsed=%.1fs\n', ...
            it, params.nSteps, it*params.dt, char(p.phase), p.compression, p.yTop, ...
            finite_or_nan(stepDiag.wallInfo.pressureTopWall), rows(isamp).Pkin, toc(wallClock));
    end

    if isfinite(params.maxWallClockSeconds) && toc(wallClock) > params.maxWallClockSeconds
        warning('run_q9_piston_eos_cycle:TimeLimit', ...
            'Stopping at step %d/%d after %.1fs.', it, params.nSteps, toc(wallClock));
        break;
    end
end

elapsedScript = toc(wallClock);
rows = rows(1:isamp);
T = struct2table(rows);
T = add_smoothed_pressures(T, params.smoothWindowSamples);
summary = summarize_eos_table(T, params, actualLastStep, elapsedScript);

%% === Save ===
matFile = fullfile(outputDir, 'q9_piston_eos_cycle.mat');
csvFile = fullfile(outputDir, 'q9_piston_eos_cycle_timeseries.csv');
summaryCsvFile = fullfile(outputDir, 'q9_piston_eos_cycle_summary.csv');
txtFile = fullfile(outputDir, 'q9_piston_eos_cycle_summary.txt');

save(matFile, 'params', 'initInfo', 'T', 'summary', '-v7.3');
writetable(T, csvFile);
writetable(summary, summaryCsvFile);
write_summary_text(txtFile, params, initInfo, summary, matFile, csvFile, summaryCsvFile, txtFile);

fprintf('\nTimeseries written to:\n%s\n', csvFile);
fprintf('Summary written to:\n%s\n', txtFile);
fprintf('MAT written to:\n%s\n', matFile);

if params.makeFigures
    make_eos_figures(T, summary, outputDir);
end

diary off

%% ========================================================================
function row = empty_row()
row = struct();
row.step = NaN;
row.t = NaN;
row.phase = "";
row.phaseIndex = NaN;
row.yTop = NaN;
row.compression = NaN;
row.volume = NaN;
row.volumeRatio = NaN;
row.rhoMean = NaN;
row.rhoRatio = NaN;
row.Pkin = NaN;
row.Pwall = NaN;
row.PexcessWallKinetic = NaN;
row.PprojectionWorkRaw = NaN;
row.PprojectionWorkThermostatted = NaN;
row.wallImpulseTopY = NaN;
row.wallHitsTop = NaN;
row.dVInterval = NaN;
row.dVdtMean = NaN;
row.projectionEnergyRaw = NaN;
row.projectionEnergyThermostatted = NaN;
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
acc.dV = 0.0;
acc.projEnergyRaw = 0.0;
acc.projEnergyThermo = 0.0;
end

function acc = accumulate_step(acc, stepDiag, params)
acc.n = acc.n + 1;
acc.dt = acc.dt + params.dt;
W = getf(stepDiag, 'wallInfo', struct());
acc.wallImpulseTopY = acc.wallImpulseTopY + getf(W, 'impulseOnTopWallY', signed_top_impulse_fallback(W));
acc.wallHitsTop = acc.wallHitsTop + getf(W, 'nTop', 0);
P = getf(stepDiag, 'piston', struct());
acc.dV = acc.dV + getf(P, 'dV', 0.0);
acc.projEnergyRaw = acc.projEnergyRaw + finite_zero(getf(stepDiag, 'totalKEAfterProjectionRaw') - getf(stepDiag, 'totalKEBeforeProjection'));
acc.projEnergyThermo = acc.projEnergyThermo + finite_zero(getf(stepDiag, 'totalKEAfterProjection') - getf(stepDiag, 'totalKEBeforeProjection'));
end

function row = make_sample_row(it, t, state, params, stepDiag, acc)
activeParams = params;
if isfield(state, 'piston') && isfield(state.piston, 'yTop')
    activeParams.Ly = state.piston.yTop;
    yTop = state.piston.yTop;
    compression = state.piston.compression;
    phase = string(state.piston.phase);
    phaseIndex = state.piston.phaseIndex;
else
    yTop = params.Ly;
    compression = 0.0;
    phase = "initial";
    phaseIndex = 0;
end

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
rho0 = params.Nx * params.Ny * params.gamma / max(params.Lx * params.Ly, eps);

if acc.n > 0 && acc.dt > 0
    Pwall = acc.wallImpulseTopY / max(acc.dt * params.Lx, eps);
    dVdt = acc.dV / max(acc.dt, eps);
    if abs(acc.dV) > eps(volume)
        PprojRaw = -acc.projEnergyRaw / acc.dV;
        PprojThermo = -acc.projEnergyThermo / acc.dV;
    else
        PprojRaw = NaN;
        PprojThermo = NaN;
    end
else
    Pwall = NaN;
    dVdt = NaN;
    PprojRaw = NaN;
    PprojThermo = NaN;
end

row = empty_row();
row.step = it;
row.t = t;
row.phase = phase;
row.phaseIndex = phaseIndex;
row.yTop = yTop;
row.compression = compression;
row.volume = volume;
row.volumeRatio = volume / max(params.Lx * params.Ly, eps);
row.rhoMean = rhoMean;
row.rhoRatio = rhoMean / max(rho0, eps);
row.Pkin = Pkin;
row.Pwall = Pwall;
row.PexcessWallKinetic = Pwall - Pkin;
row.PprojectionWorkRaw = PprojRaw;
row.PprojectionWorkThermostatted = PprojThermo;
row.wallImpulseTopY = acc.wallImpulseTopY;
row.wallHitsTop = acc.wallHitsTop;
row.dVInterval = acc.dV;
row.dVdtMean = dVdt;
row.projectionEnergyRaw = acc.projEnergyRaw;
row.projectionEnergyThermostatted = acc.projEnergyThermo;
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
T.PprojectionWorkRawSmooth = movmean(T.PprojectionWorkRaw, w, 'omitnan');
T.PprojectionWorkThermostattedSmooth = movmean(T.PprojectionWorkThermostatted, w, 'omitnan');
end

function summary = summarize_eos_table(T, params, actualLastStep, elapsedScript)
phases = ["compression"; "hold_compressed"; "decompression"; "hold_relaxed"];
summary = table();
for i = 1:numel(phases)
    ph = phases(i);
    idx = T.phase == ph;
    idx = idx & isfinite(T.rhoMean) & isfinite(T.PwallSmooth);
    S = table();
    S.phase = ph;
    S.nSamples = nnz(idx);
    S.actualLastStep = actualLastStep;
    S.elapsedScript = elapsedScript;
    S.meanRhoRatio = mean(T.rhoRatio(idx), 'omitnan');
    S.meanCompression = mean(T.compression(idx), 'omitnan');
    S.meanPkin = mean(T.PkinSmooth(idx), 'omitnan');
    S.meanPwall = mean(T.PwallSmooth(idx), 'omitnan');
    S.meanPexcessWallKinetic = mean(T.PexcessWallKineticSmooth(idx), 'omitnan');
    S.meanPprojectionWorkRaw = mean(T.PprojectionWorkRawSmooth(idx), 'omitnan');
    S.meanStdN = mean(T.stdN(idx), 'omitnan');
    S.meanLowKProxy = mean(T.lowKProxy(idx), 'omitnan');
    S.meanDivUAfter = mean(T.rmsDivParticleAfter(idx), 'omitnan');
    S.meanMassFluxResidual = mean(T.rmsMassFluxResidual(idx), 'omitnan');
    S.KeffWallFit = fit_keff(T.rhoMean(idx), T.PwallSmooth(idx));
    S.KeffKineticFit = fit_keff(T.rhoMean(idx), T.PkinSmooth(idx));
    S.KeffExcessFit = fit_keff(T.rhoMean(idx), T.PexcessWallKineticSmooth(idx));
    summary = [summary; S]; %#ok<AGROW>
end
summary.beta = repmat(params.massFluxDensityRelaxationBeta, height(summary), 1);
summary.lowKMaxIndex = repmat(params.massFluxLowKMaxIndex, height(summary), 1);
summary.cleanupStrength = repmat(params.massFluxFinalVelocityProjectionStrength, height(summary), 1);
end

function K = fit_keff(rho, P)
mask = isfinite(rho) & isfinite(P);
rho = rho(mask);
P = P(mask);
if numel(rho) < 5 || max(rho) - min(rho) <= 10*eps(max(abs(rho)))
    K = NaN;
    return;
end
coef = polyfit(rho, P, 1);
rhoBar = mean(rho, 'omitnan');
K = rhoBar * coef(1);
end

function write_summary_text(filename, params, initInfo, summary, matFile, csvFile, summaryCsvFile, txtFile)
fid = fopen(filename, 'w');
fprintf(fid, '=== Q9 piston EOS cycle ===\n');
fprintf(fid, 'seed                              : %d\n', params.seed);
fprintf(fid, 'Np                                : %d\n', initInfo.Np);
fprintf(fid, 'steps/sampleEvery                 : %d / %d\n', params.nSteps, params.sampleEvery);
fprintf(fid, 'grid, gamma                       : %d x %d, %g\n', params.Nx, params.Ny, params.gamma);
fprintf(fid, 'compression target                : %.12g\n', params.pistonCompressionTarget);
fprintf(fid, 'cycle C/H/D/H steps               : %d / %d / %d / %d\n', ...
    params.pistonCycleCompressSteps, params.pistonCycleHoldSteps, ...
    params.pistonCycleDecompressSteps, params.pistonCycleFinalHoldSteps);
fprintf(fid, 'dt, kBT, alphaDeg                 : %.12g, %.12g, %.12g\n', params.dt, params.kBT, params.alphaDeg);
fprintf(fid, '\n--- Q9 parameters ---\n');
fprintf(fid, 'projectionStrength                : %.6g\n', params.projectionStrength);
fprintf(fid, 'massFluxProjectionMode            : %s\n', params.massFluxProjectionMode);
fprintf(fid, 'massFluxDensityRelaxationBeta     : %.6g\n', params.massFluxDensityRelaxationBeta);
fprintf(fid, 'massFluxFinalVelocityProjectionCleanup : %d\n', params.massFluxFinalVelocityProjectionCleanup);
fprintf(fid, 'massFluxFinalVelocityProjectionStrength : %.6g\n', params.massFluxFinalVelocityProjectionStrength);
fprintf(fid, 'massFluxLowKMaxIndex              : %d\n', params.massFluxLowKMaxIndex);
fprintf(fid, '\n--- Pressure definitions ---\n');
fprintf(fid, 'Pkin              = mean_cell(rho_cell * kBT_cell)\n');
fprintf(fid, 'Pwall             = impulse_on_top_wall_y / (Delta_t_sample * Lx)\n');
fprintf(fid, 'PexcessWallKinetic= Pwall - Pkin\n');
fprintf(fid, 'PprojectionWorkRaw= -Delta E_projection_raw / Delta V, only during ramps\n');
fprintf(fid, '\n--- Phase summary ---\n');
for i = 1:height(summary)
    fprintf(fid, '\nphase: %s\n', char(summary.phase(i)));
    fprintf(fid, 'nSamples                         : %d\n', summary.nSamples(i));
    fprintf(fid, 'mean rho/rho0                    : %.12g\n', summary.meanRhoRatio(i));
    fprintf(fid, 'mean compression                 : %.12g\n', summary.meanCompression(i));
    fprintf(fid, 'mean Pkin                        : %.12g\n', summary.meanPkin(i));
    fprintf(fid, 'mean Pwall                       : %.12g\n', summary.meanPwall(i));
    fprintf(fid, 'mean PexcessWallKinetic          : %.12g\n', summary.meanPexcessWallKinetic(i));
    fprintf(fid, 'mean PprojectionWorkRaw          : %.12g\n', summary.meanPprojectionWorkRaw(i));
    fprintf(fid, 'Keff wall fit                    : %.12g\n', summary.KeffWallFit(i));
    fprintf(fid, 'Keff kinetic fit                 : %.12g\n', summary.KeffKineticFit(i));
    fprintf(fid, 'Keff excess fit                  : %.12g\n', summary.KeffExcessFit(i));
    fprintf(fid, 'mean std(N)                      : %.12g\n', summary.meanStdN(i));
    fprintf(fid, 'mean low-k proxy                 : %.12g\n', summary.meanLowKProxy(i));
end
fprintf(fid, '\n--- Files ---\n');
fprintf(fid, '%s\n', matFile);
fprintf(fid, '%s\n', csvFile);
fprintf(fid, '%s\n', summaryCsvFile);
fprintf(fid, '%s\n', txtFile);
fclose(fid);
end

function make_eos_figures(T, summary, outputDir)
fig1 = figure('Name', 'Piston EOS cycle: pressure versus density');
hold on; grid on;
plot(T.rhoRatio, T.PkinSmooth, '-', 'DisplayName', 'Pkin');
plot(T.rhoRatio, T.PwallSmooth, '-', 'DisplayName', 'Pwall');
plot(T.rhoRatio, T.PexcessWallKineticSmooth, '-', 'DisplayName', 'Pwall-Pkin');
xlabel('\rho / \rho_0'); ylabel('pressure-like quantity'); legend('Location', 'best');
title('Effective EOS diagnostics');
saveas(fig1, fullfile(outputDir, 'piston_eos_pressure_vs_density.png'));

fig2 = figure('Name', 'Piston EOS cycle: time series');
tiledlayout(4,1);
nexttile; plot(T.t, T.compression, '-'); grid on; ylabel('compression');
nexttile; plot(T.t, T.PkinSmooth, '-', T.t, T.PwallSmooth, '-'); grid on; ylabel('P'); legend('Pkin','Pwall');
nexttile; plot(T.t, T.PexcessWallKineticSmooth, '-'); grid on; ylabel('Pwall-Pkin');
nexttile; plot(T.t, T.stdN, '-', T.t, T.lowKProxy, '-'); grid on; ylabel('density diag'); xlabel('t'); legend('stdN','low-k proxy');
saveas(fig2, fullfile(outputDir, 'piston_eos_timeseries.png'));

fig3 = figure('Name', 'Piston EOS cycle: fitted Keff');
bar(categorical(summary.phase), [summary.KeffKineticFit, summary.KeffWallFit, summary.KeffExcessFit]);
grid on; ylabel('K_{eff} fit'); legend('kinetic','wall','excess','Location','best');
title('Phase-wise effective compressibility modulus');
saveas(fig3, fullfile(outputDir, 'piston_eos_keff_summary.png'));
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

function v = signed_top_impulse_fallback(W)
if isstruct(W) && isfield(W, 'dPyTopSigned') && isfinite(W.dPyTopSigned)
    v = -W.dPyTopSigned;
elseif isstruct(W) && isfield(W, 'dPyTop') && isfinite(W.dPyTop)
    v = W.dPyTop;
else
    v = 0.0;
end
end

function y = finite_zero(x)
if isempty(x) || ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
    y = 0.0;
else
    y = x;
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
