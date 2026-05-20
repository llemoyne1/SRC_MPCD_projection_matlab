function out = run_q9_piston_smooth_wallbulk_validation(varargin)
%RUN_Q9_PISTON_SMOOTH_WALLBULK_VALIDATION Smooth-ramp wall/bulk pressure test.
%
%   out = run_q9_piston_smooth_wallbulk_validation(...)
%
% This validation replaces staircase pressure plateaux by a smooth quasi-static
% piston compression and monitors rolling-window residuals between wall impulse
% pressure and local bulk pressure near the two walls:
%
%   R_top    = P_wall_top^Tw    - <Pkin+Pvir>_topLayer^Tw
%   R_bottom = P_wall_bottom^Tw - <Pkin+Pvir>_bottomLayer^Tw
%   R_sym    = 0.5*(R_top + R_bottom)
%   R_anti   = 0.5*(R_top - R_bottom)
%
% The wall pressure is reconstructed from cumulative wall impulses over the
% rolling window, not from instantaneous sampled pressures.

params = parse_inputs(varargin{:});
params = set_defaults(params);

params.pistonMotionMode = 'smooth_ramp';
params.pistonSmoothRampCompressionFinal = params.compressionFinal;
params.pistonSmoothRampInitialHoldSteps = params.initialHoldSteps;
params.pistonSmoothRampSteps = params.rampSteps;
params.pistonSmoothRampFinalHoldSteps = params.finalHoldSteps;
params.pistonSmoothRampProfile = params.rampProfile;
params.pistonCompressionTarget = params.compressionFinal;
params.compressionTarget = params.compressionFinal;
params.pistonYMin = params.pistonY0 * (1.0 - params.compressionFinal);
params.nSteps = params.initialHoldSteps + params.rampSteps + params.finalHoldSteps;
params.wallBulkResidualEnable = true;
params.wallBulkResidualWindowSteps = params.wallBulkResidualWindowSteps;
params.wallBulkResidualLayerCells = params.wallBulkResidualLayerCells;

if isempty(params.outputRoot)
    params.outputRoot = fullfile(pwd, ['q9_piston_smooth_wallbulk_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if ~exist(params.outputRoot, 'dir')
    mkdir(params.outputRoot);
end

runParams = rmfield_if_exists(params, {'compressionFinal','initialHoldSteps','rampSteps','finalHoldSteps', ...
    'rampProfile','outputRoot','saveOutput','plotEnable'});

logFile = fullfile(params.outputRoot, 'console_log.txt');
diary(logFile);
diary on;
cleanupObj = onCleanup(@() diary('off'));

fprintf('=== Q9 piston smooth wall/bulk validation ===\n');
fprintf('compression final : %.6g\n', params.compressionFinal);
fprintf('initial/ramp/final: %d / %d / %d steps\n', params.initialHoldSteps, params.rampSteps, params.finalHoldSteps);
fprintf('window/layer      : %d steps / %d cells\n', params.wallBulkResidualWindowSteps, params.wallBulkResidualLayerCells);
fprintf('Kvirial/beta      : %.6g / %.6g\n', params.Kvirial, params.virialBeta);
fprintf('nSteps            : %d\n\n', params.nSteps);

baseOut = run_projection_piston_wallvp_demo(runParams);
T = baseOut.diagTable;
summaryTable = build_wallbulk_summary_table(T);

out = baseOut;
out.wallBulk = struct();
out.wallBulk.outputRoot = params.outputRoot;
out.wallBulk.summaryTable = summaryTable;
out.wallBulk.windowSteps = params.wallBulkResidualWindowSteps;
out.wallBulk.layerCells = params.wallBulkResidualLayerCells;
out.wallBulk.compressionFinal = params.compressionFinal;
out.wallBulk.rampSteps = params.rampSteps;

fprintf('\n=== Rolling wall/bulk residual summary ===\n');
disp(summaryTable);

writetable(T, fullfile(params.outputRoot, 'wall_bulk_timeseries.csv'));
writetable(summaryTable, fullfile(params.outputRoot, 'wall_bulk_summary.csv'));
if logical(params.plotEnable)
    try
        make_wallbulk_plots(T, params.outputRoot);
    catch ME
        warning('run_q9_piston_smooth_wallbulk_validation:PlotFailed', ...
            'Could not create wall/bulk plots: %s', ME.message);
    end
end
if params.saveOutput
    save(fullfile(params.outputRoot, 'smooth_wallbulk_output.mat'), 'out', '-v7.3');
end

clear cleanupObj;
diary off;
end

function params = set_defaults(params)
params = set_default(params, 'method', 'q9_full');
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 32);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'Lx', 1.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Ly0', params.Ly);
params = set_default(params, 'pistonY0', params.Ly0);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 0.01);
params = set_default(params, 'seed', 11);
params = set_default(params, 'compressionFinal', 0.01);
params = set_default(params, 'initialHoldSteps', 1000);
params = set_default(params, 'rampSteps', 6500);
params = set_default(params, 'finalHoldSteps', 1500);
params = set_default(params, 'rampProfile', 'smoothstep');
params = set_default(params, 'wallBulkResidualWindowSteps', 1500);
params = set_default(params, 'wallBulkResidualLayerCells', 3);
params = set_default(params, 'sampleEvery', 100);
params = set_default(params, 'progressEvery', 500);
params = set_default(params, 'visualEnable', false);
params = set_default(params, 'storeDensityMaps', false);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'densityFactor', 0.8);
params = set_default(params, 'wallVirtualParticlesDensityFactor', params.densityFactor);
params = set_default(params, 'massFluxProjectionOperator', 'general_bc');
params = set_default(params, 'massFluxTargetFilter', 'elliptic_lowpass');
params = set_default(params, 'massFluxDensityRelaxationBeta', 5e-4);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', false);
params = set_default(params, 'pistonQ9TargetMode', 'reference_gamma');
params = set_default(params, 'pistonQ9FreezeEllipticMetric', false);
params = set_default(params, 'virialDiagnosticsEnable', true);
params = set_default(params, 'virialKickEnable', true);
params = set_default(params, 'Kvirial', 0.01);
params = set_default(params, 'virialBeta', 0.02);
params = set_default(params, 'virialLimiterEnable', true);
params = set_default(params, 'virialMaxDuFractionThermal', 0.025);
params = set_default(params, 'virialDriveTargetMode', 'current_uniform');
params = set_default(params, 'virialRhoEOSRefMode', 'initial_physical_density');
params = set_default(params, 'virialRhoUniformMode', 'reference_gamma_current_volume');
params = set_default(params, 'virialRhoKickMode', 'uniform_now');
params = set_default(params, 'thermostatAfterStep', true);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'plotEnable', true);
params = set_default(params, 'outputRoot', '');
params = set_default(params, 'saveOutput', true);
params = set_default(params, 'maxWallClockSeconds', Inf);
end

function S = build_wallbulk_summary_table(T)
vars = {'compression','rhoPhysicalMean','PkinMean','PvirMean','PtotMean', ...
    'pressureWallSymRollingAverage','PtotLayerSymRollingAverage', ...
    'wallBulkResidualSymPtotRollingAverage','wallBulkResidualAntiPtotRollingAverage', ...
    'meanAbsVy','rmsVy','pistonPowerOnFluid'};
labels = {'all_finite','last_25_percent','final_sample'};
rows = cell(numel(labels),1);
for k = 1:numel(labels)
    switch labels{k}
        case 'all_finite'
            idx = isfinite(column_or_nan(T, 'wallBulkResidualSymPtotRollingAverage'));
        case 'last_25_percent'
            idx0 = isfinite(column_or_nan(T, 'wallBulkResidualSymPtotRollingAverage'));
            ids0 = find(idx0);
            nKeep = max(1, ceil(0.25*numel(ids0)));
            idx = false(height(T),1);
            if ~isempty(ids0), idx(ids0(end-nKeep+1:end)) = true; end
        case 'final_sample'
            idx = false(height(T),1);
            ids0 = find(isfinite(column_or_nan(T, 'wallBulkResidualSymPtotRollingAverage')));
            if ~isempty(ids0), idx(ids0(end)) = true; end
    end
    row = table();
    row.window = string(labels{k});
    row.nSamples = nnz(idx);
    for iv = 1:numel(vars)
        x = column_or_nan(T, vars{iv});
        row.(vars{iv}) = mean(x(idx), 'omitnan');
    end
    rows{k} = row;
end
S = vertcat(rows{:});
end

function make_wallbulk_plots(T, outputRoot)
idx = isfinite(column_or_nan(T, 'wallBulkResidualSymPtotRollingAverage'));
if ~any(idx)
    return;
end
fig = figure('Visible','off');
plot(T.compression(idx), T.wallBulkResidualSymPtotRollingAverage(idx), '-o');
hold on;
plot(T.compression(idx), T.wallBulkResidualAntiPtotRollingAverage(idx), '-x');
yline(0, '--');
xlabel('compression');
ylabel('rolling residual');
legend({'sym: Pwall-Ptot','anti top-bottom'}, 'Location', 'best');
grid on;
saveas(fig, fullfile(outputRoot, 'wall_bulk_residual_vs_compression.png'));
close(fig);

fig = figure('Visible','off');
plot(T.t(idx), T.pressureWallSymRollingAverage(idx), '-o');
hold on;
plot(T.t(idx), T.PtotLayerSymRollingAverage(idx), '-x');
xlabel('time');
ylabel('pressure');
legend({'Pwall sym rolling','Ptot layer sym rolling'}, 'Location', 'best');
grid on;
saveas(fig, fullfile(outputRoot, 'wall_vs_bulk_pressure_timeseries.png'));
close(fig);
end

function x = column_or_nan(T, name)
if istable(T) && ismember(name, T.Properties.VariableNames)
    x = T.(name);
else
    x = nan(height(T), 1);
end
end

function params = parse_inputs(varargin)
if nargin == 1 && isstruct(varargin{1})
    params = varargin{1};
    return;
end
params = struct();
if mod(nargin, 2) ~= 0
    error('Use name-value pairs or a single struct.');
end
for k = 1:2:nargin
    params.(char(varargin{k})) = varargin{k+1};
end
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end

function s = rmfield_if_exists(s, names)
for i = 1:numel(names)
    if isfield(s, names{i})
        s = rmfield(s, names{i});
    end
end
end
