function out = run_q9_piston_staircase_eos_validation(varargin)
%RUN_Q9_PISTON_STAIRCASE_EOS_VALIDATION Staircase EOS test for Q6/Q9+virial piston.
%
%   out = run_q9_piston_staircase_eos_validation(...)
%
% This runner treats Q6/Q9 as the fluid model and tests the virial EOS/kick
% closure on a sequence of compression plateaux.  The piston moves from one
% compression level to the next over a short ramp, then holds so that bulk
% diagnostics can be averaged on the plateau.
%
% Typical quick limiter/EOS test:
%   out = run_q9_piston_staircase_eos_validation( ...
%       'compressionLevels', [0 0.0025 0.005 0.01], ...
%       'moveSteps', 50, 'holdSteps', 250, 'averageHoldFraction', 0.5, ...
%       'Kvirial', 0.05, 'virialBeta', 0.1, ...
%       'virialMaxDuFractionThermal', 0.02, 'visualEnable', false);

params = parse_inputs(varargin{:});
params = set_defaults(params);

levels = double(params.compressionLevels(:).');
levels(~isfinite(levels)) = [];
if isempty(levels)
    levels = [0 0.0025 0.005 0.01];
end
levels = min(max(levels, 0), 0.95);
if levels(1) ~= 0
    levels = [0 levels];
end
levels = unique(levels, 'stable');
params.compressionLevels = levels;

params.pistonMotionMode = 'staircase';
params.pistonStaircaseCompressionList = levels;
params.pistonStaircaseInitialHoldSteps = params.initialHoldSteps;
params.pistonStaircaseMoveSteps = params.moveSteps;
params.pistonStaircaseHoldSteps = params.holdSteps;
params.pistonYMin = params.pistonY0 * (1.0 - max(levels));
params.pistonCompressionTarget = max(levels);
params.compressionTarget = max(levels);
params.nSteps = params.initialHoldSteps + max(numel(levels)-1, 0) * (params.moveSteps + params.holdSteps);
if params.nSteps <= 0
    params.nSteps = 1;
end

if isempty(params.outputRoot)
    params.outputRoot = fullfile(pwd, ['q9_piston_staircase_eos_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if ~exist(params.outputRoot, 'dir')
    mkdir(params.outputRoot);
end

runParams = rmfield_if_exists(params, {'compressionLevels','moveSteps','holdSteps','initialHoldSteps','averageHoldFraction','outputRoot','saveOutput'});

logFile = fullfile(params.outputRoot, 'console_log.txt');
diary(logFile);
diary on;
cleanupObj = onCleanup(@() diary('off'));

fprintf('=== Q9 piston staircase EOS validation ===\n');
fprintf('compression levels : %s\n', mat2str(levels, 6));
fprintf('move/hold steps    : %d / %d\n', params.moveSteps, params.holdSteps);
fprintf('average hold frac  : %.3g\n', params.averageHoldFraction);
fprintf('Kvirial/beta       : %.6g / %.6g\n', params.Kvirial, params.virialBeta);
fprintf('limiter frac       : %.6g\n', params.virialMaxDuFractionThermal);
fprintf('nSteps             : %d\n\n', params.nSteps);

baseOut = run_projection_piston_wallvp_demo(runParams);

plateauTable = build_plateau_table(baseOut.diagTable, levels, params.averageHoldFraction);
fitTable = build_fit_table(plateauTable);

out = baseOut;
out.staircase = struct();
out.staircase.levels = levels;
out.staircase.moveSteps = params.moveSteps;
out.staircase.holdSteps = params.holdSteps;
out.staircase.averageHoldFraction = params.averageHoldFraction;
out.staircase.plateauTable = plateauTable;
out.staircase.fitTable = fitTable;
out.staircase.outputRoot = params.outputRoot;

fprintf('\n=== Staircase plateau averages ===\n');
disp(plateauTable);
fprintf('\n=== Staircase EOS fits ===\n');
disp(fitTable);

writetable(plateauTable, fullfile(params.outputRoot, 'staircase_plateaux.csv'));
writetable(fitTable, fullfile(params.outputRoot, 'staircase_fits.csv'));
if params.saveOutput
    save(fullfile(params.outputRoot, 'staircase_output.mat'), 'out', '-v7.3');
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
params = set_default(params, 'compressionLevels', [0 0.0025 0.005 0.01]);
params = set_default(params, 'initialHoldSteps', 0);
params = set_default(params, 'moveSteps', 50);
params = set_default(params, 'holdSteps', 250);
params = set_default(params, 'averageHoldFraction', 0.5);
params = set_default(params, 'sampleEvery', 50);
params = set_default(params, 'progressEvery', 100);
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
params = set_default(params, 'outputRoot', '');
params = set_default(params, 'saveOutput', true);
params = set_default(params, 'maxWallClockSeconds', Inf);
end

function plateauTable = build_plateau_table(T, levels, avgFrac)
vars = {'rhoPhysicalMean','PkinMean','PvirMean','PtotMean','PdriveMean','PkinIdealRatio','PtotIdealRatio', ...
    'stdN','lowKDensity','kBTCell','virialDuRms','virialDuOverThermalRms', ...
    'virialLimitedCellFraction','virialLimitedCellCount','virialLimiterDuMax', ...
    'virialResidualMomentumKickNorm','massFluxBefore','massFluxAfter','massFluxResidual'};
rows = cell(numel(levels), 1);
for j = 1:numel(levels)
    if ismember('pistonPhaseIndex', T.Properties.VariableNames)
        idx = T.pistonPhaseIndex == j;
    else
        idx = abs(T.compression - levels(j)) < 1e-12;
    end
    if ismember('pistonPhaseName', T.Properties.VariableNames)
        phase = string(T.pistonPhaseName);
        if j == 1
            idx = idx & (phase == "staircase_initial_hold" | phase == "staircase_hold" | abs(T.compression - levels(j)) < 1e-12);
        else
            idx = idx & phase == "staircase_hold";
        end
    end
    ids = find(idx);
    if isempty(ids)
        ids = find(abs(T.compression - levels(j)) < 10*eps(max(1, levels(j))));
    end
    if isempty(ids)
        ids = find(abs(T.compression - levels(j)) == min(abs(T.compression - levels(j))), 1, 'last');
    end
    nKeep = max(1, ceil(numel(ids) * min(max(avgFrac, 0), 1)));
    idsAvg = ids(end-nKeep+1:end);
    R = T(idsAvg, :);
    row = table();
    row.levelIndex = j;
    row.targetCompression = levels(j);
    row.nSamples = height(R);
    row.stepMin = min(R.step);
    row.stepMax = max(R.step);
    row.tMin = min(R.t);
    row.tMax = max(R.t);
    row.yTopMean = mean(R.yTop, 'omitnan');
    row.compressionMean = mean(R.compression, 'omitnan');
    for iv = 1:numel(vars)
        name = vars{iv};
        if ismember(name, R.Properties.VariableNames)
            row.(name) = mean(R.(name), 'omitnan');
        else
            row.(name) = NaN;
        end
    end
    rows{j} = row;
end
plateauTable = vertcat(rows{:});
end

function fitTable = build_fit_table(P)
quantities = {'PkinMean','PvirMean','PtotMean','PdriveMean'};
rows = cell(numel(quantities), 1);
for i = 1:numel(quantities)
    q = quantities{i};
    x = P.rhoPhysicalMean;
    y = P.(q);
    ok = isfinite(x) & isfinite(y);
    row = table();
    row.quantity = string(q);
    row.n = nnz(ok);
    if nnz(ok) >= 2 && range(x(ok)) > 0
        c = polyfit(x(ok), y(ok), 1);
        yhat = polyval(c, x(ok));
        ssRes = sum((y(ok)-yhat).^2);
        ssTot = sum((y(ok)-mean(y(ok))).^2);
        row.slope_dPdrho = c(1);
        row.intercept = c(2);
        row.Keff = c(1) * mean(x(ok));
        if ssTot > 0
            row.R2 = 1 - ssRes/ssTot;
        else
            row.R2 = NaN;
        end
        row.Pmin = min(y(ok));
        row.Pmax = max(y(ok));
    else
        row.slope_dPdrho = NaN;
        row.intercept = NaN;
        row.Keff = NaN;
        row.R2 = NaN;
        row.Pmin = NaN;
        row.Pmax = NaN;
    end
    rows{i} = row;
end
fitTable = vertcat(rows{:});
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
