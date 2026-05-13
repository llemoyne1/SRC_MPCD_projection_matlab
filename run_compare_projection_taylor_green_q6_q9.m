function cmp = run_compare_projection_taylor_green_q6_q9(params)
%RUN_COMPARE_PROJECTION_TAYLOR_GREEN_Q6_Q9 Compare Q6 and Q9 on Taylor-Green.
%
%   cmp = run_compare_projection_taylor_green_q6_q9(params)
%
% Uses the same seed for Q6 and Q9, so both start from the same particle
% population and initial Taylor-Green velocity field.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);

fprintf('\n=== Taylor-Green comparison: Q6 projection-only vs Q9 low-k mass flux ===\n');

paramsQ6 = params;
paramsQ6.massFluxProjectionMode = 'off';
paramsQ6.massFluxDensityRelaxationBeta = 0.0;
paramsQ6.massFluxFinalVelocityProjectionCleanup = false;
paramsQ6.methodName = 'Q6_projection';
paramsQ6.makeFigures = false;
fprintf('\n--- Running Q6 Taylor-Green ---\n');
outQ6 = run_projection_taylor_green_demo(paramsQ6);

paramsQ9 = params;
paramsQ9.massFluxProjectionMode = 'relax_to_uniform_lowk';
paramsQ9.massFluxProjectionStrength = getf(params, 'massFluxProjectionStrength', 1.0);
paramsQ9.massFluxDensityRelaxationBeta = getf(params, 'massFluxDensityRelaxationBeta', 0.002);
paramsQ9.massFluxApplyAfterVelocityProjection = true;
paramsQ9.massFluxTargetFilter = 'lowpass_fft';
paramsQ9.massFluxLowKMaxIndex = getf(params, 'massFluxLowKMaxIndex', 2);
paramsQ9.massFluxFinalVelocityProjectionCleanup = getf(params, 'massFluxFinalVelocityProjectionCleanup', true);
paramsQ9.massFluxFinalVelocityProjectionStrength = getf(params, 'massFluxFinalVelocityProjectionStrength', 0.5);
paramsQ9.methodName = 'Q9_lowk_mass_flux';
fprintf('\n--- Running Q9 Taylor-Green ---\n');
outQ9 = run_projection_taylor_green_demo(paramsQ9);

summaryTable = build_summary_table(outQ6, outQ9);

cmp = struct();
cmp.params = params;
cmp.Q6 = outQ6;
cmp.Q9 = outQ9;
cmp.summaryTable = summaryTable;
cmp.ratios = table_to_ratio_struct(summaryTable);

fprintf('\n--- Taylor-Green Q6/Q9 summary ---\n');
disp(summaryTable);

if getf(params, 'makeComparisonFigures', false)
    make_comparison_figures(cmp);
end
end

function T = build_summary_table(outQ6, outQ9)
methods = ["Q6_projection"; "Q9_lowk_mass_flux"];
outs = {outQ6; outQ9};
rows = cell(2, 1);
for i = 1:2
    s = outs{i}.summary;
    p = outs{i}.params;
    rows{i} = {methods(i), p.nSteps, p.seed, outs{i}.actualLastStep, outs{i}.stoppedEarly, ...
        p.Nx, p.Ny, p.gamma, p.taylorGreenAmplitude, p.taylorGreenModeX, p.taylorGreenModeY, ...
        s.meanStdN, s.meanOutBand, s.meanLowKEnergy, s.timeAvgRelRms, ...
        s.meanRmsDivParticleAfter, s.meanMassFluxDivResidual, s.meanMassFluxDivParticleAfter, ...
        s.meanTGAmplitude, s.finalTGAmplitude, s.amplitudeRetention, ...
        s.meanTGModeEnergy, s.finalTGModeEnergy, s.modeEnergyRetention, ...
        s.meanTGCoherence, s.finalTGCoherence, ...
        s.meanEnstrophy, s.finalEnstrophy, ...
        s.meanHighKVelocityFraction, s.finalHighKVelocityFraction, ...
        s.meanKBTCell, s.finalKBTCell};
end
T = cell2table(vertcat(rows{:}), 'VariableNames', { ...
    'method','nSteps','seed','actualLastStep','stoppedEarly', ...
    'Nx','Ny','gamma','initialAmplitude','modeX','modeY', ...
    'meanStdN','meanOutBand','meanLowK','timeAvgRelRms', ...
    'meanDivUAfterParticles','meanMassFluxResidual','meanMassFluxAfterParticles', ...
    'meanTGAmplitude','finalTGAmplitude','amplitudeRetention', ...
    'meanTGModeEnergy','finalTGModeEnergy','modeEnergyRetention', ...
    'meanTGCoherence','finalTGCoherence', ...
    'meanEnstrophy','finalEnstrophy', ...
    'meanHighKVelocityFraction','finalHighKVelocityFraction', ...
    'meanKBTCell','finalKBTCell'});

% Append Q9/Q6 ratios where meaningful.  NaN is kept when the Q6 reference is
% absent or non-finite.
ratioNames = {'meanStdN','meanOutBand','meanLowK','timeAvgRelRms', ...
    'meanDivUAfterParticles','meanTGAmplitude','finalTGAmplitude','amplitudeRetention', ...
    'meanTGModeEnergy','finalTGModeEnergy','modeEnergyRetention', ...
    'meanTGCoherence','finalTGCoherence','meanEnstrophy','finalEnstrophy', ...
    'meanHighKVelocityFraction','finalHighKVelocityFraction','meanKBTCell'};
for k = 1:numel(ratioNames)
    name = ratioNames{k};
    ratio = safe_ratio(T{2, name}, T{1, name});
    varName = ['ratio_' name '_vs_Q6'];
    T.(varName) = nan(height(T), 1);
    T.(varName)(1) = 1.0;
    T.(varName)(2) = ratio;
end
end

function r = table_to_ratio_struct(T)
r = struct();
vars = T.Properties.VariableNames;
for k = 1:numel(vars)
    name = vars{k};
    if startsWith(name, 'ratio_')
        r.(name) = T{2, name};
    end
end
end

function make_comparison_figures(cmp)
figure('Name','Taylor-Green Q6/Q9 comparison','Color','w');
subplot(2,2,1);
plot(cmp.Q6.sampleTimes, cmp.Q6.tgAmplitudeSeries, 'LineWidth', 1.2); hold on;
plot(cmp.Q9.sampleTimes, cmp.Q9.tgAmplitudeSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('TG amplitude'); legend('Q6','Q9'); grid on;
subplot(2,2,2);
plot(cmp.Q6.sampleTimes, cmp.Q6.tgCoherenceSeries, 'LineWidth', 1.2); hold on;
plot(cmp.Q9.sampleTimes, cmp.Q9.tgCoherenceSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('coherence'); legend('Q6','Q9'); grid on;
subplot(2,2,3);
plot(cmp.Q6.sampleTimes, cmp.Q6.enstrophySeries, 'LineWidth', 1.2); hold on;
plot(cmp.Q9.sampleTimes, cmp.Q9.enstrophySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('enstrophy'); legend('Q6','Q9'); grid on;
subplot(2,2,4);
plot(cmp.Q6.sampleTimes, cmp.Q6.highKVelocityFractionSeries, 'LineWidth', 1.2); hold on;
plot(cmp.Q9.sampleTimes, cmp.Q9.highKVelocityFractionSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('high-k velocity fraction'); legend('Q6','Q9'); grid on;
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
params = set_default(params, 'dt', 2.0e-3);
params = set_default(params, 'kBT', 0.05);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'taylorGreenAmplitude', 0.20);
params = set_default(params, 'taylorGreenModeX', 1);
params = set_default(params, 'taylorGreenModeY', 1);
params = set_default(params, 'taylorGreenThermalNoise', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'massFluxProjectionStrength', 1.0);
params = set_default(params, 'massFluxDensityRelaxationBeta', 0.002);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', true);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', true);
params = set_default(params, 'massFluxFinalVelocityProjectionStrength', 0.5);
params = set_default(params, 'massFluxTargetFilter', 'lowpass_fft');
params = set_default(params, 'massFluxLowKMaxIndex', 2);
params = set_default(params, 'thermostatAfterProjection', true);
params = set_default(params, 'thermostatTargetKBT', params.kBT);
params = set_default(params, 'thermostatStrength', 1.0);
params = set_default(params, 'densityBandFraction', 0.20);
params = set_default(params, 'lowKMaxIndex', 2);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'storeVorticityMaps', false);
params = set_default(params, 'computeFullDiagnosticsEveryStep', false);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'makeComparisonFigures', false);
params = set_default(params, 'seed', 11);
params = set_default(params, 'maxWallClockSeconds', Inf);
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

function r = safe_ratio(a, b)
if isfinite(a) && isfinite(b) && abs(b) > eps
    r = a / b;
else
    r = NaN;
end
end
