function cmp = run_compare_projection_step_q6_q9(params)
%RUN_COMPARE_PROJECTION_STEP_Q6_Q9 Compare Q6 and Q9 on grid-aligned step flow.
%
%   cmp = run_compare_projection_step_q6_q9(params)
%
% Runs two simulations with identical seed and geometry:
%   Q6 : velocity projection only
%   Q9 : velocity projection + low-k mass-flux relaxation

if nargin < 1 || isempty(params)
    params = struct();
end

params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'storeVorticityMaps', true);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'massFluxFinalVelocityProjectionCleanup', true);
params = set_default(params, 'massFluxFinalVelocityProjectionStrength', 0.5);

paramsQ6 = params;
paramsQ6.massFluxProjectionMode = 'off';
paramsQ6.massFluxProjectionStrength = 1.0;
paramsQ6.massFluxDensityRelaxationBeta = 0.0;
paramsQ6.massFluxApplyAfterVelocityProjection = false;
paramsQ6.massFluxFinalVelocityProjectionCleanup = false;
paramsQ6.massFluxFinalVelocityProjectionStrength = 0.0;
paramsQ6.makeFigures = false;
paramsQ6.visualTitleSuffix = 'Q6 projection';
paramsQ6.visualFigureId = get_field(params, 'visualFigureId', 101);

paramsQ9 = params;
paramsQ9.projectionStrength = get_field(params, 'projectionStrength', 1.0);
paramsQ9.massFluxProjectionMode = get_field(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk');
paramsQ9.massFluxProjectionStrength = get_field(params, 'massFluxProjectionStrength', 1.0);
paramsQ9.massFluxDensityRelaxationBeta = get_field(params, 'massFluxDensityRelaxationBeta', 0.002);
paramsQ9.massFluxApplyAfterVelocityProjection = get_field(params, 'massFluxApplyAfterVelocityProjection', true);
paramsQ9.massFluxFinalVelocityProjectionCleanup = get_field(params, 'massFluxFinalVelocityProjectionCleanup', true);
paramsQ9.massFluxFinalVelocityProjectionStrength = get_field(params, 'massFluxFinalVelocityProjectionStrength', 0.5);
paramsQ9.massFluxTargetFilter = get_field(params, 'massFluxTargetFilter', 'lowpass_fft');
paramsQ9.massFluxLowKMaxIndex = get_field(params, 'massFluxLowKMaxIndex', 2);
paramsQ9.makeFigures = false;
paramsQ9.visualTitleSuffix = 'Q9 low-k mass-flux';
paramsQ9.visualFigureId = get_field(params, 'visualFigureId', 101);

fprintf('\n=== Step-channel comparison: Q6 projection-only vs Q9 low-k mass flux ===\n');
fprintf('\n--- Running Q6 step-channel ---\n');
outQ6 = run_projection_step_channel_demo(paramsQ6);

fprintf('\n--- Running Q9 step-channel ---\n');
outQ9 = run_projection_step_channel_demo(paramsQ9);

summaryTable = build_summary_table(outQ6, outQ9);

cmp = struct();
cmp.paramsQ6 = paramsQ6;
cmp.paramsQ9 = paramsQ9;
cmp.outQ6 = outQ6;
cmp.outQ9 = outQ9;
cmp.summaryTable = summaryTable;

disp(summaryTable);

if logical(get_field(params, 'makeComparisonFigures', false))
    make_comparison_figures(outQ6, outQ9);
end
end

function summary = build_summary_table(outQ6, outQ9)
s6 = outQ6.summary;
s9 = outQ9.summary;
method = string({'Q6_projection'; 'Q9_lowk_mass_flux'});
nSteps = [outQ6.params.nSteps; outQ9.params.nSteps];
seed = [outQ6.params.seed; outQ9.params.seed];
actualLastStep = [s6.actualLastStep; s9.actualLastStep];
stoppedEarly = [s6.stoppedEarly; s9.stoppedEarly];
meanStdN = [s6.meanStdN; s9.meanStdN];
meanOutBand = [s6.meanOutBand; s9.meanOutBand];
meanLowK = [s6.meanLowKEnergy; s9.meanLowKEnergy];
timeAvgRelRms = [s6.timeAvgRelRms; s9.timeAvgRelRms];
meanRhoTransport = [s6.meanDensityTransportProjectedRms; s9.meanDensityTransportProjectedRms];
meanDivUAfter = [s6.meanRmsDivParticleAfter; s9.meanRmsDivParticleAfter];
meanMassFluxResidual = [s6.meanMassFluxDivResidual; s9.meanMassFluxDivResidual];
meanMassFluxAfter = [s6.meanMassFluxDivParticleAfter; s9.meanMassFluxDivParticleAfter];
meanUx = [s6.meanMeanUx; s9.meanMeanUx];
finalUx = [s6.finalMeanUx; s9.finalMeanUx];
meanEnstrophy = [s6.meanEnstrophy; s9.meanEnstrophy];
finalEnstrophy = [s6.finalEnstrophy; s9.finalEnstrophy];
meanShearOmegaRms = [s6.meanShearOmegaRms; s9.meanShearOmegaRms];
finalShearOmegaRms = [s6.finalShearOmegaRms; s9.finalShearOmegaRms];
meanRecirculationLength = [s6.meanRecirculationLength; s9.meanRecirculationLength];
finalRecirculationLength = [s6.finalRecirculationLength; s9.finalRecirculationLength];
meanRecirculationArea = [s6.meanRecirculationArea; s9.meanRecirculationArea];
finalRecirculationArea = [s6.finalRecirculationArea; s9.finalRecirculationArea];
meanRecirculationMinUx = [s6.meanRecirculationMinUx; s9.meanRecirculationMinUx];
probeOmegaRms = [s6.probeOmegaRms; s9.probeOmegaRms];
probeUyRms = [s6.probeUyRms; s9.probeUyRms];
meanStepHits = [s6.meanStepHits; s9.meanStepHits];
meanKBTCell = [s6.meanKBTCell; s9.meanKBTCell];
finalKBTCell = [s6.finalKBTCell; s9.finalKBTCell];

ratioStd_vs_Q6 = ratio_to_first(meanStdN);
ratioOutBand_vs_Q6 = ratio_to_first(meanOutBand);
ratioLowK_vs_Q6 = ratio_to_first(meanLowK);
ratioRhoTransport_vs_Q6 = ratio_to_first(meanRhoTransport);
ratioDivUAfter_vs_Q6 = ratio_to_first(meanDivUAfter);
ratioEnstrophy_vs_Q6 = ratio_to_first(meanEnstrophy);
ratioShearOmegaRms_vs_Q6 = ratio_to_first(meanShearOmegaRms);
ratioRecircLength_vs_Q6 = ratio_to_first(meanRecirculationLength);
ratioRecircArea_vs_Q6 = ratio_to_first(meanRecirculationArea);
ratioProbeOmegaRms_vs_Q6 = ratio_to_first(probeOmegaRms);
ratioProbeUyRms_vs_Q6 = ratio_to_first(probeUyRms);

summary = table(method, nSteps, seed, actualLastStep, stoppedEarly, ...
    meanStdN, meanOutBand, meanLowK, timeAvgRelRms, meanRhoTransport, ...
    meanDivUAfter, meanMassFluxResidual, meanMassFluxAfter, meanUx, finalUx, ...
    meanEnstrophy, finalEnstrophy, meanShearOmegaRms, finalShearOmegaRms, ...
    meanRecirculationLength, finalRecirculationLength, meanRecirculationArea, ...
    finalRecirculationArea, meanRecirculationMinUx, probeOmegaRms, probeUyRms, ...
    meanStepHits, meanKBTCell, finalKBTCell, ratioStd_vs_Q6, ratioOutBand_vs_Q6, ...
    ratioLowK_vs_Q6, ratioRhoTransport_vs_Q6, ratioDivUAfter_vs_Q6, ...
    ratioEnstrophy_vs_Q6, ratioShearOmegaRms_vs_Q6, ratioRecircLength_vs_Q6, ...
    ratioRecircArea_vs_Q6, ratioProbeOmegaRms_vs_Q6, ratioProbeUyRms_vs_Q6);
end

function make_comparison_figures(outQ6, outQ9)
figure('Name','Step-channel Q6/Q9 comparison','Color','w');
subplot(2,2,1);
plot(outQ6.sampleTimes, outQ6.meanEnstrophySeries, 'LineWidth', 1.2); hold on;
plot(outQ9.sampleTimes, outQ9.meanEnstrophySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('mean enstrophy'); grid on; legend('Q6','Q9');
subplot(2,2,2);
plot(outQ6.sampleTimes, outQ6.shearOmegaRmsSeries, 'LineWidth', 1.2); hold on;
plot(outQ9.sampleTimes, outQ9.shearOmegaRmsSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('shear \omega RMS'); grid on; legend('Q6','Q9');
subplot(2,2,3);
plot(outQ6.sampleTimes, outQ6.recirculationLengthSeries, 'LineWidth', 1.2); hold on;
plot(outQ9.sampleTimes, outQ9.recirculationLengthSeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('recirculation length'); grid on; legend('Q6','Q9');
subplot(2,2,4);
plot(outQ6.sampleTimes, outQ6.probeUySeries, 'LineWidth', 1.2); hold on;
plot(outQ9.sampleTimes, outQ9.probeUySeries, 'LineWidth', 1.2);
xlabel('t'); ylabel('probe Uy'); grid on; legend('Q6','Q9');
end

function r = ratio_to_first(v)
r = nan(size(v));
if isempty(v) || ~isfinite(v(1)) || abs(v(1)) < eps
    return;
end
r(1) = 1;
for i = 2:numel(v)
    if isfinite(v(i))
        r(i) = v(i) / v(1);
    end
end
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function val = get_field(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end
