function cmp = run_compare_projection_cylinder_q6_q9(params)
%RUN_COMPARE_PROJECTION_CYLINDER_Q6_Q9 Compare Q6 and Q9 on cylinder wake.
%
%   cmp = run_compare_projection_cylinder_q6_q9(params)
%
% Runs two simulations with identical seed, geometry and forcing:
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

fprintf('\n=== Cylinder comparison: Q6 projection-only vs Q9 low-k mass flux ===\n');
fprintf('\n--- Running Q6 cylinder ---\n');
outQ6 = run_projection_cylinder_demo(paramsQ6);

fprintf('\n--- Running Q9 cylinder ---\n');
outQ9 = run_projection_cylinder_demo(paramsQ9);

summary = build_summary_table(outQ6, outQ9);

cmp = struct();
cmp.paramsQ6 = paramsQ6;
cmp.paramsQ9 = paramsQ9;
cmp.outQ6 = outQ6;
cmp.outQ9 = outQ9;
cmp.summary = summary;

disp(summary);
end

function summary = build_summary_table(outQ6, outQ9)
s6 = outQ6.summary;
s9 = outQ9.summary;
method = string({'Q6_projection'; 'Q9_lowk_mass_flux'});
nSteps = [outQ6.params.nSteps; outQ9.params.nSteps];
seed = [outQ6.params.seed; outQ9.params.seed];
meanStdN = [s6.meanStdN; s9.meanStdN];
meanOutBand = [s6.meanOutBand; s9.meanOutBand];
meanLowK = [s6.meanLowKEnergy; s9.meanLowKEnergy];
timeAvgRelRms = [s6.timeAvgRelRms; s9.timeAvgRelRms];
meanRhoTransport = [s6.meanDensityTransportProjectedRms; s9.meanDensityTransportProjectedRms];
meanDivUAfter = [s6.meanRmsDivParticleAfter; s9.meanRmsDivParticleAfter];
meanMassFluxResidual = [s6.meanMassFluxDivResidual; s9.meanMassFluxDivResidual];
meanMassFluxParticleAfter = [s6.meanMassFluxDivParticleAfter; s9.meanMassFluxDivParticleAfter];
meanUx = [s6.meanMeanUx; s9.meanMeanUx];
finalUx = [s6.finalMeanUx; s9.finalMeanUx];
meanEnstrophy = [s6.meanEnstrophy; s9.meanEnstrophy];
finalEnstrophy = [s6.finalEnstrophy; s9.finalEnstrophy];
wakeOmegaRms = [s6.wakeOmegaRms; s9.wakeOmegaRms];
wakeUyRms = [s6.wakeUyRms; s9.wakeUyRms];
wakeOmegaUpperRms = [s6.wakeOmegaUpperRms; s9.wakeOmegaUpperRms];
wakeUyUpperRms = [s6.wakeUyUpperRms; s9.wakeUyUpperRms];
wakeOmegaAntiSymRms = [s6.wakeOmegaAntiSymRms; s9.wakeOmegaAntiSymRms];
wakeUyAntiSymRms = [s6.wakeUyAntiSymRms; s9.wakeUyAntiSymRms];
wakeOmegaSignChanges = [s6.wakeOmegaSignChanges; s9.wakeOmegaSignChanges];
wakeUyUpperSignChanges = [s6.wakeUyUpperSignChanges; s9.wakeUyUpperSignChanges];
wakeOmegaAntiSymSignChanges = [s6.wakeOmegaAntiSymSignChanges; s9.wakeOmegaAntiSymSignChanges];
sheddingFrequencyWakeUyUpper = [s6.sheddingFrequencyWakeUyUpper; s9.sheddingFrequencyWakeUyUpper];
strouhalWakeUyUpper = [s6.strouhalWakeUyUpper; s9.strouhalWakeUyUpper];
sheddingSNRWakeUyUpper = [s6.sheddingSNRWakeUyUpper; s9.sheddingSNRWakeUyUpper];
sheddingFrequencyWakeOmegaAntiSym = [s6.sheddingFrequencyWakeOmegaAntiSym; s9.sheddingFrequencyWakeOmegaAntiSym];
strouhalWakeOmegaAntiSym = [s6.strouhalWakeOmegaAntiSym; s9.strouhalWakeOmegaAntiSym];
sheddingSNROmegaAntiSym = [s6.sheddingSNROmegaAntiSym; s9.sheddingSNROmegaAntiSym];
reEffMeanUx = [s6.reEffMeanUx; s9.reEffMeanUx];
reEffFinalUx = [s6.reEffFinalUx; s9.reEffFinalUx];
meanCylinderHits = [s6.meanCylinderHits; s9.meanCylinderHits];
meanKBTCell = [s6.kBTCellMean; s9.kBTCellMean];

ratioStd_vs_Q6 = ratio_to_first(meanStdN);
ratioOutBand_vs_Q6 = ratio_to_first(meanOutBand);
ratioLowK_vs_Q6 = ratio_to_first(meanLowK);
ratioRhoTransport_vs_Q6 = ratio_to_first(meanRhoTransport);
ratioDivUAfter_vs_Q6 = ratio_to_first(meanDivUAfter);
ratioEnstrophy_vs_Q6 = ratio_to_first(meanEnstrophy);
ratioWakeOmegaRms_vs_Q6 = ratio_to_first(wakeOmegaRms);
ratioWakeUyRms_vs_Q6 = ratio_to_first(wakeUyRms);
ratioWakeOmegaUpperRms_vs_Q6 = ratio_to_first(wakeOmegaUpperRms);
ratioWakeUyUpperRms_vs_Q6 = ratio_to_first(wakeUyUpperRms);
ratioWakeOmegaAntiSymRms_vs_Q6 = ratio_to_first(wakeOmegaAntiSymRms);
ratioWakeUyAntiSymRms_vs_Q6 = ratio_to_first(wakeUyAntiSymRms);

summary = table(method, nSteps, seed, meanStdN, meanOutBand, meanLowK, ...
    timeAvgRelRms, meanRhoTransport, meanDivUAfter, meanMassFluxResidual, ...
    meanMassFluxParticleAfter, meanUx, finalUx, meanEnstrophy, finalEnstrophy, ...
    wakeOmegaRms, wakeUyRms, wakeOmegaUpperRms, wakeUyUpperRms, ...
    wakeOmegaAntiSymRms, wakeUyAntiSymRms, wakeOmegaSignChanges, wakeUyUpperSignChanges, ...
    wakeOmegaAntiSymSignChanges, sheddingFrequencyWakeUyUpper, strouhalWakeUyUpper, ...
    sheddingSNRWakeUyUpper, sheddingFrequencyWakeOmegaAntiSym, strouhalWakeOmegaAntiSym, ...
    sheddingSNROmegaAntiSym, reEffMeanUx, reEffFinalUx, meanCylinderHits, meanKBTCell, ...
    ratioStd_vs_Q6, ratioOutBand_vs_Q6, ratioLowK_vs_Q6, ratioRhoTransport_vs_Q6, ...
    ratioDivUAfter_vs_Q6, ratioEnstrophy_vs_Q6, ratioWakeOmegaRms_vs_Q6, ratioWakeUyRms_vs_Q6, ...
    ratioWakeOmegaUpperRms_vs_Q6, ratioWakeUyUpperRms_vs_Q6, ...
    ratioWakeOmegaAntiSymRms_vs_Q6, ratioWakeUyAntiSymRms_vs_Q6);
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
