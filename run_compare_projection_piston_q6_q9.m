function cmp = run_compare_projection_piston_q6_q9(params)
%RUN_COMPARE_PROJECTION_PISTON_Q6_Q9 Compare Q6 and Q9 on moving piston.
%
%   cmp = run_compare_projection_piston_q6_q9(params)
%
% Runs two piston simulations with identical seed and kinematics:
%   Q6 : velocity projection only
%   Q9 : velocity projection + low-k mass-flux relaxation

if nargin < 1 || isempty(params)
    params = struct();
end

params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'storeDensityMaps', true);
params = set_default(params, 'makeFigures', false);
params = set_default(params, 'massFluxLowKMaxIndex', get_param(params, 'lowKMaxIndex', 2));

fprintf('\n=== Piston comparison: Q6 projection-only vs Q9 low-k mass flux ===\n');

paramsQ6 = params;
paramsQ6.massFluxProjectionMode = 'off';
paramsQ6.massFluxProjectionStrength = paramsQ6.projectionStrength;
paramsQ6.massFluxDensityRelaxationBeta = 0.0;
paramsQ6.massFluxApplyAfterVelocityProjection = false;
paramsQ6.massFluxFinalVelocityProjectionCleanup = false;
paramsQ6.massFluxFinalVelocityProjectionStrength = 0.0;
paramsQ6.massFluxTargetFilter = 'none';

paramsQ9 = params;
paramsQ9.massFluxProjectionMode = 'relax_to_uniform_lowk';
paramsQ9.massFluxProjectionStrength = 1.0;
paramsQ9.massFluxDensityRelaxationBeta = get_param(paramsQ9, 'massFluxDensityRelaxationBeta', 0.002);
paramsQ9.massFluxApplyAfterVelocityProjection = true;
paramsQ9.massFluxFinalVelocityProjectionCleanup = logical(get_param(paramsQ9, 'massFluxFinalVelocityProjectionCleanup', false));
paramsQ9.massFluxFinalVelocityProjectionStrength = get_param(paramsQ9, 'massFluxFinalVelocityProjectionStrength', 1.0);
paramsQ9.massFluxTargetFilter = 'lowpass_fft';
paramsQ9.massFluxLowKMaxIndex = get_param(paramsQ9, 'massFluxLowKMaxIndex', 2);

fprintf('\n--- Running Q6 piston ---\n');
outQ6 = run_projection_piston_demo(paramsQ6);

fprintf('\n--- Running Q9 piston ---\n');
outQ9 = run_projection_piston_demo(paramsQ9);

summary = build_summary_table(outQ6, outQ9);

disp(summary);

cmp = struct();
cmp.paramsBase = params;
cmp.outQ6 = outQ6;
cmp.outQ9 = outQ9;
cmp.summary = summary;
cmp.description = 'Q6 velocity projection versus Q9 low-k mass-flux relaxation on moving piston';
end

function T = build_summary_table(outQ6, outQ9)
S6 = outQ6.summary;
S9 = outQ9.summary;
T = table();
T.method = ["Q6_projection"; "Q9_lowk_mass_flux"];
T.nSteps = [outQ6.params.nSteps; outQ9.params.nSteps];
T.seed = [outQ6.params.seed; outQ9.params.seed];
T.finalCompression = [S6.finalCompression; S9.finalCompression];
T.meanStdN = [S6.meanStdN; S9.meanStdN];
T.meanOutBand = [S6.meanOutBandFraction; S9.meanOutBandFraction];
T.meanLowK = [S6.meanLowKEnergy; S9.meanLowKEnergy];
T.timeAvgRelRms = [S6.timeAvgRelRms; S9.timeAvgRelRms];
T.meanDivUAfter = [S6.meanRmsDivParticleAfter; S9.meanRmsDivParticleAfter];
T.meanDivReductionParticle = [S6.meanDivReductionParticle; S9.meanDivReductionParticle];
T.meanRhoTransport = [S6.meanDensityTransportProjectedRms; S9.meanDensityTransportProjectedRms];
T.meanMassFluxResidual = [S6.meanMassFluxDivResidual; S9.meanMassFluxDivResidual];
T.meanMassFluxParticleAfter = [S6.meanMassFluxDivParticleAfter; S9.meanMassFluxDivParticleAfter];
T.meanPkin = [S6.meanPkinMean; S9.meanPkinMean];
T.finalPkin = [S6.finalPkinMean; S9.finalPkinMean];
T.meanKBTCell = [S6.meanKBTCell; S9.meanKBTCell];

T.ratioStd_vs_Q6 = [1; safe_ratio(S9.meanStdN, S6.meanStdN)];
T.ratioOutBand_vs_Q6 = [1; safe_ratio(S9.meanOutBandFraction, S6.meanOutBandFraction)];
T.ratioLowK_vs_Q6 = [1; safe_ratio(S9.meanLowKEnergy, S6.meanLowKEnergy)];
T.ratioDivUAfter_vs_Q6 = [1; safe_ratio(S9.meanRmsDivParticleAfter, S6.meanRmsDivParticleAfter)];
T.ratioRhoTransport_vs_Q6 = [1; safe_ratio(S9.meanDensityTransportProjectedRms, S6.meanDensityTransportProjectedRms)];
T.ratioMassFluxResidual_vs_Q6 = [1; safe_ratio(S9.meanMassFluxDivResidual, S6.meanMassFluxDivResidual)];
T.ratioMassFluxParticleAfter_vs_Q6 = [1; safe_ratio(S9.meanMassFluxDivParticleAfter, S6.meanMassFluxDivParticleAfter)];
end

function r = safe_ratio(num, den)
if isempty(num) || isempty(den) || ~isnumeric(num) || ~isnumeric(den) || ...
        ~isscalar(num) || ~isscalar(den) || ~isfinite(num) || ~isfinite(den) || abs(den) <= eps
    r = NaN;
else
    r = num / den;
end
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
