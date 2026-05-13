function cmp = run_compare_density_projection_lowk_mass_flux_poiseuille(params)
%RUN_COMPARE_DENSITY_PROJECTION_LOWK_MASS_FLUX_POISEUILLE Q9 hybrid comparison.
%
%   cmp = run_compare_density_projection_lowk_mass_flux_poiseuille(params)
%
% Q9 is a hybrid low-k density correction:
%   C : SRC classic
%   P : Q6 velocity projection, div(u)=0
%   H : Q6 velocity projection followed by low-k mass-flux relaxation
%
% The low-k target is beta/dt * lowpass(N-gamma).  It aims to keep the
% spectral benefit of Q8 on coherent density modes while avoiding corrections
% driven by cell-scale/high-k particle noise.

if nargin < 1 || isempty(params)
    params = struct();
end

params = set_default(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk');
params = set_default(params, 'massFluxProjectionStrength', 1.0);
params = set_default(params, 'massFluxDensityRelaxationBeta', 0.002);
params = set_default(params, 'massFluxApplyAfterVelocityProjection', true);
params = set_default(params, 'massFluxTargetFilter', 'lowpass_fft');
params = set_default(params, 'massFluxLowKMaxIndex', get_param(params, 'lowKMaxIndex', 2));

fprintf('\n=== Q9 low-k mass-flux hybrid comparison ===\n');
fprintf('mode=%s, beta=%.6g, lowK=%d, afterVelocity=%d\n', ...
    char(params.massFluxProjectionMode), params.massFluxDensityRelaxationBeta, ...
    params.massFluxLowKMaxIndex, params.massFluxApplyAfterVelocityProjection);

cmp = run_compare_density_projection_mass_flux_poiseuille(params);
cmp.q9Description = 'velocity projection + low-k mass-flux density relaxation';
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
