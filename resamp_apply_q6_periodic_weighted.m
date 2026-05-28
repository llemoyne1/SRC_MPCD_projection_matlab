function [stateOut, diag] = resamp_apply_q6_periodic_weighted(stateClassic, params)
%RESAMP_APPLY_Q6_PERIODIC_WEIGHTED Apply periodic Q6 velocity projection to weighted particles.
%
% The projected grid velocity is reconstructed as U = sum(m v)/sum(m), but
% the projection correction remains a velocity correction.  The particle
% kick is then corrected with a weighted global momentum correction so that
% the projection stage does not introduce a net weighted momentum drift.

validate_state(stateClassic);
projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', ...
    get_param(params, 'thermostatAfterStep', false)));
computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));

Gbefore = resamp_deposit_weighted_to_grid(stateClassic.x, stateClassic.v, stateClassic.m, params, ...
    'periodicX', true, 'periodicY', true, 'minMass', eps);
proj = projection_project_grid_periodic_fft(Gbefore.Ux, Gbefore.Uy, params);

stateOut = stateClassic;
if projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    dvApplied = projectionStrength * dv;
    [stateOut.v, momCorr] = resamp_apply_global_momentum_correction_weighted( ...
        stateClassic.v, dvApplied, stateClassic.m, params, 'stage', 'weighted_q6');
else
    dvApplied = zeros(size(stateClassic.v));
    momCorr = resamp_apply_global_momentum_correction_weighted( ...
        stateClassic.v, dvApplied, stateClassic.m, params, 'stage', 'weighted_q6_disabled');
end

if thermostatAfterProjection
    [stateOut.v, thermostatInfo] = resamp_apply_cell_thermostat_weighted(stateOut.x, stateOut.v, stateOut.m, params, ...
        'periodicX', true, 'periodicY', true);
else
    thermostatInfo = empty_thermostat_info();
end

Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', true, 'periodicY', true, 'minMass', eps);
projAfter = projection_project_grid_periodic_fft(Gafter.Ux, Gafter.Uy, params);
mdAfter = resamp_population_mass_diagnostics(stateOut, params, 'periodicX', true, 'periodicY', true);

diag = struct();
diag.kind = 'weighted_q6';
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.thermostatAfterProjection = thermostatAfterProjection;
diag.thermostat = thermostatInfo;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivAfterGrid = proj.rmsDivAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivAfterGrid = proj.maxAbsDivAfter;
diag.rmsDivParticleAfter = projAfter.rmsDivBefore;
diag.divReductionGrid = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
diag.divReductionParticle = diag.rmsDivParticleAfter / max(proj.rmsDivBefore, eps);
diag.dvAppliedRms = sqrt(mean(sum(dvApplied.^2, 2), 'omitnan'));
diag.momentumCorrection = momCorr;
diag.Gbefore = [];
diag.Gafter = [];
diag.proj = [];
diag.projAfterParticles = [];
if computeDiagnostics
    diag.Gbefore = Gbefore;
    diag.Gafter = Gafter;
    diag.proj = proj;
    diag.projAfterParticles = projAfter;
end
diag.NStdAfterProjection = mdAfter.NStd;
diag.MStdAfterProjection = mdAfter.MStd;
diag.MRelRmsAfterProjection = mdAfter.MRelRms;
diag.totalMassAfterProjection = mdAfter.totalMass;
diag.totalMomentumAfterProjection = mdAfter.totalMomentum;
diag.massDiagnosticsAfterProjection = mdAfter;
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must be a struct with fields x, v and m.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
if numel(state.m) ~= size(state.x,1)
    error('state.m must have one entry per particle.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end

function info = empty_thermostat_info()
info = struct();
info.enabled = false;
info.kBTTarget = NaN;
info.kBTBefore = NaN;
info.kBTAfter = NaN;
info.meanScale = NaN;
info.maxScale = NaN;
end
