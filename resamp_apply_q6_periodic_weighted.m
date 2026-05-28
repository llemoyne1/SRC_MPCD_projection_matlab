function [stateOut, diag] = resamp_apply_q6_periodic_weighted(stateClassic, params)
%RESAMP_APPLY_Q6_PERIODIC_WEIGHTED Apply periodic Q6 velocity projection to weighted active particles.

validate_state(stateClassic);
projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', ...
    get_param(params, 'thermostatAfterStep', false)));
computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));

activeMask = resamp_active_mask(stateClassic);
Gbefore = resamp_deposit_weighted_to_grid(stateClassic.x, stateClassic.v, stateClassic.m, params, ...
    'periodicX', true, 'periodicY', true, 'minMass', eps, 'activeMask', activeMask);
proj = projection_project_grid_periodic_fft(Gbefore.Ux, Gbefore.Uy, params);

stateOut = stateClassic;
dvApplied = zeros(size(stateClassic.v));
if projectionEnable && projectionStrength ~= 0 && any(activeMask)
    dvA = projection_interpolate_grid_delta_to_particles(stateClassic.x(activeMask,:), proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    dvApplied(activeMask,:) = projectionStrength * dvA;
    [stateOut.v, momCorr] = resamp_apply_global_momentum_correction_weighted( ...
        stateClassic.v, dvApplied, stateClassic.m, params, 'stage', 'weighted_q6', 'activeMask', activeMask);
else
    momCorr = resamp_apply_global_momentum_correction_weighted( ...
        stateClassic.v, dvApplied, stateClassic.m, params, 'stage', 'weighted_q6_disabled', 'activeMask', activeMask);
end

if thermostatAfterProjection
    [stateOut.v, thermostatInfo] = resamp_apply_cell_thermostat_weighted(stateOut.x, stateOut.v, stateOut.m, params, ...
        'periodicX', true, 'periodicY', true, 'activeMask', resamp_active_mask(stateOut));
else
    thermostatInfo = empty_thermostat_info();
end
stateOut.Nactive = nnz(resamp_active_mask(stateOut));
stateOut.Ncapacity = size(stateOut.x, 1);

Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', true, 'periodicY', true, 'minMass', eps, 'activeMask', resamp_active_mask(stateOut));
projAfter = projection_project_grid_periodic_fft(Gafter.Ux, Gafter.Uy, params);
mdAfter = resamp_population_mass_diagnostics(stateOut, params, 'periodicX', true, 'periodicY', true);
pool = resamp_particle_pool_info(stateOut);

diag = struct();
diag.kind = 'weighted_q6';
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.thermostatAfterProjection = thermostatAfterProjection;
diag.thermostat = thermostatInfo;
diag.NpActive = pool.Nactive;
diag.Ncapacity = pool.Ncapacity;
diag.Nfree = pool.Nfree;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivAfterGrid = proj.rmsDivAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivAfterGrid = proj.maxAbsDivAfter;
diag.rmsDivParticleAfter = projAfter.rmsDivBefore;
diag.divReductionGrid = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
diag.divReductionParticle = diag.rmsDivParticleAfter / max(proj.rmsDivBefore, eps);
if any(activeMask)
    diag.dvAppliedRms = sqrt(mean(sum(dvApplied(activeMask,:).^2, 2), 'omitnan'));
else
    diag.dvAppliedRms = 0;
end
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
    error('state.m must have one entry per particle slot.');
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
