function [stateOut, diag] = mpcd_apply_q9_projection_periodic(stateClassic, params)
%MPCD_APPLY_Q9_PROJECTION_PERIODIC Apply Q6/Q9 projection in a periodic box.
%
%   [stateOut, diag] = mpcd_apply_q9_projection_periodic(stateClassic, params)
%
% Periodic counterpart of mpcd_apply_q9_projection_channel.  It is intended
% for structure-preservation tests such as Taylor-Green and shear layers.

projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
massFluxProjectionMode = char(string(get_param(params, 'massFluxProjectionMode', 'off')));
massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', projectionStrength);
massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 0.0);
massFluxApplyAfterVelocityProjection = logical(get_param(params, 'massFluxApplyAfterVelocityProjection', false));
massFluxFinalVelocityProjectionCleanup = logical(get_param(params, 'massFluxFinalVelocityProjectionCleanup', false));
massFluxFinalVelocityProjectionStrength = get_param(params, 'massFluxFinalVelocityProjectionStrength', 1.0);
useMassFluxProjection = ~strcmpi(strrep(massFluxProjectionMode, '-', '_'), 'off');
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', false));
computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));

validate_state(stateClassic);

Gbefore = projection_deposit_particles_to_grid(stateClassic.x, stateClassic.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
proj = projection_project_grid_periodic_fft(Gbefore.Ux, Gbefore.Uy, params);
massFluxProj = empty_mass_flux_projection();
projectionKind = 'velocity';
appliedProjectionStrength = projectionStrength;
stateOut = stateClassic;
finalVelocityCleanup = empty_velocity_cleanup();

if computeDiagnostics
    popAfterClassic = projection_population_diagnostics(stateClassic.x, params, ...
        'periodicX', true, 'periodicY', true);
    thermalBeforeProjection = projection_thermal_diagnostics(stateClassic.x, stateClassic.v, params, ...
        'periodicX', true, 'periodicY', true, 'minCount', 1);
end

if useMassFluxProjection
    appliedProjectionStrength = massFluxProjectionStrength;
    if massFluxApplyAfterVelocityProjection
        projectionKind = 'velocity_plus_mass_flux';
        stateForMassFlux = stateClassic;
        if projectionEnable && projectionStrength ~= 0
            dvVelocity = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
                'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
            stateForMassFlux.v = stateClassic.v + projectionStrength * dvVelocity;
        end
        GforMassFlux = projection_deposit_particles_to_grid(stateForMassFlux.x, stateForMassFlux.v, params, ...
            'periodicX', true, 'periodicY', true, 'minCount', 1);
        massFluxProj = projection_project_mass_flux_periodic_fft(GforMassFlux.N, GforMassFlux.Ux, GforMassFlux.Uy, params, ...
            'mode', massFluxProjectionMode, ...
            'relaxationBeta', massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dvMassFlux = projection_interpolate_grid_delta_to_particles(stateForMassFlux.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
            stateOut.v = stateForMassFlux.v + massFluxProjectionStrength * dvMassFlux;
        else
            stateOut = stateForMassFlux;
        end
    else
        projectionKind = 'mass_flux';
        massFluxProj = projection_project_mass_flux_periodic_fft(Gbefore.N, Gbefore.Ux, Gbefore.Uy, params, ...
            'mode', massFluxProjectionMode, ...
            'relaxationBeta', massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
            stateOut.v = stateClassic.v + massFluxProjectionStrength * dv;
        end
    end
elseif projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    stateOut.v = stateClassic.v + projectionStrength * dv;
end

if useMassFluxProjection && massFluxFinalVelocityProjectionCleanup && massFluxFinalVelocityProjectionStrength ~= 0
    GcleanupBefore = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', true, 'minCount', 1);
    projCleanup = projection_project_grid_periodic_fft(GcleanupBefore.Ux, GcleanupBefore.Uy, params);
    dvCleanup = projection_interpolate_grid_delta_to_particles(stateOut.x, projCleanup.dUx, projCleanup.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    stateOut.v = stateOut.v + massFluxFinalVelocityProjectionStrength * dvCleanup;
    finalVelocityCleanup = struct();
    finalVelocityCleanup.enabled = true;
    finalVelocityCleanup.strength = massFluxFinalVelocityProjectionStrength;
    finalVelocityCleanup.rmsDivBefore = projCleanup.rmsDivBefore;
    finalVelocityCleanup.rmsDivAfter = projCleanup.rmsDivAfter;
    finalVelocityCleanup.maxAbsDivBefore = projCleanup.maxAbsDivBefore;
    finalVelocityCleanup.maxAbsDivAfter = projCleanup.maxAbsDivAfter;
    finalVelocityCleanup.divReduction = projCleanup.rmsDivAfter / max(projCleanup.rmsDivBefore, eps);
    finalVelocityCleanup.dvRms = sqrt(mean(sum((massFluxFinalVelocityProjectionStrength * dvCleanup).^2, 2), 'omitnan'));
    finalVelocityCleanup.proj = projCleanup;
end

dvTotal = stateOut.v - stateClassic.v;

if computeDiagnostics
    thermalAfterProjectionRaw = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', true, 'minCount', 1);
end

if thermostatAfterProjection
    [stateOut.v, thermostatInfo] = projection_apply_cell_thermostat(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', true);
else
    thermostatInfo = empty_thermostat_info();
end

if ~computeDiagnostics
    diag = minimal_diag(proj, massFluxProj, finalVelocityCleanup, thermostatInfo, projectionEnable, ...
        projectionStrength, massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, ...
        massFluxApplyAfterVelocityProjection, massFluxFinalVelocityProjectionCleanup, ...
        massFluxFinalVelocityProjectionStrength, projectionKind, appliedProjectionStrength, projectionInterpolationMethod, dvTotal);
    return;
end

Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
projAfterParticles = projection_project_grid_periodic_fft(Gafter.Ux, Gafter.Uy, params);
if useMassFluxProjection
    massFluxAfterParticles = projection_project_mass_flux_periodic_fft(Gafter.N, Gafter.Ux, Gafter.Uy, params, ...
        'mode', massFluxProjectionMode, ...
        'relaxationBeta', massFluxDensityRelaxationBeta);
else
    massFluxAfterParticles = empty_mass_flux_projection();
end
thermalAfterProjection = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
popAfterProjection = projection_population_diagnostics(stateOut.x, params, ...
    'periodicX', true, 'periodicY', true);
densityTransport = density_transport_periodic(Gbefore, Gafter, params);

momBefore = mean(stateClassic.v, 1, 'omitnan');
momAfter = mean(stateOut.v, 1, 'omitnan');

if useMassFluxProjection
    appliedRmsDivProjectedAfter = projAfterParticles.rmsDivBefore;
    appliedMaxAbsDivProjectedAfter = projAfterParticles.maxAbsDivBefore;
    appliedDivReductionProjected = massFluxProj.divMassReduction;
else
    appliedRmsDivProjectedAfter = proj.rmsDivAfter;
    appliedMaxAbsDivProjectedAfter = proj.maxAbsDivAfter;
    appliedDivReductionProjected = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
end

diag = struct();
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.appliedProjectionStrength = appliedProjectionStrength;
diag.projectionKind = projectionKind;
diag.massFluxProjectionMode = massFluxProjectionMode;
diag.massFluxProjectionStrength = massFluxProjectionStrength;
diag.massFluxDensityRelaxationBeta = massFluxDensityRelaxationBeta;
diag.massFluxApplyAfterVelocityProjection = massFluxApplyAfterVelocityProjection;
diag.massFluxFinalVelocityProjectionCleanup = massFluxFinalVelocityProjectionCleanup;
diag.massFluxFinalVelocityProjectionStrength = massFluxFinalVelocityProjectionStrength;
diag.finalVelocityCleanup = finalVelocityCleanup;
diag.finalVelocityCleanupRmsDivBefore = finalVelocityCleanup.rmsDivBefore;
diag.finalVelocityCleanupRmsDivAfter = finalVelocityCleanup.rmsDivAfter;
diag.finalVelocityCleanupDivReduction = finalVelocityCleanup.divReduction;
diag.finalVelocityCleanupDvRms = finalVelocityCleanup.dvRms;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivProjectedAfter = appliedRmsDivProjectedAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivProjectedAfter = appliedMaxAbsDivProjectedAfter;
diag.rmsDivParticleAfter = projAfterParticles.rmsDivBefore;
diag.maxAbsDivParticleAfter = projAfterParticles.maxAbsDivBefore;
diag.divReductionProjected = appliedDivReductionProjected;
diag.divReductionParticle = diag.rmsDivParticleAfter / max(proj.rmsDivBefore, eps);
diag.rmsMassFluxDivBefore = massFluxProj.rmsDivMassBefore;
diag.rmsMassFluxDivTarget = massFluxProj.rmsTargetDivMass;
diag.rmsMassFluxDivProjectedAfter = massFluxProj.rmsDivMassAfter;
diag.rmsMassFluxDivResidual = massFluxProj.rmsDivMassResidual;
diag.rmsMassFluxDivParticleAfter = massFluxAfterParticles.rmsDivMassBefore;
diag.maxAbsMassFluxDivResidual = massFluxProj.maxAbsDivMassResidual;
diag.massFluxDivReduction = massFluxProj.divMassReduction;
diag.massFluxProj = massFluxProj;
diag.massFluxAfterParticles = massFluxAfterParticles;
diag.nEmptyCellsBefore = nnz(Gbefore.N(:) == 0);
diag.nEmptyCellsAfter = nnz(Gafter.N(:) == 0);
diag.thermalBeforeProjection = thermalBeforeProjection;
diag.thermalAfterProjectionRaw = thermalAfterProjectionRaw;
diag.thermalAfterProjection = thermalAfterProjection;
diag.kBTCellBeforeProjection = thermalBeforeProjection.kBTCellRelative;
diag.kBTCellAfterProjectionRaw = thermalAfterProjectionRaw.kBTCellRelative;
diag.kBTCellAfterProjection = thermalAfterProjection.kBTCellRelative;
diag.kBTGlobalBeforeProjection = thermalBeforeProjection.kBTGlobal;
diag.kBTGlobalAfterProjection = thermalAfterProjection.kBTGlobal;
diag.hydroKEBeforeProjection = thermalBeforeProjection.hydroKineticEnergy;
diag.hydroKEAfterProjection = thermalAfterProjection.hydroKineticEnergy;
diag.totalKEBeforeProjection = thermalBeforeProjection.totalKineticEnergy;
diag.totalKEAfterProjection = thermalAfterProjection.totalKineticEnergy;
diag.thermostatAfterProjection = thermostatAfterProjection;
diag.thermostatInfo = thermostatInfo;
diag.thermostatRmsVelocityChange = thermostatInfo.rmsVelocityChange;
diag.thermostatMeanScale = thermostatInfo.meanScale;
diag.thermostatNCells = thermostatInfo.nThermostattedCells;
diag.meanVxBeforeProjection = momBefore(1);
diag.meanVyBeforeProjection = momBefore(2);
diag.meanVxAfterProjection = momAfter(1);
diag.meanVyAfterProjection = momAfter(2);
diag.deltaMeanVx = momAfter(1) - momBefore(1);
diag.deltaMeanVy = momAfter(2) - momBefore(2);
diag.dvRms = sqrt(mean(sum(dvTotal.^2, 2), 'omitnan'));
diag.Gbefore = Gbefore;
diag.Gafter = Gafter;
diag.proj = proj;
diag.populationAfterClassic = popAfterClassic;
diag.populationAfterProjection = popAfterProjection;
diag.densityTransport = densityTransport;
diag.densityTransportClassicRms = densityTransport.classic.rms;
diag.densityTransportProjectedRms = densityTransport.projected.rms;
diag.densityTransportProjectedMinusClassicRms = densityTransport.projectedMinusClassic.rms;
end

function tr = density_transport_periodic(Gclassic, Gprojected, params)
N = double(Gclassic.N);
fluxClassicX = N .* Gclassic.Ux;
fluxClassicY = N .* Gclassic.Uy;
fluxProjectedX = N .* Gprojected.Ux;
fluxProjectedY = N .* Gprojected.Uy;
divClassic = divergence_periodic(fluxClassicX, fluxClassicY, params);
divProjected = divergence_periodic(fluxProjectedX, fluxProjectedY, params);
deltaClassic = -params.dt * divClassic;
deltaProjected = -params.dt * divProjected;
deltaDiff = deltaProjected - deltaClassic;
tr = struct();
tr.deltaClassic = deltaClassic;
tr.deltaProjected = deltaProjected;
tr.deltaProjectedMinusClassic = deltaDiff;
tr.divFluxClassic = divClassic;
tr.divFluxProjected = divProjected;
tr.classic = summarize_delta(deltaClassic, N);
tr.projected = summarize_delta(deltaProjected, N);
tr.projectedMinusClassic = summarize_delta(deltaDiff, N);
tr.massDeltaClassic = sum(deltaClassic(:));
tr.massDeltaProjected = sum(deltaProjected(:));
tr.massDeltaProjectedMinusClassic = sum(deltaDiff(:));
end

function div = divergence_periodic(Fx, Fy, params)
[Nx, Ny] = size(Fx);
dx = params.Lx / Nx;
dy = params.Ly / Ny;
ixp = [2:Nx, 1];
ixm = [Nx, 1:Nx-1];
iyp = [2:Ny, 1];
iym = [Ny, 1:Ny-1];
dFx = (Fx(ixp, :) - Fx(ixm, :)) / (2*dx);
dFy = (Fy(:, iyp) - Fy(:, iym)) / (2*dy);
div = dFx + dFy;
end

function s = summarize_delta(delta, N)
d = delta(:);
s = struct();
s.rms = sqrt(mean(d.^2, 'omitnan'));
s.maxAbs = max(abs(d));
s.meanAbs = mean(abs(d), 'omitnan');
s.std = std(d, 'omitnan');
s.sum = sum(d, 'omitnan');
s.rmsOverGamma = s.rms / max(mean(N(:), 'omitnan'), eps);
s.maxAbsOverGamma = s.maxAbs / max(mean(N(:), 'omitnan'), eps);
s.predictedStdAfter = std(N(:) + d, 'omitnan');
s.predictedStdDelta = s.predictedStdAfter - std(N(:), 'omitnan');
end

function diag = minimal_diag(proj, massFluxProj, finalVelocityCleanup, thermostatInfo, projectionEnable, projectionStrength, massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, massFluxApplyAfterVelocityProjection, massFluxFinalVelocityProjectionCleanup, massFluxFinalVelocityProjectionStrength, projectionKind, appliedProjectionStrength, projectionInterpolationMethod, dvTotal)
diag = struct();
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.appliedProjectionStrength = appliedProjectionStrength;
diag.projectionKind = projectionKind;
diag.massFluxProjectionMode = massFluxProjectionMode;
diag.massFluxProjectionStrength = massFluxProjectionStrength;
diag.massFluxDensityRelaxationBeta = massFluxDensityRelaxationBeta;
diag.massFluxApplyAfterVelocityProjection = massFluxApplyAfterVelocityProjection;
diag.massFluxFinalVelocityProjectionCleanup = massFluxFinalVelocityProjectionCleanup;
diag.massFluxFinalVelocityProjectionStrength = massFluxFinalVelocityProjectionStrength;
diag.finalVelocityCleanup = finalVelocityCleanup;
diag.finalVelocityCleanupRmsDivBefore = finalVelocityCleanup.rmsDivBefore;
diag.finalVelocityCleanupRmsDivAfter = finalVelocityCleanup.rmsDivAfter;
diag.finalVelocityCleanupDivReduction = finalVelocityCleanup.divReduction;
diag.finalVelocityCleanupDvRms = finalVelocityCleanup.dvRms;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivProjectedAfter = proj.rmsDivAfter;
diag.rmsDivParticleAfter = NaN;
diag.rmsMassFluxDivBefore = massFluxProj.rmsDivMassBefore;
diag.rmsMassFluxDivResidual = massFluxProj.rmsDivMassResidual;
diag.rmsMassFluxDivParticleAfter = NaN;
diag.massFluxDivReduction = massFluxProj.divMassReduction;
diag.dvRms = sqrt(mean(sum(dvTotal.^2, 2), 'omitnan'));
diag.thermostatAfterProjection = thermostatInfo.enabled;
diag.thermostatInfo = thermostatInfo;
diag.thermostatRmsVelocityChange = thermostatInfo.rmsVelocityChange;
diag.thermostatMeanScale = thermostatInfo.meanScale;
diag.thermostatNCells = thermostatInfo.nThermostattedCells;
diag.kBTCellAfterProjection = NaN;
diag.densityTransportProjectedRms = NaN;
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v')
    error('state must be a struct with fields x and v.');
end
if size(state.x, 2) ~= 2 || size(state.v, 2) ~= 2 || size(state.x, 1) ~= size(state.v, 1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
end

function info = empty_velocity_cleanup()
info = struct();
info.enabled = false;
info.strength = 0.0;
info.rmsDivBefore = NaN;
info.rmsDivAfter = NaN;
info.maxAbsDivBefore = NaN;
info.maxAbsDivAfter = NaN;
info.divReduction = NaN;
info.dvRms = NaN;
info.proj = [];
end

function info = empty_mass_flux_projection()
info = struct();
info.mode = 'off';
info.beta = 0.0;
info.rmsDivMassBefore = NaN;
info.rmsTargetDivMass = NaN;
info.rmsDivMassAfter = NaN;
info.rmsDivMassResidual = NaN;
info.rmsDivMassBeforeFull = NaN;
info.rmsDivMassAfterFull = NaN;
info.rmsDivMassResidualFull = NaN;
info.maxAbsDivMassBefore = NaN;
info.maxAbsDivMassAfter = NaN;
info.maxAbsDivMassResidual = NaN;
info.divMassReduction = NaN;
end

function info = empty_thermostat_info()
info = struct('enabled', false, 'targetKBT', 0.0, 'strength', 0.0, ...
    'minParticlesPerCell', 0, 'maxScale', 1.0, 'nThermostattedCells', 0, ...
    'meanScale', 1.0, 'minScale', 1.0, 'maxScaleApplied', 1.0, ...
    'meanKBTBefore', 0.0, 'meanKBTAfter', 0.0, 'rmsVelocityChange', 0.0);
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
