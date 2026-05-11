function [stateOut, diag] = mpcd_step_projection_poiseuille(state, params)
%MPCD_STEP_PROJECTION_POISEUILLE MPCD channel step plus pressure projection.
%
%   [stateOut, diag] = mpcd_step_projection_poiseuille(state, params)
%
% Prototype for Poiseuille flow:
%   x periodic, y walls, bodyForceX, classic SRD collision, then algebraic
%   pressure projection on a periodic-x / bounded-y grid.
%
% The projection changes particle velocities only. It must not alter the cell
% population distribution; diagnostics report population before the step,
% after classic MPCD, and after projection.

projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
massFluxProjectionMode = char(string(get_param(params, 'massFluxProjectionMode', 'off')));
massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', projectionStrength);
massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 0.0);
massFluxApplyAfterVelocityProjection = logical(get_param(params, 'massFluxApplyAfterVelocityProjection', false));
useMassFluxProjection = ~strcmpi(strrep(massFluxProjectionMode, '-', '_'), 'off');
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
projectionTransportDiagnosticsEnable = logical(get_param(params, 'projectionTransportDiagnosticsEnable', true));
densityTransportDiagnosticsEnable = logical(get_param(params, 'densityTransportDiagnosticsEnable', true));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', false));
computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));

if computeDiagnostics
    popBeforeStep = projection_population_diagnostics(state.x, params, ...
        'periodicX', true, 'periodicY', false);
end

[stateClassic, classicDiag] = mpcd_step_classic_poiseuille(state, params);

if computeDiagnostics
    popAfterClassic = projection_population_diagnostics(stateClassic.x, params, ...
        'periodicX', true, 'periodicY', false);
end

Gbefore = projection_deposit_particles_to_grid(stateClassic.x, stateClassic.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
if computeDiagnostics
    thermalBeforeProjection = projection_thermal_diagnostics(stateClassic.x, stateClassic.v, params, ...
        'periodicX', true, 'periodicY', false, 'minCount', 1);
end

proj = projection_project_grid_periodic_x_neumann_y(Gbefore.Ux, Gbefore.Uy, params);
massFluxProj = empty_mass_flux_projection();
projectionKind = 'velocity';
appliedProjectionStrength = projectionStrength;

stateOut = stateClassic;
if useMassFluxProjection
    appliedProjectionStrength = massFluxProjectionStrength;

    if massFluxApplyAfterVelocityProjection
        projectionKind = 'velocity_plus_mass_flux';
        stateForMassFlux = stateClassic;
        if projectionEnable && projectionStrength ~= 0
            dvVelocity = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
                'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
            stateForMassFlux.v = stateClassic.v + projectionStrength * dvVelocity;
        else
            dvVelocity = zeros(size(stateClassic.v));
        end
        GforMassFlux = projection_deposit_particles_to_grid(stateForMassFlux.x, stateForMassFlux.v, params, ...
            'periodicX', true, 'periodicY', false, 'minCount', 1);
        massFluxProj = projection_project_mass_flux_periodic_x_neumann_y(GforMassFlux.N, GforMassFlux.Ux, GforMassFlux.Uy, params, ...
            'mode', massFluxProjectionMode, ...
            'relaxationBeta', massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dvMassFlux = projection_interpolate_grid_delta_to_particles(stateForMassFlux.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
            stateOut.v = stateForMassFlux.v + massFluxProjectionStrength * dvMassFlux;
        else
            dvMassFlux = zeros(size(stateClassic.v));
            stateOut = stateForMassFlux;
        end
        dv = stateOut.v - stateClassic.v;
    else
        projectionKind = 'mass_flux';
        massFluxProj = projection_project_mass_flux_periodic_x_neumann_y(Gbefore.N, Gbefore.Ux, Gbefore.Uy, params, ...
            'mode', massFluxProjectionMode, ...
            'relaxationBeta', massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
            stateOut.v = stateClassic.v + massFluxProjectionStrength * dv;
        else
            dv = zeros(size(stateClassic.v));
        end
    end
elseif projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
    stateOut.v = stateClassic.v + projectionStrength * dv;
else
    dv = zeros(size(stateClassic.v));
end

if computeDiagnostics
    thermalAfterProjectionRaw = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', false, 'minCount', 1);
end

if thermostatAfterProjection
    [stateOut.v, thermostatInfo] = projection_apply_cell_thermostat(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', false);
else
    thermostatInfo = empty_thermostat_info();
end

if ~computeDiagnostics
    diag = minimal_projection_diag(classicDiag, proj, massFluxProj, thermostatInfo, projectionEnable, ...
        projectionStrength, projectionInterpolationMethod, dv, projectionKind, ...
        massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, appliedProjectionStrength, massFluxApplyAfterVelocityProjection);
    return;
end

thermalAfterProjection = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);

popAfterProjection = projection_population_diagnostics(stateOut.x, params, ...
    'periodicX', true, 'periodicY', false);

Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
projAfterParticles = projection_project_grid_periodic_x_neumann_y(Gafter.Ux, Gafter.Uy, params);
if useMassFluxProjection
    massFluxAfterParticles = projection_project_mass_flux_periodic_x_neumann_y(Gafter.N, Gafter.Ux, Gafter.Uy, params, ...
        'mode', massFluxProjectionMode, ...
        'relaxationBeta', massFluxDensityRelaxationBeta);
else
    massFluxAfterParticles = empty_mass_flux_projection();
end

if densityTransportDiagnosticsEnable
    densityTransport = projection_density_transport_continuous_diagnostics(Gbefore, Gafter, params);
else
    densityTransport = empty_density_transport_diag();
end

momBefore = mean(stateClassic.v, 1);
momAfter = mean(stateOut.v, 1);

popDeltaProjection = double(popAfterProjection.N) - double(popAfterClassic.N);

if projectionTransportDiagnosticsEnable
    populationTransport = projection_population_transport_diagnostics(stateClassic.x, stateClassic.v, stateOut.v, params, ...
        'periodicX', true, 'periodicY', false, 'includeBodyForce', true);
else
    populationTransport = empty_population_transport_diag();
end

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
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.appliedProjectionStrength = appliedProjectionStrength;
diag.projectionKind = projectionKind;
diag.massFluxProjectionMode = massFluxProjectionMode;
diag.massFluxProjectionStrength = massFluxProjectionStrength;
diag.massFluxDensityRelaxationBeta = massFluxDensityRelaxationBeta;
diag.massFluxApplyAfterVelocityProjection = massFluxApplyAfterVelocityProjection;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.computeDiagnostics = computeDiagnostics;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivProjectedAfter = appliedRmsDivProjectedAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivProjectedAfter = appliedMaxAbsDivProjectedAfter;
diag.rmsDivParticleAfter = projAfterParticles.rmsDivBefore;
diag.maxAbsDivParticleAfter = projAfterParticles.maxAbsDivBefore;
diag.divReductionProjected = appliedDivReductionProjected;
diag.divReductionParticle = diag.rmsDivParticleAfter / max(proj.rmsDivBefore, eps);
diag.meanDivBefore = proj.meanDivBefore;
diag.rmsMassFluxDivBefore = massFluxProj.rmsDivMassBefore;
diag.rmsMassFluxDivTarget = massFluxProj.rmsTargetDivMass;
diag.rmsMassFluxDivProjectedAfter = massFluxProj.rmsDivMassAfter;
diag.rmsMassFluxDivResidual = massFluxProj.rmsDivMassResidual;
diag.rmsMassFluxDivParticleAfter = massFluxAfterParticles.rmsDivMassBefore;
diag.maxAbsMassFluxDivBefore = massFluxProj.maxAbsDivMassBefore;
diag.maxAbsMassFluxDivProjectedAfter = massFluxProj.maxAbsDivMassAfter;
diag.maxAbsMassFluxDivResidual = massFluxProj.maxAbsDivMassResidual;
diag.massFluxDivReduction = massFluxProj.divMassReduction;
diag.massFluxProj = massFluxProj;
diag.massFluxAfterParticles = massFluxAfterParticles;
diag.nEmptyCellsBefore = nnz(Gbefore.N(:) == 0);
diag.nEmptyCellsAfter = nnz(Gafter.N(:) == 0);
diag.kBTBeforeProjection = estimate_kBT(stateClassic.v);
diag.kBTAfterProjection = estimate_kBT(stateOut.v);
diag.thermalBeforeProjection = thermalBeforeProjection;
diag.thermalAfterProjectionRaw = thermalAfterProjectionRaw;
diag.thermalAfterProjection = thermalAfterProjection;
diag.kBTGlobalBeforeProjection = thermalBeforeProjection.kBTGlobal;
diag.kBTGlobalAfterProjectionRaw = thermalAfterProjectionRaw.kBTGlobal;
diag.kBTGlobalAfterProjection = thermalAfterProjection.kBTGlobal;
diag.kBTCellBeforeProjection = thermalBeforeProjection.kBTCellRelative;
diag.kBTCellAfterProjectionRaw = thermalAfterProjectionRaw.kBTCellRelative;
diag.kBTCellAfterProjection = thermalAfterProjection.kBTCellRelative;
diag.hydroKEBeforeProjection = thermalBeforeProjection.hydroKineticEnergy;
diag.hydroKEAfterProjectionRaw = thermalAfterProjectionRaw.hydroKineticEnergy;
diag.hydroKEAfterProjection = thermalAfterProjection.hydroKineticEnergy;
diag.thermalKEBeforeProjection = thermalBeforeProjection.thermalKineticEnergy;
diag.thermalKEAfterProjectionRaw = thermalAfterProjectionRaw.thermalKineticEnergy;
diag.thermalKEAfterProjection = thermalAfterProjection.thermalKineticEnergy;
diag.totalKEBeforeProjection = thermalBeforeProjection.totalKineticEnergy;
diag.totalKEAfterProjectionRaw = thermalAfterProjectionRaw.totalKineticEnergy;
diag.totalKEAfterProjection = thermalAfterProjection.totalKineticEnergy;
diag.projectionThermalDeltaRaw = thermalAfterProjectionRaw.kBTCellRelative - thermalBeforeProjection.kBTCellRelative;
diag.projectionThermalDeltaAfterThermostat = thermalAfterProjection.kBTCellRelative - thermalBeforeProjection.kBTCellRelative;
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
diag.dvRms = sqrt(mean(sum(dv.^2, 2)));
diag.Gbefore = Gbefore;
diag.Gafter = Gafter;
diag.proj = proj;
diag.populationBeforeStep = popBeforeStep;
diag.populationAfterClassic = popAfterClassic;
diag.populationAfterProjection = popAfterProjection;
diag.populationProjectionDeltaMaxAbs = max(abs(popDeltaProjection(:)));
diag.populationProjectionDeltaRms = sqrt(mean(popDeltaProjection(:).^2));
diag.populationStepDeltaStd = popAfterClassic.stdN - popBeforeStep.stdN;
diag.populationStepDeltaEmpty = popAfterClassic.nEmptyCells - popBeforeStep.nEmptyCells;
diag.populationTransportDiagnosticsEnable = projectionTransportDiagnosticsEnable;
diag.populationTransport = populationTransport;
diag.populationTransportClassicDeltaRms = populationTransport.classic.rms;
diag.populationTransportProjectedDeltaRms = populationTransport.projected.rms;
diag.populationTransportProjectedMinusClassicRms = populationTransport.projectedMinusClassic.rms;
diag.populationTransportClassicDeltaMaxAbs = populationTransport.classic.maxAbs;
diag.populationTransportProjectedDeltaMaxAbs = populationTransport.projected.maxAbs;
diag.populationTransportProjectedMinusClassicMaxAbs = populationTransport.projectedMinusClassic.maxAbs;
diag.populationTransportClassicStdAfter = populationTransport.classic.stdAfter;
diag.populationTransportProjectedStdAfter = populationTransport.projected.stdAfter;
diag.populationTransportProjectedVsClassicStdDelta = populationTransport.projectedMinusClassic.stdDelta;
diag.populationTransportClassicEmptyAfter = populationTransport.classic.emptyAfter;
diag.populationTransportProjectedEmptyAfter = populationTransport.projected.emptyAfter;
diag.populationTransportProjectedVsClassicEmptyDelta = populationTransport.projectedMinusClassic.emptyDelta;
diag.populationTransportMassErrorClassic = populationTransport.massErrorClassic;
diag.populationTransportMassErrorProjected = populationTransport.massErrorProjected;
diag.populationTransportMassErrorProjectedMinusClassic = populationTransport.massErrorProjectedMinusClassic;
diag.densityTransportDiagnosticsEnable = densityTransportDiagnosticsEnable;
diag.densityTransport = densityTransport;
diag.densityTransportClassicRms = densityTransport.classic.rms;
diag.densityTransportProjectedRms = densityTransport.projected.rms;
diag.densityTransportProjectedMinusClassicRms = densityTransport.projectedMinusClassic.rms;
diag.densityTransportClassicMaxAbs = densityTransport.classic.maxAbs;
diag.densityTransportProjectedMaxAbs = densityTransport.projected.maxAbs;
diag.densityTransportProjectedMinusClassicMaxAbs = densityTransport.projectedMinusClassic.maxAbs;
diag.densityTransportClassicPredictedStdAfter = densityTransport.classic.predictedStdAfter;
diag.densityTransportProjectedPredictedStdAfter = densityTransport.projected.predictedStdAfter;
diag.densityTransportProjectedVsClassicStdDelta = densityTransport.projected.predictedStdAfter - densityTransport.classic.predictedStdAfter;
diag.densityTransportMassDeltaClassic = densityTransport.massDeltaClassic;
diag.densityTransportMassDeltaProjected = densityTransport.massDeltaProjected;
diag.densityTransportMassDeltaProjectedMinusClassic = densityTransport.massDeltaProjectedMinusClassic;
end



function diag = minimal_projection_diag(classicDiag, proj, massFluxProj, thermostatInfo, projectionEnable, projectionStrength, projectionInterpolationMethod, dv, projectionKind, massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, appliedProjectionStrength, massFluxApplyAfterVelocityProjection)
% Minimal finite diagnostics for fast steps. Expensive diagnostics such as a
% second particle-grid projection, thermal diagnostics, and population-transport
% forecasts are intentionally skipped. run_projection_poiseuille_demo enables
% full diagnostics only on sampled/progress/final steps.
diag = struct();
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.appliedProjectionStrength = appliedProjectionStrength;
diag.projectionKind = projectionKind;
diag.massFluxProjectionMode = massFluxProjectionMode;
diag.massFluxProjectionStrength = massFluxProjectionStrength;
diag.massFluxDensityRelaxationBeta = massFluxDensityRelaxationBeta;
diag.massFluxApplyAfterVelocityProjection = massFluxApplyAfterVelocityProjection;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.computeDiagnostics = false;
diag.rmsDivBefore = proj.rmsDivBefore;
if ~isempty(strfind(projectionKind, 'mass_flux'))
    diag.rmsDivProjectedAfter = NaN;
    diag.maxAbsDivProjectedAfter = NaN;
    diag.divReductionProjected = massFluxProj.divMassReduction;
else
    diag.rmsDivProjectedAfter = proj.rmsDivAfter;
    diag.maxAbsDivProjectedAfter = proj.maxAbsDivAfter;
    diag.divReductionProjected = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
end
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.rmsDivParticleAfter = NaN;
diag.maxAbsDivParticleAfter = NaN;
diag.divReductionParticle = NaN;
diag.meanDivBefore = proj.meanDivBefore;
diag.rmsMassFluxDivBefore = massFluxProj.rmsDivMassBefore;
diag.rmsMassFluxDivTarget = massFluxProj.rmsTargetDivMass;
diag.rmsMassFluxDivProjectedAfter = massFluxProj.rmsDivMassAfter;
diag.rmsMassFluxDivResidual = massFluxProj.rmsDivMassResidual;
diag.rmsMassFluxDivParticleAfter = NaN;
diag.maxAbsMassFluxDivBefore = massFluxProj.maxAbsDivMassBefore;
diag.maxAbsMassFluxDivProjectedAfter = massFluxProj.maxAbsDivMassAfter;
diag.maxAbsMassFluxDivResidual = massFluxProj.maxAbsDivMassResidual;
diag.massFluxDivReduction = massFluxProj.divMassReduction;
diag.massFluxProj = massFluxProj;
diag.dvRms = sqrt(mean(sum(dv.^2, 2)));
diag.thermostatAfterProjection = thermostatInfo.enabled;
diag.thermostatInfo = thermostatInfo;
diag.thermostatRmsVelocityChange = thermostatInfo.rmsVelocityChange;
diag.thermostatMeanScale = thermostatInfo.meanScale;
diag.thermostatNCells = thermostatInfo.nThermostattedCells;
diag.kBTCellAfterProjection = NaN;
end


function info = empty_mass_flux_projection()
info = struct();
info.mode = 'off';
info.beta = 0.0;
info.rmsDivMassBefore = NaN;
info.rmsTargetDivMass = NaN;
info.rmsDivMassAfter = NaN;
info.rmsDivMassResidual = NaN;
info.maxAbsDivMassBefore = NaN;
info.maxAbsDivMassAfter = NaN;
info.maxAbsDivMassResidual = NaN;
info.divMassReduction = NaN;
end

function info = empty_thermostat_info()
% Disabled thermostat = identity operation. Keep all scalar fields finite so
% diagnostic histories can be checked with isfinite(H(:)) even when the
% optional thermostat is off.
info = struct('enabled', false, 'targetKBT', 0.0, 'strength', 0.0, ...
    'minParticlesPerCell', 0, 'maxScale', 1.0, 'nThermostattedCells', 0, ...
    'meanScale', 1.0, 'minScale', 1.0, 'maxScaleApplied', 1.0, ...
    'meanKBTBefore', 0.0, 'meanKBTAfter', 0.0, 'rmsVelocityChange', 0.0);
end

function tr = empty_density_transport_diag()
emptySummary = struct('rms', NaN, 'maxAbs', NaN, 'meanAbs', NaN, 'std', NaN, ...
    'sum', NaN, 'rmsOverGamma', NaN, 'maxAbsOverGamma', NaN, ...
    'predictedStdAfter', NaN, 'predictedStdDelta', NaN);
tr = struct();
tr.classic = emptySummary;
tr.projected = emptySummary;
tr.projectedMinusClassic = emptySummary;
tr.massDeltaClassic = NaN;
tr.massDeltaProjected = NaN;
tr.massDeltaProjectedMinusClassic = NaN;
end

function tr = empty_population_transport_diag()
emptySummary = struct('rms', NaN, 'maxAbs', NaN, 'meanAbs', NaN, 'std', NaN, 'sum', NaN, ...
    'rmsOverGamma', NaN, 'maxAbsOverGamma', NaN, 'stdBefore', NaN, 'stdAfter', NaN, ...
    'stdDelta', NaN, 'emptyBefore', NaN, 'emptyAfter', NaN, 'emptyDelta', NaN, ...
    'outBandBefore', NaN, 'outBandAfter', NaN, 'outBandDelta', NaN, ...
    'cvBefore', NaN, 'cvAfter', NaN, 'cvDelta', NaN);
tr = struct();
tr.classic = emptySummary;
tr.projected = emptySummary;
tr.projectedMinusClassic = emptySummary;
tr.massErrorClassic = NaN;
tr.massErrorProjected = NaN;
tr.massErrorProjectedMinusClassic = NaN;
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end

function kBT = estimate_kBT(v)
u = mean(v, 1);
c = v - u;
kBT = 0.5 * mean(sum(c.^2, 2));
end
