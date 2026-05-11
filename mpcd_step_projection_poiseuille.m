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
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
projectionTransportDiagnosticsEnable = logical(get_param(params, 'projectionTransportDiagnosticsEnable', true));
densityTransportDiagnosticsEnable = logical(get_param(params, 'densityTransportDiagnosticsEnable', true));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', false));
useVirialDensityKick = logical(get_param(params, 'useVirialDensityKick', false));
useVirialDensityRepair = logical(get_param(params, 'useVirialDensityRepair', false));
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

stateOut = stateClassic;
if projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
    stateOut.v = stateClassic.v + projectionStrength * dv;
else
    dv = zeros(size(stateClassic.v));
end
stateAfterPressureProjection = stateOut;

GafterPressureForRepair = projection_deposit_particles_to_grid(stateAfterPressureProjection.x, stateAfterPressureProjection.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
if useVirialDensityRepair
    [stateOut, densityRepairInfo] = projection_apply_virial_density_repair(stateOut, params, ...
        'periodicX', true, 'periodicY', false, 'targetGrid', GafterPressureForRepair);
else
    densityRepairInfo = empty_density_repair_info();
end
stateAfterDensityRepair = stateOut;

if useVirialDensityKick
    [stateOut.v, virialInfo] = projection_apply_virial_density_kick(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', false);
else
    virialInfo = empty_virial_info();
end
stateAfterVirialRaw = stateOut;

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
    diag = minimal_projection_diag(classicDiag, proj, thermostatInfo, virialInfo, densityRepairInfo, projectionEnable, ...
        projectionStrength, projectionInterpolationMethod, dv);
    return;
end

thermalAfterProjection = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);

popAfterProjection = projection_population_diagnostics(stateOut.x, params, ...
    'periodicX', true, 'periodicY', false);

GafterPressureProjection = GafterPressureForRepair;
projAfterPressureParticles = projection_project_grid_periodic_x_neumann_y(GafterPressureProjection.Ux, GafterPressureProjection.Uy, params);

GafterDensityRepair = projection_deposit_particles_to_grid(stateAfterDensityRepair.x, stateAfterDensityRepair.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
projAfterDensityRepair = projection_project_grid_periodic_x_neumann_y(GafterDensityRepair.Ux, GafterDensityRepair.Uy, params);

GafterVirialRaw = projection_deposit_particles_to_grid(stateAfterVirialRaw.x, stateAfterVirialRaw.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
projAfterVirialRaw = projection_project_grid_periodic_x_neumann_y(GafterVirialRaw.Ux, GafterVirialRaw.Uy, params);

Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
projAfterParticles = projection_project_grid_periodic_x_neumann_y(Gafter.Ux, Gafter.Uy, params);

if densityTransportDiagnosticsEnable
    densityTransportPressureOnly = projection_density_transport_continuous_diagnostics(Gbefore, GafterPressureProjection, params);
    densityTransport = projection_density_transport_continuous_diagnostics(Gbefore, Gafter, params);
else
    densityTransportPressureOnly = empty_density_transport_diag();
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

diag = struct();
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.computeDiagnostics = computeDiagnostics;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivProjectedAfter = proj.rmsDivAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivProjectedAfter = proj.maxAbsDivAfter;
diag.rmsDivParticleAfter = projAfterParticles.rmsDivBefore;
diag.maxAbsDivParticleAfter = projAfterParticles.maxAbsDivBefore;
diag.divReductionProjected = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
diag.divReductionParticle = diag.rmsDivParticleAfter / max(proj.rmsDivBefore, eps);
diag.meanDivBefore = proj.meanDivBefore;
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
diag.virialDensityKickEnabled = virialInfo.enabled;
diag.virialInfo = virialInfo;
diag.virialDensityKickStrength = virialInfo.strength;
diag.virialK = virialInfo.Kvir;
diag.virialSmoothPasses = virialInfo.smoothPasses;
diag.virialParticleDVRms = virialInfo.particleDVRms;
diag.virialParticleDVMaxAbs = virialInfo.particleDVMaxAbs;
diag.virialPvirRms = virialInfo.pVirRms;
diag.virialPvirMaxAbs = virialInfo.pVirMaxAbs;
diag.virialDensityRepairEnabled = densityRepairInfo.enabled;
diag.densityRepairInfo = densityRepairInfo;
diag.densityRepairStrength = densityRepairInfo.strength;
diag.densityRepairK = densityRepairInfo.Kvir;
diag.densityRepairSmoothPasses = densityRepairInfo.smoothPasses;
diag.densityRepairDisplacementRms = densityRepairInfo.particleDisplacementRms;
diag.densityRepairDisplacementMaxAbs = densityRepairInfo.particleDisplacementMaxAbs;
diag.densityRepairStdBefore = densityRepairInfo.stdNBefore;
diag.densityRepairStdAfter = densityRepairInfo.stdNAfter;
diag.densityRepairOutBand20Before = densityRepairInfo.outBand20Before;
diag.densityRepairOutBand20After = densityRepairInfo.outBand20After;
diag.densityRepairVelocityRestoreDeltaRms = densityRepairInfo.velocityRestoreDeltaRms;
diag.densityRepairVelocityRestoreResidualBeforeRms = densityRepairInfo.velocityRestoreResidualBeforeRms;
diag.densityRepairVelocityRestoreResidualAfterRms = densityRepairInfo.velocityRestoreResidualAfterRms;
diag.rmsDivAfterPressureParticleRaw = projAfterPressureParticles.rmsDivBefore;
diag.rmsDivAfterDensityRepair = projAfterDensityRepair.rmsDivBefore;
diag.maxAbsDivAfterDensityRepair = projAfterDensityRepair.maxAbsDivBefore;
diag.maxAbsDivAfterPressureParticleRaw = projAfterPressureParticles.maxAbsDivBefore;
diag.rmsDivAfterVirialRaw = projAfterVirialRaw.rmsDivBefore;
diag.maxAbsDivAfterVirialRaw = projAfterVirialRaw.maxAbsDivBefore;
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
diag.GafterPressureProjection = GafterPressureProjection;
diag.GafterDensityRepair = GafterDensityRepair;
diag.GafterVirialRaw = GafterVirialRaw;
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
diag.densityTransportPressureOnly = densityTransportPressureOnly;
diag.densityTransport = densityTransport;
diag.densityTransportPressureOnlyProjectedRms = densityTransportPressureOnly.projected.rms;
diag.densityTransportPressureOnlyProjectedMaxAbs = densityTransportPressureOnly.projected.maxAbs;
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



function diag = minimal_projection_diag(classicDiag, proj, thermostatInfo, virialInfo, densityRepairInfo, projectionEnable, projectionStrength, projectionInterpolationMethod, dv)
% Minimal finite diagnostics for fast steps. Expensive diagnostics such as a
% second particle-grid projection, thermal diagnostics, and population-transport
% forecasts are intentionally skipped. run_projection_poiseuille_demo enables
% full diagnostics only on sampled/progress/final steps.
diag = struct();
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
diag.computeDiagnostics = false;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivProjectedAfter = proj.rmsDivAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivProjectedAfter = proj.maxAbsDivAfter;
diag.rmsDivParticleAfter = NaN;
diag.maxAbsDivParticleAfter = NaN;
diag.divReductionProjected = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
diag.divReductionParticle = NaN;
diag.meanDivBefore = proj.meanDivBefore;
diag.dvRms = sqrt(mean(sum(dv.^2, 2)));
diag.thermostatAfterProjection = thermostatInfo.enabled;
diag.virialDensityKickEnabled = virialInfo.enabled;
diag.virialInfo = virialInfo;
diag.virialDensityKickStrength = virialInfo.strength;
diag.virialK = virialInfo.Kvir;
diag.virialParticleDVRms = virialInfo.particleDVRms;
diag.virialParticleDVMaxAbs = virialInfo.particleDVMaxAbs;
diag.virialPvirRms = virialInfo.pVirRms;
diag.virialPvirMaxAbs = virialInfo.pVirMaxAbs;
diag.virialDensityRepairEnabled = densityRepairInfo.enabled;
diag.densityRepairInfo = densityRepairInfo;
diag.densityRepairStrength = densityRepairInfo.strength;
diag.densityRepairK = densityRepairInfo.Kvir;
diag.densityRepairSmoothPasses = densityRepairInfo.smoothPasses;
diag.densityRepairDisplacementRms = densityRepairInfo.particleDisplacementRms;
diag.densityRepairDisplacementMaxAbs = densityRepairInfo.particleDisplacementMaxAbs;
diag.densityRepairStdBefore = densityRepairInfo.stdNBefore;
diag.densityRepairStdAfter = densityRepairInfo.stdNAfter;
diag.densityRepairOutBand20Before = densityRepairInfo.outBand20Before;
diag.densityRepairOutBand20After = densityRepairInfo.outBand20After;
diag.densityRepairVelocityRestoreDeltaRms = densityRepairInfo.velocityRestoreDeltaRms;
diag.densityRepairVelocityRestoreResidualBeforeRms = densityRepairInfo.velocityRestoreResidualBeforeRms;
diag.densityRepairVelocityRestoreResidualAfterRms = densityRepairInfo.velocityRestoreResidualAfterRms;
diag.rmsDivAfterPressureParticleRaw = NaN;
diag.rmsDivAfterDensityRepair = NaN;
diag.maxAbsDivAfterPressureParticleRaw = NaN;
diag.rmsDivAfterVirialRaw = NaN;
diag.maxAbsDivAfterVirialRaw = NaN;
diag.densityTransportPressureOnlyProjectedRms = NaN;
diag.densityTransportPressureOnlyProjectedMaxAbs = NaN;
diag.thermostatInfo = thermostatInfo;
diag.thermostatRmsVelocityChange = thermostatInfo.rmsVelocityChange;
diag.thermostatMeanScale = thermostatInfo.meanScale;
diag.thermostatNCells = thermostatInfo.nThermostattedCells;
diag.kBTCellAfterProjection = NaN;
end

function info = empty_virial_info()
info = struct();
info.enabled = false;
info.requested = false;
info.strength = 0.0;
info.Kvir = 0.0;
info.smoothPasses = 0;
info.minCellCount = 1.0;
info.maxParticleKick = Inf;
info.method = '';
info.meanN = NaN;
info.stdN = NaN;
info.minN = NaN;
info.maxN = NaN;
info.pVirMean = 0.0;
info.pVirRms = 0.0;
info.pVirMaxAbs = 0.0;
info.gradPvirRms = 0.0;
info.gradPvirMaxAbs = 0.0;
info.gridDVRms = 0.0;
info.gridDVMaxAbs = 0.0;
info.particleDVRms = 0.0;
info.particleDVMaxAbs = 0.0;
info.particleDVMeanX = 0.0;
info.particleDVMeanY = 0.0;
info.nLimitedParticles = 0;
info.N = [];
info.Pvir = [];
info.gradPvirX = [];
info.gradPvirY = [];
info.dUx = [];
info.dUy = [];
end


function info = empty_density_repair_info()
info = struct();
info.enabled = false;
info.requested = false;
info.strength = 0.0;
info.Kvir = 0.0;
info.smoothPasses = 0;
info.minCellCount = 1.0;
info.maxDisplacementFraction = 0.0;
info.restoreVelocity = true;
info.method = '';
info.meanNBefore = NaN;
info.stdNBefore = NaN;
info.minNBefore = NaN;
info.maxNBefore = NaN;
info.outBand20Before = NaN;
info.meanNAfter = NaN;
info.stdNAfter = NaN;
info.minNAfter = NaN;
info.maxNAfter = NaN;
info.outBand20After = NaN;
info.deltaStdN = NaN;
info.deltaOutBand20 = NaN;
info.pVirMean = 0.0;
info.pVirRms = 0.0;
info.pVirMaxAbs = 0.0;
info.gradPvirRms = 0.0;
info.gradPvirMaxAbs = 0.0;
info.gridDisplacementRms = 0.0;
info.gridDisplacementMaxAbs = 0.0;
info.particleDisplacementRms = 0.0;
info.particleDisplacementMaxAbs = 0.0;
info.particleDisplacementMeanX = 0.0;
info.particleDisplacementMeanY = 0.0;
info.nLimitedParticles = 0;
info.nWallClamped = 0;
info.velocityRestored = true;
info.velocityRestoreDeltaRms = 0.0;
info.velocityRestoreDeltaMaxAbs = 0.0;
info.velocityRestoreResidualBeforeRms = NaN;
info.velocityRestoreResidualAfterRms = NaN;
info.NBefore = [];
info.NAfter = [];
info.Pvir = [];
info.gradPvirX = [];
info.gradPvirY = [];
info.dXGrid = [];
info.dYGrid = [];
info.displacement = [];
info.dvRestore = [];
info.Gtarget = [];
info.GafterMoveRaw = [];
info.GafterRestore = [];
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
