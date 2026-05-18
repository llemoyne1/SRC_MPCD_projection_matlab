function [stateOut, diag] = mpcd_apply_q9_projection_channel(stateClassic, params)
%MPCD_APPLY_Q9_PROJECTION_CHANNEL Apply Q6/Q9 projection to a channel state.
%
%   [stateOut, diag] = mpcd_apply_q9_projection_channel(stateClassic, params)
%
% Generic post-SRC projection block for x-periodic / bounded-y domains.
% It does not perform the classic MPCD step itself.  This allows the same
% projection block to be used after Poiseuille, piston, and later obstacle
% kernels.
%
% Q6 mode:
%   massFluxProjectionMode = 'off'
%
% Q9 mode:
%   projectionStrength = 1
%   massFluxProjectionMode = 'relax_to_uniform_lowk'
%   massFluxProjectionStrength = 1
%   massFluxDensityRelaxationBeta = 0.002
%   massFluxApplyAfterVelocityProjection = true
%   massFluxTargetFilter = 'lowpass_fft'
%   massFluxLowKMaxIndex = 2

projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
massFluxProjectionMode = char(string(get_param(params, 'massFluxProjectionMode', 'off')));
massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', projectionStrength);
massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 0.0);
massFluxApplyAfterVelocityProjection = logical(get_param(params, 'massFluxApplyAfterVelocityProjection', false));
% Optional diagnostic/experimental cleanup: the Q9 mass-flux correction is
% applied after the Q6 velocity projection and can reintroduce a small
% velocity divergence.  Keep this disabled for the current Q9 reference, but
% allow a Q9-clean variant for piston/cylinder sensitivity tests.
massFluxFinalVelocityProjectionCleanup = logical(get_param(params, 'massFluxFinalVelocityProjectionCleanup', false));
massFluxFinalVelocityProjectionStrength = get_param(params, 'massFluxFinalVelocityProjectionStrength', 1.0);
useMassFluxProjection = ~strcmpi(strrep(massFluxProjectionMode, '-', '_'), 'off');
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
projectionTransportDiagnosticsEnable = logical(get_param(params, 'projectionTransportDiagnosticsEnable', true));
densityTransportDiagnosticsEnable = logical(get_param(params, 'densityTransportDiagnosticsEnable', true));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', false));
computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));

% Exact global momentum correction for numerical projection stages.  This
% removes the spatially uniform component of the particle velocity increment
% introduced by Q6/Q9/cleanup without touching physical impulses from the
% classic SRC step, walls, body forces, or the mean-flow controller applied
% by the runner.
momentumCorrectionInfo = init_momentum_correction_info(params);

validate_state(stateClassic);

if computeDiagnostics
    popAfterClassic = projection_population_diagnostics(stateClassic.x, params, ...
        'periodicX', true, 'periodicY', false);
end

Gbefore = projection_deposit_particles_to_grid(stateClassic.x, stateClassic.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
Gbefore = apply_solid_mask_to_grid(Gbefore, params);
if computeDiagnostics
    thermalBeforeProjection = projection_thermal_diagnostics(stateClassic.x, stateClassic.v, params, ...
        'periodicX', true, 'periodicY', false, 'minCount', 1);
end

proj = projection_project_grid_periodic_x_neumann_y(Gbefore.Ux, Gbefore.Uy, params);
massFluxProj = empty_mass_flux_projection();
projectionKind = 'velocity';
appliedProjectionStrength = projectionStrength;

stateOut = stateClassic;
finalVelocityCleanup = empty_velocity_cleanup();
if useMassFluxProjection
    appliedProjectionStrength = massFluxProjectionStrength;

    if massFluxApplyAfterVelocityProjection
        projectionKind = 'velocity_plus_mass_flux';
        stateForMassFlux = stateClassic;
        if projectionEnable && projectionStrength ~= 0
            dvVelocity = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
                'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
            dvVelocityApplied = projectionStrength * dvVelocity;
        else
            dvVelocity = zeros(size(stateClassic.v)); %#ok<NASGU>
            dvVelocityApplied = zeros(size(stateClassic.v));
        end
        [stateForMassFlux.v, momentumCorrectionInfo.q6] = ...
            projection_apply_global_momentum_correction(stateClassic.v, dvVelocityApplied, params, 'stage', 'q6');
        GforMassFlux = projection_deposit_particles_to_grid(stateForMassFlux.x, stateForMassFlux.v, params, ...
            'periodicX', true, 'periodicY', false, 'minCount', 1);
        GforMassFlux = apply_solid_mask_to_grid(GforMassFlux, params);
        massFluxProj = projection_project_mass_flux_channel_operator(GforMassFlux.N, GforMassFlux.Ux, GforMassFlux.Uy, params, ...
            'mode', massFluxProjectionMode, ...
            'relaxationBeta', massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dvMassFlux = projection_interpolate_grid_delta_to_particles(stateForMassFlux.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
            dvMassFluxApplied = massFluxProjectionStrength * dvMassFlux;
        else
            dvMassFlux = zeros(size(stateClassic.v)); %#ok<NASGU>
            dvMassFluxApplied = zeros(size(stateClassic.v));
        end
        [stateOut.v, momentumCorrectionInfo.q9] = ...
            projection_apply_global_momentum_correction(stateForMassFlux.v, dvMassFluxApplied, params, 'stage', 'q9');
        dv = stateOut.v - stateClassic.v;
    else
        projectionKind = 'mass_flux';
        massFluxProj = projection_project_mass_flux_channel_operator(Gbefore.N, Gbefore.Ux, Gbefore.Uy, params, ...
            'mode', massFluxProjectionMode, ...
            'relaxationBeta', massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
            dvMassFluxApplied = massFluxProjectionStrength * dv;
        else
            dv = zeros(size(stateClassic.v));
            dvMassFluxApplied = zeros(size(stateClassic.v));
        end
        [stateOut.v, momentumCorrectionInfo.q9] = ...
            projection_apply_global_momentum_correction(stateClassic.v, dvMassFluxApplied, params, 'stage', 'q9');
        dv = stateOut.v - stateClassic.v;
    end
elseif projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
    dvVelocityApplied = projectionStrength * dv;
    [stateOut.v, momentumCorrectionInfo.q6] = ...
        projection_apply_global_momentum_correction(stateClassic.v, dvVelocityApplied, params, 'stage', 'q6');
    dv = stateOut.v - stateClassic.v;
else
    dv = zeros(size(stateClassic.v));
end

if useMassFluxProjection && massFluxFinalVelocityProjectionCleanup && massFluxFinalVelocityProjectionStrength ~= 0
    GcleanupBefore = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', false, 'minCount', 1);
    GcleanupBefore = apply_solid_mask_to_grid(GcleanupBefore, params);
    projCleanup = projection_project_grid_periodic_x_neumann_y(GcleanupBefore.Ux, GcleanupBefore.Uy, params);
    dvCleanup = projection_interpolate_grid_delta_to_particles(stateOut.x, projCleanup.dUx, projCleanup.dUy, params, ...
        'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
    dvCleanupApplied = massFluxFinalVelocityProjectionStrength * dvCleanup;
    [stateOut.v, momentumCorrectionInfo.cleanup] = ...
        projection_apply_global_momentum_correction(stateOut.v, dvCleanupApplied, params, 'stage', 'cleanup');
    dv = stateOut.v - stateClassic.v;
    finalVelocityCleanup = struct();
    finalVelocityCleanup.enabled = true;
    finalVelocityCleanup.strength = massFluxFinalVelocityProjectionStrength;
    finalVelocityCleanup.rmsDivBefore = projCleanup.rmsDivBefore;
    finalVelocityCleanup.rmsDivAfter = projCleanup.rmsDivAfter;
    finalVelocityCleanup.maxAbsDivBefore = projCleanup.maxAbsDivBefore;
    finalVelocityCleanup.maxAbsDivAfter = projCleanup.maxAbsDivAfter;
    finalVelocityCleanup.divReduction = projCleanup.rmsDivAfter / max(projCleanup.rmsDivBefore, eps);
    finalVelocityCleanup.dvRms = sqrt(mean(sum((massFluxFinalVelocityProjectionStrength * dvCleanup).^2, 2)));
    finalVelocityCleanup.proj = projCleanup;
end

momentumCorrectionInfo = finalize_momentum_correction_info(momentumCorrectionInfo);

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
    diag = minimal_projection_diag(proj, massFluxProj, thermostatInfo, projectionEnable, ...
        projectionStrength, projectionInterpolationMethod, dv, projectionKind, ...
        massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, ...
        appliedProjectionStrength, massFluxApplyAfterVelocityProjection, ...
        massFluxFinalVelocityProjectionCleanup, massFluxFinalVelocityProjectionStrength, finalVelocityCleanup, ...
        momentumCorrectionInfo);
    return;
end

thermalAfterProjection = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);

popAfterProjection = projection_population_diagnostics(stateOut.x, params, ...
    'periodicX', true, 'periodicY', false);

Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
Gafter = apply_solid_mask_to_grid(Gafter, params);
projAfterParticles = projection_project_grid_periodic_x_neumann_y(Gafter.Ux, Gafter.Uy, params);
if useMassFluxProjection
    massFluxAfterParticles = projection_project_mass_flux_channel_operator(Gafter.N, Gafter.Ux, Gafter.Uy, params, ...
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
diag.projectionMomentumCorrectionEnable = momentumCorrectionInfo.enabled;
diag.projectionMomentumCorrectionMode = momentumCorrectionInfo.mode;
diag.projectionMomentumCorrection = momentumCorrectionInfo;
diag.projectionMomentumQ6RawMeanDVNorm = momentumCorrectionInfo.q6.rawDeltaMeanVNorm;
diag.projectionMomentumQ9RawMeanDVNorm = momentumCorrectionInfo.q9.rawDeltaMeanVNorm;
diag.projectionMomentumCleanupRawMeanDVNorm = momentumCorrectionInfo.cleanup.rawDeltaMeanVNorm;
diag.projectionMomentumTotalRawMeanDVNorm = momentumCorrectionInfo.totalRawMeanVNorm;
diag.projectionMomentumTotalResidualMeanDVNorm = momentumCorrectionInfo.totalResidualMeanVNorm;
diag.finalVelocityCleanup = finalVelocityCleanup;
diag.finalVelocityCleanupRmsDivBefore = finalVelocityCleanup.rmsDivBefore;
diag.finalVelocityCleanupRmsDivAfter = finalVelocityCleanup.rmsDivAfter;
diag.finalVelocityCleanupDivReduction = finalVelocityCleanup.divReduction;
diag.finalVelocityCleanupDvRms = finalVelocityCleanup.dvRms;
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
diag.populationAfterClassic = popAfterClassic;
diag.populationAfterProjection = popAfterProjection;
diag.populationProjectionDeltaMaxAbs = max(abs(popDeltaProjection(:)));
diag.populationProjectionDeltaRms = sqrt(mean(popDeltaProjection(:).^2));
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


function G = apply_solid_mask_to_grid(G, params)
%APPLY_SOLID_MASK_TO_GRID Fill solid cells before rectangular Q6/Q9 projection.
%
% The current projection operators are rectangular channel operators.  For
% obstacle tests, empty cells inside the cylinder would otherwise be
% interpreted as a large density deficit and dominate the low-k correction.
% This helper keeps the obstacle cells neutral for the projection algebra:
%   N = gamma, U = 0, P = 0, valid = true inside the solid mask.
% Particle positions are not changed and diagnostics can still exclude the
% same mask explicitly.
mask = [];
if isfield(params, 'solidMaskForProjection') && ~isempty(params.solidMaskForProjection)
    mask = logical(params.solidMaskForProjection);
elseif isfield(params, 'solidMask') && ~isempty(params.solidMask)
    mask = logical(params.solidMask);
end
if isempty(mask)
    return;
end
if ~isequal(size(mask), size(G.N))
    error('solidMaskForProjection size must match the grid size [%d %d].', size(G.N,1), size(G.N,2));
end
gamma = get_param(params, 'gamma', mean(double(G.N(~mask)), 'omitnan'));
G.N(mask) = gamma;
G.rho(mask) = gamma / max(G.dx * G.dy, eps);
G.Px(mask) = 0;
G.Py(mask) = 0;
G.Ux(mask) = 0;
G.Uy(mask) = 0;
G.valid(mask) = true;
end


function info = init_momentum_correction_info(params)
enabled = logical(get_param(params, 'projectionMomentumCorrectionEnable', true));
mode = char(string(get_param(params, 'projectionMomentumCorrectionMode', 'particle_global_exact')));
info = struct();
info.enabled = enabled;
info.mode = mode;
info.q6 = empty_momentum_stage('q6', enabled, mode);
info.q9 = empty_momentum_stage('q9', enabled, mode);
info.cleanup = empty_momentum_stage('cleanup', enabled, mode);
info.totalRawDeltaP = [0 0];
info.totalRawDeltaPNorm = 0;
info.totalRawMeanVNorm = 0;
info.totalCorrectionDeltaP = [0 0];
info.totalCorrectionDeltaPNorm = 0;
info.totalResidualDeltaP = [0 0];
info.totalResidualDeltaPNorm = 0;
info.totalResidualMeanVNorm = 0;
info.nActiveParticles = 0;
end

function info = finalize_momentum_correction_info(info)
stages = {info.q6, info.q9, info.cleanup};
totalRaw = [0 0];
totalCorr = [0 0];
totalRes = [0 0];
nActive = 0;
for k = 1:numel(stages)
    st = stages{k};
    if isfield(st, 'nActiveParticles') && st.nActiveParticles > 0
        nActive = max(nActive, st.nActiveParticles);
    end
    if isfield(st, 'rawDeltaP'), totalRaw = totalRaw + st.rawDeltaP; end
    if isfield(st, 'correctionDeltaP'), totalCorr = totalCorr + st.correctionDeltaP; end
    if isfield(st, 'residualDeltaP'), totalRes = totalRes + st.residualDeltaP; end
end
info.totalRawDeltaP = totalRaw;
info.totalRawDeltaPNorm = norm(totalRaw);
info.totalCorrectionDeltaP = totalCorr;
info.totalCorrectionDeltaPNorm = norm(totalCorr);
info.totalResidualDeltaP = totalRes;
info.totalResidualDeltaPNorm = norm(totalRes);
info.nActiveParticles = nActive;
if nActive > 0
    info.totalRawMeanVNorm = norm(totalRaw ./ nActive);
    info.totalResidualMeanVNorm = norm(totalRes ./ nActive);
else
    info.totalRawMeanVNorm = NaN;
    info.totalResidualMeanVNorm = NaN;
end
end

function info = empty_momentum_stage(stage, enabled, mode)
info = struct();
info.stage = stage;
info.enabled = enabled;
info.applied = false;
info.mode = mode;
info.nParticles = 0;
info.nActiveParticles = 0;
info.rawDeltaP = [0 0];
info.rawDeltaPNorm = 0;
info.rawDeltaMeanV = [0 0];
info.rawDeltaMeanVNorm = 0;
info.correctionVector = [0 0];
info.correctionVectorNorm = 0;
info.correctionDeltaP = [0 0];
info.correctionDeltaPNorm = 0;
info.residualDeltaP = [0 0];
info.residualDeltaPNorm = 0;
info.residualMeanV = [0 0];
info.residualMeanVNorm = 0;
info.dvRawRms = 0;
info.dvCorrectedRms = 0;
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v')
    error('state must be a struct with fields x and v.');
end
if size(state.x, 2) ~= 2 || size(state.v, 2) ~= 2 || size(state.x, 1) ~= size(state.v, 1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
end

function diag = minimal_projection_diag(proj, massFluxProj, thermostatInfo, projectionEnable, projectionStrength, projectionInterpolationMethod, dv, projectionKind, massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, appliedProjectionStrength, massFluxApplyAfterVelocityProjection, massFluxFinalVelocityProjectionCleanup, massFluxFinalVelocityProjectionStrength, finalVelocityCleanup, momentumCorrectionInfo)
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
diag.projectionMomentumCorrectionEnable = momentumCorrectionInfo.enabled;
diag.projectionMomentumCorrectionMode = momentumCorrectionInfo.mode;
diag.projectionMomentumCorrection = momentumCorrectionInfo;
diag.projectionMomentumQ6RawMeanDVNorm = momentumCorrectionInfo.q6.rawDeltaMeanVNorm;
diag.projectionMomentumQ9RawMeanDVNorm = momentumCorrectionInfo.q9.rawDeltaMeanVNorm;
diag.projectionMomentumCleanupRawMeanDVNorm = momentumCorrectionInfo.cleanup.rawDeltaMeanVNorm;
diag.projectionMomentumTotalRawMeanDVNorm = momentumCorrectionInfo.totalRawMeanVNorm;
diag.projectionMomentumTotalResidualMeanDVNorm = momentumCorrectionInfo.totalResidualMeanVNorm;
diag.finalVelocityCleanup = finalVelocityCleanup;
diag.finalVelocityCleanupRmsDivBefore = finalVelocityCleanup.rmsDivBefore;
diag.finalVelocityCleanupRmsDivAfter = finalVelocityCleanup.rmsDivAfter;
diag.finalVelocityCleanupDivReduction = finalVelocityCleanup.divReduction;
diag.finalVelocityCleanupDvRms = finalVelocityCleanup.dvRms;
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
