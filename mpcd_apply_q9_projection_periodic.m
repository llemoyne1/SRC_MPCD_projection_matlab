function [stateOut, diag] = mpcd_apply_q9_projection_periodic(stateClassic, params)
%MPCD_APPLY_Q9_PROJECTION_PERIODIC Apply Q6/Q9 projection in a fully periodic box.
%
% Q6: massFluxProjectionMode='off'.
% Q9: massFluxProjectionMode='relax_to_uniform_lowk', usually after Q6.
% The optional global momentum correction is applied stage-by-stage to Q6,
% Q9 and cleanup particle kicks.

projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
massFluxProjectionMode = char(string(get_param(params, 'massFluxProjectionMode', 'off')));
massFluxProjectionStrength = get_param(params, 'massFluxProjectionStrength', projectionStrength);
massFluxDensityRelaxationBeta = get_param(params, 'massFluxDensityRelaxationBeta', 0.0);
massFluxApplyAfterVelocityProjection = logical(get_param(params, 'massFluxApplyAfterVelocityProjection', false));
massFluxFinalVelocityProjectionCleanup = logical(get_param(params, 'massFluxFinalVelocityProjectionCleanup', false));
massFluxFinalVelocityProjectionStrength = get_param(params, 'massFluxFinalVelocityProjectionStrength', 1.0);
useMassFluxProjection = ~strcmpi(strrep(massFluxProjectionMode, '-', '_'), 'off');
massFluxProjectionOperator = char(string(get_param(params, 'massFluxProjectionOperator', 'periodic_fft')));
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));
thermostatAfterProjection = logical(get_param(params, 'thermostatAfterProjection', ...
    get_param(params, 'thermostatAfterStep', false))); % common TG post-step thermostat flag
computeDiagnostics = logical(get_param(params, 'computeDiagnostics', true));

validate_state(stateClassic);
momentumCorrectionInfo = init_momentum_correction_info(params);

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
else
    popAfterClassic = [];
    thermalBeforeProjection = [];
end

if useMassFluxProjection
    appliedProjectionStrength = massFluxProjectionStrength;
    if massFluxApplyAfterVelocityProjection
        projectionKind = 'velocity_plus_mass_flux';
        stateForMassFlux = stateClassic;
        if projectionEnable && projectionStrength ~= 0
            dvVelocity = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
                'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
            dvVelocityApplied = projectionStrength * dvVelocity;
        else
            dvVelocityApplied = zeros(size(stateClassic.v));
        end
        [stateForMassFlux.v, momentumCorrectionInfo.q6] = ...
            projection_apply_global_momentum_correction(stateClassic.v, dvVelocityApplied, params, 'stage', 'q6');

        GforMassFlux = projection_deposit_particles_to_grid(stateForMassFlux.x, stateForMassFlux.v, params, ...
            'periodicX', true, 'periodicY', true, 'minCount', 1);
        massFluxProj = project_periodic_mass_flux_with_selected_operator(GforMassFlux.N, GforMassFlux.Ux, GforMassFlux.Uy, params, ...
            massFluxProjectionOperator, massFluxProjectionMode, massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dvMassFlux = projection_interpolate_grid_delta_to_particles(stateForMassFlux.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
            dvMassFluxApplied = massFluxProjectionStrength * dvMassFlux;
        else
            dvMassFluxApplied = zeros(size(stateClassic.v));
        end
        [stateOut.v, momentumCorrectionInfo.q9] = ...
            projection_apply_global_momentum_correction(stateForMassFlux.v, dvMassFluxApplied, params, 'stage', 'q9');
    else
        projectionKind = 'mass_flux';
        massFluxProj = project_periodic_mass_flux_with_selected_operator(Gbefore.N, Gbefore.Ux, Gbefore.Uy, params, ...
            massFluxProjectionOperator, massFluxProjectionMode, massFluxDensityRelaxationBeta);
        if projectionEnable && massFluxProjectionStrength ~= 0
            dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, massFluxProj.dUx, massFluxProj.dUy, params, ...
                'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
            dvMassFluxApplied = massFluxProjectionStrength * dv;
        else
            dvMassFluxApplied = zeros(size(stateClassic.v));
        end
        [stateOut.v, momentumCorrectionInfo.q9] = ...
            projection_apply_global_momentum_correction(stateClassic.v, dvMassFluxApplied, params, 'stage', 'q9');
    end
elseif projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    dvVelocityApplied = projectionStrength * dv;
    [stateOut.v, momentumCorrectionInfo.q6] = ...
        projection_apply_global_momentum_correction(stateClassic.v, dvVelocityApplied, params, 'stage', 'q6');
else
    stateOut = stateClassic;
end

if useMassFluxProjection && massFluxFinalVelocityProjectionCleanup && massFluxFinalVelocityProjectionStrength ~= 0
    GcleanupBefore = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', true, 'minCount', 1);
    projCleanup = projection_project_grid_periodic_fft(GcleanupBefore.Ux, GcleanupBefore.Uy, params);
    dvCleanup = projection_interpolate_grid_delta_to_particles(stateOut.x, projCleanup.dUx, projCleanup.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    dvCleanupApplied = massFluxFinalVelocityProjectionStrength * dvCleanup;
    [stateOut.v, momentumCorrectionInfo.cleanup] = ...
        projection_apply_global_momentum_correction(stateOut.v, dvCleanupApplied, params, 'stage', 'cleanup');
    finalVelocityCleanup = struct();
    finalVelocityCleanup.enabled = true;
    finalVelocityCleanup.strength = massFluxFinalVelocityProjectionStrength;
    finalVelocityCleanup.rmsDivBefore = projCleanup.rmsDivBefore;
    finalVelocityCleanup.rmsDivAfter = projCleanup.rmsDivAfter;
    finalVelocityCleanup.maxAbsDivBefore = projCleanup.maxAbsDivBefore;
    finalVelocityCleanup.maxAbsDivAfter = projCleanup.maxAbsDivAfter;
    finalVelocityCleanup.divReduction = projCleanup.rmsDivAfter / max(projCleanup.rmsDivBefore, eps);
    finalVelocityCleanup.dvRms = sqrt(mean(sum(dvCleanupApplied.^2, 2), 'omitnan'));
    finalVelocityCleanup.proj = projCleanup;
end

momentumCorrectionInfo = finalize_momentum_correction_info(momentumCorrectionInfo);

if thermostatAfterProjection
    [stateOut.v, thermostatInfo] = projection_apply_cell_thermostat(stateOut.x, stateOut.v, params, ...
        'periodicX', true, 'periodicY', true);
    %thermCheck = projection_thermal_diagnostics(x, v, params,     'periodicX', true, 'periodicY', true, 'minCount', 1);
thermCheck = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);

methodLabel = char(string(get_param(params, 'method', 'projection')));

% fprintf(['THERMO %s: target=%.6g, strength=%.3g, cells=%d, ', ...
%          'info before=%.6g, info after=%.6g, diag cellRel=%.6g, diag localMean=%.6g, ', ...
%          'scale mean/min/max=%.6g/%.6g/%.6g\n'], ...
%     upper(methodLabel), ...
%     thermostatInfo.targetKBT, thermostatInfo.strength, ...
%     thermostatInfo.nThermostattedCells, ...
%     thermostatInfo.meanKBTBefore, thermostatInfo.meanKBTAfter, ...
%     thermCheck.kBTCellRelative, thermCheck.localKBTMean, ...
%     thermostatInfo.meanScale, thermostatInfo.minScale, thermostatInfo.maxScaleApplied);
% fprintf(['THERMO classic: target=%.6g, strength=%.3g, cells=%d, ', ...
%          'info before=%.6g, info after=%.6g, diag cellRel=%.6g, diag localMean=%.6g, ', ...
%          'scale mean/min/max=%.6g/%.6g/%.6g\n'], ...
%     thermostatInfo.targetKBT, thermostatInfo.strength, ...
%     thermostatInfo.nThermostattedCells, ...
%     thermostatInfo.meanKBTBefore, thermostatInfo.meanKBTAfter, ...
%     thermCheck.kBTCellRelative, thermCheck.localKBTMean, ...
%     thermostatInfo.meanScale, thermostatInfo.minScale, thermostatInfo.maxScaleApplied);
else
    thermostatInfo = empty_thermostat_info();
end

if ~computeDiagnostics
    diag = minimal_diag(proj, massFluxProj, thermostatInfo, projectionEnable, projectionStrength, projectionKind, ...
        massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, appliedProjectionStrength, ...
        massFluxApplyAfterVelocityProjection, massFluxFinalVelocityProjectionCleanup, massFluxFinalVelocityProjectionStrength, ...
        finalVelocityCleanup, momentumCorrectionInfo, stateOut.v - stateClassic.v);
    diag.massFluxProjectionOperator = canonical_mass_flux_operator(massFluxProjectionOperator);
    return;
end

thermalAfterProjection = projection_thermal_diagnostics(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
popAfterProjection = projection_population_diagnostics(stateOut.x, params, ...
    'periodicX', true, 'periodicY', true);
Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
projAfterParticles = projection_project_grid_periodic_fft(Gafter.Ux, Gafter.Uy, params);
if useMassFluxProjection
    massFluxAfterParticles = project_periodic_mass_flux_with_selected_operator(Gafter.N, Gafter.Ux, Gafter.Uy, params, ...
        massFluxProjectionOperator, massFluxProjectionMode, massFluxDensityRelaxationBeta);
else
    massFluxAfterParticles = empty_mass_flux_projection();
end

dvTotal = stateOut.v - stateClassic.v;
diag = minimal_diag(proj, massFluxProj, thermostatInfo, projectionEnable, projectionStrength, projectionKind, ...
    massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, appliedProjectionStrength, ...
    massFluxApplyAfterVelocityProjection, massFluxFinalVelocityProjectionCleanup, massFluxFinalVelocityProjectionStrength, ...
    finalVelocityCleanup, momentumCorrectionInfo, dvTotal);
diag.massFluxProjectionOperator = canonical_mass_flux_operator(massFluxProjectionOperator);
diag.computeDiagnostics = true;
diag.rmsDivParticleAfter = projAfterParticles.rmsDivBefore;
diag.maxAbsDivParticleAfter = projAfterParticles.maxAbsDivBefore;
diag.divReductionParticle = diag.rmsDivParticleAfter / max(proj.rmsDivBefore, eps);
diag.rmsMassFluxDivParticleAfter = massFluxAfterParticles.rmsDivMassBefore;
diag.kBTCellBeforeProjection = local_get_kBT_cell(thermalBeforeProjection);
diag.kBTCellAfterProjection = local_get_kBT_cell(thermalAfterProjection);
diag.thermalBeforeProjection = thermalBeforeProjection;
diag.thermalAfterProjection = thermalAfterProjection;
diag.popAfterClassic = popAfterClassic;
diag.popAfterProjection = popAfterProjection;
diag.NStdAfterProjection = popAfterProjection.stdN;
diag.NOutBandAfterProjection = popAfterProjection.outBandFraction;
diag.Gbefore = Gbefore;
diag.Gafter = Gafter;
diag.projAfterParticles = projAfterParticles;
diag.massFluxAfterParticles = massFluxAfterParticles;
diag.momBefore = mean(stateClassic.v, 1, 'omitnan');
diag.momAfter = mean(stateOut.v, 1, 'omitnan');
end


function proj = project_periodic_mass_flux_with_selected_operator(N, Ux, Uy, params, operatorName, mode, beta)
%PROJECT_PERIODIC_MASS_FLUX_WITH_SELECTED_OPERATOR Dispatch periodic Q9 mass-flux projection.
%
% Historical Taylor--Green runs use the FFT periodic projector.  This
% dispatcher also allows the finite-volume/general-BC elliptic projector to
% be exercised on a fully periodic domain, which is useful for checking the
% consistency of the future unified operator path.

op = canonical_mass_flux_operator(operatorName);
switch op
    case 'periodic_fft'
        proj = projection_project_mass_flux_periodic_fft(N, Ux, Uy, params, ...
            'mode', mode, 'relaxationBeta', beta);
        proj.operator = op;
    case 'general_bc_periodic'
        p = params;
        % Force a fully periodic descriptor, independent of the Poiseuille
        % defaults used inside projection_project_mass_flux_general_bc.
        p.massFluxBC = struct();
        p.massFluxBC.left = 'periodic';
        p.massFluxBC.right = 'periodic';
        p.massFluxBC.bottom = 'periodic';
        p.massFluxBC.top = 'periodic';
        p.massFluxBCLeft = 'periodic';
        p.massFluxBCRight = 'periodic';
        p.massFluxBCBottom = 'periodic';
        p.massFluxBCTop = 'periodic';
        p.boundary_left = 'periodic';
        p.boundary_right = 'periodic';
        p.boundary_bottom = 'periodic';
        p.boundary_top = 'periodic';
        if ~isfield(p, 'massFluxTargetFilter') || isempty(p.massFluxTargetFilter)
            p.massFluxTargetFilter = 'elliptic_lowpass';
        end
        if ~isfield(p, 'massFluxEllipticAlphaMode') || isempty(p.massFluxEllipticAlphaMode)
            p.massFluxEllipticAlphaMode = 'constant';
        end
        if ~isfield(p, 'massFluxEllipticUseSolveCache') || isempty(p.massFluxEllipticUseSolveCache)
            p.massFluxEllipticUseSolveCache = true;
        end
        proj = projection_project_mass_flux_general_bc(N, Ux, Uy, p, ...
            'mode', mode, 'relaxationBeta', beta);
        proj.operator = op;
    otherwise
        error('Unsupported massFluxProjectionOperator: %s', operatorName);
end
end

function op = canonical_mass_flux_operator(operatorName)
op = lower(strrep(char(string(operatorName)), '-', '_'));
switch op
    case {'periodic_fft','fft','fft_periodic','spectral','spectral_fft'}
        op = 'periodic_fft';
    case {'general_bc','general','elliptic','finite_volume','fv','general_bc_periodic','periodic_general_bc'}
        op = 'general_bc_periodic';
    otherwise
        error('Unknown massFluxProjectionOperator: %s', operatorName);
end
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
totalRaw = [0 0]; totalCorr = [0 0]; totalRes = [0 0]; nActive = 0;
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

function diag = minimal_diag(proj, massFluxProj, thermostatInfo, projectionEnable, projectionStrength, projectionKind, massFluxProjectionMode, massFluxProjectionStrength, massFluxDensityRelaxationBeta, appliedProjectionStrength, massFluxApplyAfterVelocityProjection, massFluxFinalVelocityProjectionCleanup, massFluxFinalVelocityProjectionStrength, finalVelocityCleanup, momentumCorrectionInfo, dv)
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
diag.computeDiagnostics = false;
diag.rmsDivBefore = proj.rmsDivBefore;
diag.rmsDivProjectedAfter = proj.rmsDivAfter;
diag.maxAbsDivBefore = proj.maxAbsDivBefore;
diag.maxAbsDivProjectedAfter = proj.maxAbsDivAfter;
diag.meanDivBefore = proj.meanDivBefore;
diag.divReductionProjected = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
diag.rmsDivParticleAfter = NaN;
diag.maxAbsDivParticleAfter = NaN;
diag.rmsMassFluxDivBefore = massFluxProj.rmsDivMassBefore;
diag.rmsMassFluxDivTarget = massFluxProj.rmsTargetDivMass;
diag.rmsMassFluxDivProjectedAfter = massFluxProj.rmsDivMassAfter;
diag.rmsMassFluxDivResidual = massFluxProj.rmsDivMassResidual;
diag.rmsMassFluxDivParticleAfter = NaN;
diag.massFluxDivReduction = massFluxProj.divMassReduction;
diag.massFluxProj = massFluxProj;
diag.dvRms = sqrt(mean(sum(dv.^2, 2), 'omitnan'));
diag.thermostatAfterProjection = thermostatInfo.enabled;
diag.thermostatInfo = thermostatInfo;
diag.thermostatRmsVelocityChange = thermostatInfo.rmsVelocityChange;
diag.thermostatMeanScale = thermostatInfo.meanScale;
diag.thermostatNCells = thermostatInfo.nThermostattedCells;
diag.kBTCellAfterProjection = NaN;
end

function info = empty_velocity_cleanup()
info = struct('enabled', false, 'strength', 0.0, 'rmsDivBefore', NaN, 'rmsDivAfter', NaN, ...
    'maxAbsDivBefore', NaN, 'maxAbsDivAfter', NaN, 'divReduction', NaN, 'dvRms', NaN, 'proj', []);
end

function info = empty_mass_flux_projection()
info = struct();
info.mode = 'off'; info.beta = 0.0;
info.rmsDivMassBefore = NaN; info.rmsTargetDivMass = NaN;
info.rmsDivMassAfter = NaN; info.rmsDivMassResidual = NaN;
info.maxAbsDivMassBefore = NaN; info.maxAbsDivMassAfter = NaN; info.maxAbsDivMassResidual = NaN;
info.divMassReduction = NaN;
end

function info = empty_thermostat_info()
info = struct('enabled', false, 'targetKBT', 0.0, 'strength', 0.0, ...
    'minParticlesPerCell', 0, 'maxScale', 1.0, 'nThermostattedCells', 0, ...
    'meanScale', 1.0, 'minScale', 1.0, 'maxScaleApplied', 1.0, ...
    'meanKBTBefore', 0.0, 'meanKBTAfter', 0.0, 'rmsVelocityChange', 0.0);
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v')
    error('state must be a struct with fields x and v.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end


function val = local_get_kBT_cell(therm)
%LOCAL_GET_KBT_CELL Compatibility helper for thermal diagnostics.
% Newer semi-parallel scripts use kBTCellRelative; older snapshots used
% kBTCellMean.
if isfield(therm,'kBTCellRelative')
    val = therm.kBTCellRelative;
elseif isfield(therm,'kBTCellMean')
    val = therm.kBTCellMean;
elseif isfield(therm,'thermalKineticEnergy')
    val = therm.thermalKineticEnergy;
elseif isfield(therm,'localKBTMean')
    val = therm.localKBTMean;
else
    val = NaN;
end
end
