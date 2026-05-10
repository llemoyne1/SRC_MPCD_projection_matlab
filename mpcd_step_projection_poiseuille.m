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

popBeforeStep = projection_population_diagnostics(state.x, params, ...
    'periodicX', true, 'periodicY', false);

[stateClassic, classicDiag] = mpcd_step_classic_poiseuille(state, params);

popAfterClassic = projection_population_diagnostics(stateClassic.x, params, ...
    'periodicX', true, 'periodicY', false);

Gbefore = projection_deposit_particles_to_grid(stateClassic.x, stateClassic.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);

proj = projection_project_grid_periodic_x_neumann_y(Gbefore.Ux, Gbefore.Uy, params);

stateOut = stateClassic;
if projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', false, 'method', projectionInterpolationMethod);
    stateOut.v = stateClassic.v + projectionStrength * dv;
else
    dv = zeros(size(stateClassic.v));
end

popAfterProjection = projection_population_diagnostics(stateOut.x, params, ...
    'periodicX', true, 'periodicY', false);

Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
projAfterParticles = projection_project_grid_periodic_x_neumann_y(Gafter.Ux, Gafter.Uy, params);

momBefore = mean(stateClassic.v, 1);
momAfter = mean(stateOut.v, 1);

popDeltaProjection = double(popAfterProjection.N) - double(popAfterClassic.N);

diag = struct();
diag.classic = classicDiag;
diag.wallInfo = classicDiag.wallInfo;
diag.projectionEnable = projectionEnable;
diag.projectionStrength = projectionStrength;
diag.projectionInterpolationMethod = projectionInterpolationMethod;
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
