function [stateOut, diag] = mpcd_step_projection_periodic(state, params)
%MPCD_STEP_PROJECTION_PERIODIC Classic MPCD step plus periodic pressure projection.
%
%   [stateOut, diag] = mpcd_step_projection_periodic(state, params)
%
% This standalone prototype couples particles to the grid-only FFT projection:
%
%   1. classic periodic SRD/MPCD step;
%   2. deposit U*(x,y) from particles to cell centers;
%   3. project U* so that div(U) = 0 on a periodic grid;
%   4. apply dU = U - U* back to particles.
%
% Required params fields are the same as mpcd_step_classic_periodic.
%
% Optional params fields:
%   projectionEnable              default true
%   projectionStrength            default 1
%   projectionInterpolationMethod default 'nearest'
%
% 'nearest' is intentionally the default for this first prototype because it
% makes the deposited grid velocity after correction exactly consistent with
% the projected cell velocity when cells are occupied. 'bilinear' is smoother
% but no longer makes the particle re-deposit an exact grid projection.

projectionEnable = logical(get_param(params, 'projectionEnable', true));
projectionStrength = get_param(params, 'projectionStrength', 1.0);
projectionInterpolationMethod = char(string(get_param(params, 'projectionInterpolationMethod', 'nearest')));

[stateClassic, classicDiag] = mpcd_step_classic_periodic(state, params);

Gbefore = projection_deposit_particles_to_grid(stateClassic.x, stateClassic.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);

proj = projection_project_grid_periodic_fft(Gbefore.Ux, Gbefore.Uy, params);

stateOut = stateClassic;
if projectionEnable && projectionStrength ~= 0
    dv = projection_interpolate_grid_delta_to_particles(stateClassic.x, proj.dUx, proj.dUy, params, ...
        'periodicX', true, 'periodicY', true, 'method', projectionInterpolationMethod);
    stateOut.v = stateClassic.v + projectionStrength * dv;
else
    dv = zeros(size(stateClassic.v));
end

Gafter = projection_deposit_particles_to_grid(stateOut.x, stateOut.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
projAfterParticles = projection_project_grid_periodic_fft(Gafter.Ux, Gafter.Uy, params);

momBefore = mean(stateClassic.v, 1);
momAfter = mean(stateOut.v, 1);

diag = struct();
diag.classic = classicDiag;
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
