function result = test_projection_virial_density_repair()
%TEST_PROJECTION_VIRIAL_DENSITY_REPAIR Unit smoke test for Q7c position repair.

params = struct();
params.Lx = 1.0;
params.Ly = 1.0;
params.Nx = 8;
params.Ny = 6;
params.gamma = 6;
params.dt = 1.0e-3;
params.kBT = 1.0;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;
params.seed = 123;
params.useVirialDensityRepair = true;
params.virialDensityRepairStrength = 0.05;
params.virialDensityRepairK = 1.0;
params.virialDensityRepairSmoothPasses = 1;
params.virialDensityRepairMaxDisplacementFraction = 0.10;
params.virialDensityRepairRestoreVelocity = true;
params.virialDensityRepairInterpolationMethod = 'nearest';
params.projectionInterpolationMethod = 'nearest';

rng(params.seed);
[state, ~] = projection_initialize_particles(params);

% Create a controlled density defect by moving a few particles from a known
% source cell into a different target cell.  Do not use the first particle IDs:
% with exact_per_cell initialization they can already belong to the target cell,
% which would leave N perfectly uniform and make the repair a no-op.
dx = params.Lx / params.Nx;
dy = params.Ly / params.Ny;
ix0 = min(max(floor(state.x(:, 1) / dx) + 1, 1), params.Nx);
iy0 = min(max(floor(state.x(:, 2) / dy) + 1, 1), params.Ny);
sourceIx = max(2, params.Nx);
sourceIy = max(2, params.Ny);
targetIx = 1;
targetIy = 1;
idsSource = find(ix0 == sourceIx & iy0 == sourceIy);
nMove = min(4, numel(idsSource));
if nMove < 1
    error('Could not find source-cell particles for repair test.');
end
ids = idsSource(1:nMove);
state.x(ids, 1) = (targetIx - 0.5) * dx;
state.x(ids, 2) = (targetIy - 0.5) * dy;

Gtarget = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);
[stateOut, info] = projection_apply_virial_density_repair(state, params, ...
    'periodicX', true, 'periodicY', false, 'targetGrid', Gtarget);

sameCount = size(stateOut.x, 1) == size(state.x, 1);
finiteState = all(isfinite(stateOut.x(:))) && all(isfinite(stateOut.v(:)));
finiteInfo = isfinite(info.particleDisplacementRms) && isfinite(info.velocityRestoreDeltaRms);
repairRan = info.enabled && info.particleDisplacementRms > 0;
inDomain = all(stateOut.x(:, 1) >= 0 & stateOut.x(:, 1) < params.Lx) && ...
           all(stateOut.x(:, 2) >= 0 & stateOut.x(:, 2) < params.Ly);
restoreImproves = info.velocityRestoreResidualAfterRms <= info.velocityRestoreResidualBeforeRms + 1e-12;

result = struct();
result.sameCount = sameCount;
result.finiteState = finiteState;
result.finiteInfo = finiteInfo;
result.repairRan = repairRan;
result.inDomain = inDomain;
result.restoreImproves = restoreImproves;
result.particleDisplacementRms = info.particleDisplacementRms;
result.particleDisplacementMaxAbs = info.particleDisplacementMaxAbs;
result.stdBefore = info.stdNBefore;
result.stdAfter = info.stdNAfter;
result.velocityRestoreResidualBeforeRms = info.velocityRestoreResidualBeforeRms;
result.velocityRestoreResidualAfterRms = info.velocityRestoreResidualAfterRms;
result.passed = sameCount && finiteState && finiteInfo && repairRan && inDomain && restoreImproves;

fprintf('\n=== test_projection_virial_density_repair ===\n');
fprintf('dx repair rms/max        : %.6g / %.6g\n', result.particleDisplacementRms, result.particleDisplacementMaxAbs);
fprintf('std before/after repair  : %.6g / %.6g\n', result.stdBefore, result.stdAfter);
fprintf('restore residual b/a     : %.6g / %.6g\n', result.velocityRestoreResidualBeforeRms, result.velocityRestoreResidualAfterRms);
fprintf('pass same count/domain   : %d / %d\n', sameCount, inDomain);
fprintf('pass finite/repair ran   : %d / %d\n', finiteState && finiteInfo, repairRan);
fprintf('pass restore improves    : %d\n', restoreImproves);
fprintf('passed                   : %d\n', result.passed);

if ~result.passed
    error('test_projection_virial_density_repair failed.');
end
end
