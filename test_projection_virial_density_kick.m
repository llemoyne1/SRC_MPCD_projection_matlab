function result = test_projection_virial_density_kick()
%TEST_PROJECTION_VIRIAL_DENSITY_KICK Unit smoke test for the local virial kick.

params = struct();
params.Lx = 1.0;
params.Ly = 1.0;
params.Nx = 4;
params.Ny = 3;
params.gamma = 4;
params.dt = 1.0e-3;
params.virialK = 1.0;
params.virialDensityKickStrength = 0.1;
params.useVirialDensityKick = true;
params.virialSmoothPasses = 0;
params.virialMinCellCount = 1;
params.virialMaxParticleKick = Inf;
params.virialInterpolationMethod = 'nearest';

rng(7);
% Deliberately nonuniform positions: left half is denser than right half.
Np = params.gamma * params.Nx * params.Ny;
x = zeros(Np, 2);
x(:, 1) = params.Lx * rand(Np, 1);
x(:, 2) = params.Ly * rand(Np, 1);
x(1:round(0.7*Np), 1) = 0.25 * params.Lx * rand(round(0.7*Np), 1);
v = zeros(Np, 2);

[vOut, info] = projection_apply_virial_density_kick(x, v, params, ...
    'periodicX', true, 'periodicY', false, 'method', 'nearest');

dv = vOut - v;
passEnabled = info.enabled;
passFinite = all(isfinite(dv(:))) && isfinite(info.particleDVRms) && isfinite(info.pVirRms);
passNonzero = info.particleDVRms > 0 && info.pVirRms > 0;
passNoPositionChange = isequal(size(vOut), size(v));

paramsOff = params;
paramsOff.useVirialDensityKick = false;
[vOff, infoOff] = projection_apply_virial_density_kick(x, v, paramsOff, ...
    'periodicX', true, 'periodicY', false, 'method', 'nearest');
passOffIdentity = max(abs(vOff(:) - v(:))) == 0 && ~infoOff.enabled;

passed = passEnabled && passFinite && passNonzero && passNoPositionChange && passOffIdentity;

result = struct();
result.info = info;
result.passEnabled = passEnabled;
result.passFinite = passFinite;
result.passNonzero = passNonzero;
result.passNoPositionChange = passNoPositionChange;
result.passOffIdentity = passOffIdentity;
result.passed = passed;

fprintf('\n=== test_projection_virial_density_kick ===\n');
fprintf('particle dv rms/max : %.6g / %.6g\n', info.particleDVRms, info.particleDVMaxAbs);
fprintf('Pvir rms/max        : %.6g / %.6g\n', info.pVirRms, info.pVirMaxAbs);
fprintf('passEnabled         : %d\n', passEnabled);
fprintf('passFinite          : %d\n', passFinite);
fprintf('passNonzero         : %d\n', passNonzero);
fprintf('passOffIdentity     : %d\n', passOffIdentity);
fprintf('passed              : %d\n', passed);

if ~passed
    error('test_projection_virial_density_kick failed.');
end
end
