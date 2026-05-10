function result = test_projection_periodic_x_neumann_y_grid_only()
%TEST_PROJECTION_PERIODIC_X_NEUMANN_Y_GRID_ONLY Validate bounded-y projection.

rng(11);
params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Nx = 24;
params.Ny = 12;
params.dt = 1e-3;
params.rho0 = 1.0;

Ux = randn(params.Nx, params.Ny);
Uy = randn(params.Nx, params.Ny);
proj = projection_project_grid_periodic_x_neumann_y(Ux, Uy, params);

result = struct();
result.rmsDivBefore = proj.rmsDivBefore;
result.rmsDivAfter = proj.rmsDivAfter;
result.reduction = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
result.passFinite = all(isfinite(proj.Ux(:))) && all(isfinite(proj.Uy(:))) && all(isfinite(proj.p(:)));
result.passProjected = result.reduction < 1e-10;
result.passed = result.passFinite && result.passProjected;

fprintf('\n=== test_projection_periodic_x_neumann_y_grid_only ===\n');
fprintf('rms div before : %.12e\n', result.rmsDivBefore);
fprintf('rms div after  : %.12e\n', result.rmsDivAfter);
fprintf('reduction      : %.12e\n', result.reduction);
fprintf('passed         : %d\n', result.passed);

if ~result.passed
    error('test_projection_periodic_x_neumann_y_grid_only failed.');
end
end
