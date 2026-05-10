function result = test_projection_periodic_grid_only()
%TEST_PROJECTION_PERIODIC_GRID_ONLY Validate periodic FFT pressure projection.
%
% Run from MATLAB with:
%   result = test_projection_periodic_grid_only();
%
% The test builds a smooth periodic velocity field with non-zero divergence,
% projects it, and checks that the spectral divergence is reduced by many
% orders of magnitude.

params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Nx = 96;
params.Ny = 48;
params.dt = 1.0e-3;
params.rho0 = 1.0;

x = ((0:params.Nx-1).' + 0.5) * params.Lx / params.Nx;
y = ((0:params.Ny-1)  + 0.5) * params.Ly / params.Ny;
[X, Y] = ndgrid(x, y);

Ux = 0.8 * sin(2*pi*X/params.Lx) .* cos(2*pi*Y/params.Ly) ...
   + 0.3 * sin(4*pi*X/params.Lx);
Uy = 0.6 * cos(2*pi*X/params.Lx) .* sin(2*pi*Y/params.Ly) ...
   - 0.2 * sin(4*pi*Y/params.Ly);

proj = projection_project_grid_periodic_fft(Ux, Uy, params);

reduction = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
passed = reduction < 1.0e-10;

result = struct();
result.passed = passed;
result.rmsDivBefore = proj.rmsDivBefore;
result.rmsDivAfter = proj.rmsDivAfter;
result.reduction = reduction;
result.maxAbsDivBefore = proj.maxAbsDivBefore;
result.maxAbsDivAfter = proj.maxAbsDivAfter;

fprintf('\n=== test_projection_periodic_grid_only ===\n');
fprintf('rms div before : %.12e\n', result.rmsDivBefore);
fprintf('rms div after  : %.12e\n', result.rmsDivAfter);
fprintf('reduction      : %.12e\n', result.reduction);
fprintf('passed         : %d\n\n', result.passed);

if ~passed
    error('Periodic projection test failed: divergence reduction is insufficient.');
end
end
