function result = test_projection_mass_flux_grid_only()
%TEST_PROJECTION_MASS_FLUX_GRID_ONLY Unit test for Q8 mass-flux projection.

rng(8101);
params = struct();
params.Lx = 2.0;
params.Ly = 1.0;
params.Nx = 16;
params.Ny = 8;
params.dt = 1.0e-3;
params.gamma = 20;
params.massFluxProjectionRegularization = 1e-12;
params.massFluxMinCellCount = 1.0;

[X, Y] = ndgrid(1:params.Nx, 1:params.Ny); %#ok<ASGLU>

% The relax_to_uniform Neumann/periodic problem is solvable only if the
% prescribed target mass-flux divergence has zero spatial mean. Use an exact
% zero-mean density perturbation around gamma for this grid-only test.
N = params.gamma + 2.0 * (-1).^X + 1.0 * cos(2*pi*Y/params.Ny);
N = N - mean(N(:)) + params.gamma;
Ux = 0.1 * randn(params.Nx, params.Ny) + 0.05 * sin(2*pi*X/params.Nx);
Uy = 0.1 * randn(params.Nx, params.Ny);

proj0 = projection_project_mass_flux_periodic_x_neumann_y(N, Ux, Uy, params, ...
    'mode', 'conservative');

params.massFluxDensityRelaxationBeta = 0.1;
projR = projection_project_mass_flux_periodic_x_neumann_y(N, Ux, Uy, params, ...
    'mode', 'relax_to_uniform', 'relaxationBeta', params.massFluxDensityRelaxationBeta);

projL = projection_project_mass_flux_periodic_x_neumann_y(N, Ux, Uy, params, ...
    'mode', 'relax_to_uniform_lowk', ...
    'relaxationBeta', params.massFluxDensityRelaxationBeta, ...
    'lowKMaxIndex', 2);

result = struct();
result.rmsBefore = proj0.rmsDivMassBefore;
result.rmsAfterConservative = proj0.rmsDivMassAfter;
result.rmsResidualConservative = proj0.rmsDivMassResidual;
result.rmsTargetRelaxRequested = projR.rmsTargetDivMassRequested;
result.rmsTargetRelax = projR.rmsTargetDivMass;
result.rmsTargetProjectionResidual = projR.rmsTargetProjectionResidual;
result.rmsResidualRelax = projR.rmsDivMassResidual;
result.rmsTargetLowK = projL.rmsTargetDivMass;
result.rmsResidualLowK = projL.rmsDivMassResidual;
result.passConservative = proj0.rmsDivMassAfter < 1.0e-9 * max(proj0.rmsDivMassBefore, 1);
result.relResidualRelax = projR.rmsDivMassResidual / max([projR.rmsDivMassBefore, projR.rmsTargetDivMass, 1]);
% The relax_to_uniform target can be large because it contains beta/dt.
% The regularized range projection leaves a small relative residual; for
% this grid-only test we require an accurate relative solve, not machine
% precision absolute cancellation.
result.passRelax = result.relResidualRelax < 1.0e-5;
result.relResidualLowK = projL.rmsDivMassResidual / max([projL.rmsDivMassBefore, projL.rmsTargetDivMass, 1]);
% The low-k target is weaker than the full relax target, so the same
% absolute regularization leaves a slightly larger relative residual.
% This still corresponds to a clean solve of the projected target.
result.passLowK = result.relResidualLowK < 3.0e-5;
result.passFinite = all(isfinite([proj0.Ux(:); proj0.Uy(:); projR.Ux(:); projR.Uy(:); projL.Ux(:); projL.Uy(:)]));
result.passed = result.passConservative && result.passRelax && result.passLowK && result.passFinite;

fprintf('\n=== test_projection_mass_flux_grid_only ===\n');
fprintf('rms div mass before       : %.12e\n', result.rmsBefore);
fprintf('rms div mass after zero   : %.12e\n', result.rmsAfterConservative);
fprintf('rms residual zero         : %.12e\n', result.rmsResidualConservative);
fprintf('rms target relax requested: %.12e\n', result.rmsTargetRelaxRequested);
fprintf('rms target relax used     : %.12e\n', result.rmsTargetRelax);
fprintf('rms target projection res : %.12e\n', result.rmsTargetProjectionResidual);
fprintf('rms residual relax        : %.12e\n', result.rmsResidualRelax);
fprintf('relative residual relax   : %.12e\n', result.relResidualRelax);
fprintf('rms target low-k used     : %.12e\n', result.rmsTargetLowK);
fprintf('rms residual low-k        : %.12e\n', result.rmsResidualLowK);
fprintf('relative residual low-k   : %.12e\n', result.relResidualLowK);
fprintf('pass conservative/relax/lowk : %d / %d / %d\n', result.passConservative, result.passRelax, result.passLowK);
fprintf('passed                    : %d\n', result.passed);

if ~result.passed
    error('test_projection_mass_flux_grid_only failed.');
end
end
