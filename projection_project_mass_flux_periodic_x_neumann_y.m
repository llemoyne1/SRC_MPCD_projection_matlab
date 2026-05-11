function proj = projection_project_mass_flux_periodic_x_neumann_y(N, Ux, Uy, params, varargin)
%PROJECTION_PROJECT_MASS_FLUX_PERIODIC_X_NEUMANN_Y Project N*u on a channel grid.
%
%   proj = projection_project_mass_flux_periodic_x_neumann_y(N, Ux, Uy, params)
%
% Projects the cell-centered mass flux M = N U on a domain periodic in x and
% bounded in y. The correction is applied to the mass flux,
%
%       M_new = M - D' lambda,
%
% with lambda solving
%
%       (D D') lambda = D M - target.
%
% The resulting velocity field is U_new = M_new / max(N, Nmin).
%
% Modes:
%   conservative     target = 0, so div(N U_new) ~= 0.
%   relax_to_uniform target = beta/dt * (N - gamma), so the continuum update
%                    N_new = N - dt div(N U_new) relaxes toward gamma.
%
% This is deliberately a constant-coefficient flux projection. It avoids the
% variable-coefficient operator div((1/N) grad(lambda)), which is more costly
% and less attractive for the future OpenMP/MPI/GPU ports.

if nargin < 4 || isempty(params)
    params = struct();
end
if ~isnumeric(N) || ~isnumeric(Ux) || ~isnumeric(Uy) || ...
        ~isequal(size(N), size(Ux)) || ~isequal(size(Ux), size(Uy))
    error('N, Ux and Uy must be numeric arrays with identical size.');
end

mode = char(string(get_param_default(params, 'massFluxProjectionMode', 'conservative')));
beta = get_param_default(params, 'massFluxDensityRelaxationBeta', 0.0);
regularization = get_param_default(params, 'massFluxProjectionRegularization', 1e-12);
minCellCount = get_param_default(params, 'massFluxMinCellCount', 1.0);
targetFilter = char(string(get_param_default(params, 'massFluxTargetFilter', 'none')));
lowKMaxIndex = get_param_default(params, 'massFluxLowKMaxIndex', get_param_default(params, 'lowKMaxIndex', 2));
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "mode"
            mode = char(string(val));
        case "relaxationbeta"
            beta = val;
        case "regularization"
            regularization = val;
        case "mincellcount"
            minCellCount = val;
        case "targetfilter"
            targetFilter = char(string(val));
        case "lowkmaxindex"
            lowKMaxIndex = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end
modeLower = lower(strrep(mode, '-', '_'));

Lx = get_param_default(params, 'Lx', size(Ux, 1));
Ly = get_param_default(params, 'Ly', size(Ux, 2));
dt = get_param_default(params, 'dt', 1.0);
gamma = get_param_default(params, 'gamma', mean(N(:)));
if dt <= 0 || Lx <= 0 || Ly <= 0
    error('dt, Lx and Ly must be strictly positive.');
end
if minCellCount <= 0
    error('massFluxMinCellCount must be strictly positive.');
end

[Nx, Ny] = size(Ux);
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

[D, A, solver, A0] = cached_mass_flux_operators(Nx, Ny, dx, dy, regularization);

N = double(N);
Ux = double(Ux);
Uy = double(Uy);
Mx = N .* Ux;
My = N .* Uy;
m = [Mx(:); My(:)];

divBeforeVec = D * m;
lowKCorrectionOnly = false;
switch modeLower
    case {"conservative", "divrho", "divrho_zero", "mass_conservative"}
        targetRequested = zeros(Nc, 1);
        target = targetRequested;
        canonicalMode = 'conservative';
    case {"relax_to_uniform", "relax", "density_relax", "rho_relax"}
        targetRaw = (beta / dt) * (N(:) - gamma);
        targetRaw = targetRaw - mean(targetRaw); % total closed/periodic mass flux must sum to zero
        targetRequested = apply_mass_flux_target_filter(targetRaw, Nx, Ny, targetFilter, lowKMaxIndex);
        targetRequested = targetRequested - mean(targetRequested);

        % The collocated D operator used in this prototype can have extra
        % high-frequency null modes on even grids. An arbitrary density
        % relaxation target is therefore not always exactly in range(D).
        % Project the requested target onto range(D) before enforcing it;
        % otherwise the Poisson solve is asked to produce an unattainable
        % checkerboard component and the residual test fails for the wrong
        % reason. The removed component is reported in
        % targetProjectionResidual below.
        target = project_onto_divergence_range(targetRequested, A0, A, solver);
        canonicalMode = 'relax_to_uniform';
    case {"relax_to_uniform_lowk", "relax_lowk", "density_relax_lowk", "rho_relax_lowk"}
        targetRaw = (beta / dt) * (N(:) - gamma);
        targetRaw = targetRaw - mean(targetRaw);
        targetFilter = 'lowpass_fft';
        targetRequested = apply_mass_flux_target_filter(targetRaw, Nx, Ny, targetFilter, lowKMaxIndex);
        targetRequested = targetRequested - mean(targetRequested);
        target = project_onto_divergence_range(targetRequested, A0, A, solver);
        canonicalMode = 'relax_to_uniform_lowk';
        lowKCorrectionOnly = true;
    otherwise
        error('Unknown massFluxProjectionMode: %s', mode);
end

rhsFull = divBeforeVec - target;
if lowKCorrectionOnly
    % Q9 must correct only the low-frequency divergence mismatch.
    % The previous prototype used rhsFull directly, which forced the
    % mass-flux projection to cancel every high-k fluctuation in div(Nu)
    % while asking only for a low-k density target. That over-constrained
    % the velocities and destroyed the local population field.
    rhsRequested = apply_mass_flux_target_filter(rhsFull, Nx, Ny, 'lowpass_fft', lowKMaxIndex);
    rhsRequested = rhsRequested - mean(rhsRequested);
    rhs = project_onto_divergence_range(rhsRequested, A0, A, solver);
else
    rhsRequested = rhsFull;
    rhs = rhsFull;
end
if isempty(solver)
    lambda = A \ rhs;
else
    lambda = solver \ rhs;
end

correction = -D' * lambda;
mProj = m + correction;
divAfterVecFull = D * mProj;
residualVecFull = divAfterVecFull - target;
if lowKCorrectionOnly
    divBeforeDiagVec = apply_mass_flux_target_filter(divBeforeVec, Nx, Ny, 'lowpass_fft', lowKMaxIndex);
    divBeforeDiagVec = divBeforeDiagVec - mean(divBeforeDiagVec);
    divAfterDiagVec = apply_mass_flux_target_filter(divAfterVecFull, Nx, Ny, 'lowpass_fft', lowKMaxIndex);
    divAfterDiagVec = divAfterDiagVec - mean(divAfterDiagVec);
    residualVec = divAfterDiagVec - target;
else
    divBeforeDiagVec = divBeforeVec;
    divAfterDiagVec = divAfterVecFull;
    residualVec = residualVecFull;
end
divAfterVec = divAfterDiagVec;

MxProj = reshape(mProj(1:Nc), [Nx, Ny]);
MyProj = reshape(mProj(Nc+1:end), [Nx, Ny]);
Nsafe = max(N, minCellCount);
UxProj = MxProj ./ Nsafe;
UyProj = MyProj ./ Nsafe;
UxProj(N <= 0) = 0;
UyProj(N <= 0) = 0;

dUx = UxProj - Ux;
dUy = UyProj - Uy;

proj = struct();
proj.mode = canonicalMode;
proj.beta = beta;
proj.N = N;
proj.gamma = gamma;
proj.Mx = Mx;
proj.My = My;
proj.MxProjected = MxProj;
proj.MyProjected = MyProj;
proj.Ux = UxProj;
proj.Uy = UyProj;
proj.dUx = dUx;
proj.dUy = dUy;
proj.lambda = reshape(lambda, [Nx, Ny]);
proj.pi = reshape(lambda / dt, [Nx, Ny]);
proj.targetDivMassRequested = reshape(targetRequested, [Nx, Ny]);
proj.targetDivMass = reshape(target, [Nx, Ny]);
proj.targetProjectionResidual = reshape(target - targetRequested, [Nx, Ny]);
proj.divMassBefore = reshape(divBeforeDiagVec, [Nx, Ny]);
proj.divMassAfter = reshape(divAfterDiagVec, [Nx, Ny]);
proj.divMassResidual = reshape(residualVec, [Nx, Ny]);
proj.divMassBeforeFull = reshape(divBeforeVec, [Nx, Ny]);
proj.divMassAfterFull = reshape(divAfterVecFull, [Nx, Ny]);
proj.divMassResidualFull = reshape(residualVecFull, [Nx, Ny]);
proj.rhsRequested = reshape(rhsRequested, [Nx, Ny]);
proj.rhsUsed = reshape(rhs, [Nx, Ny]);
proj.rmsDivMassBefore = sqrt(mean(divBeforeDiagVec.^2));
proj.rmsDivMassAfter = sqrt(mean(divAfterDiagVec.^2));
proj.rmsTargetDivMassRequested = sqrt(mean(targetRequested.^2));
proj.rmsTargetDivMass = sqrt(mean(target.^2));
proj.rmsTargetProjectionResidual = sqrt(mean((target - targetRequested).^2));
proj.rmsDivMassResidual = sqrt(mean(residualVec.^2));
proj.maxAbsDivMassBefore = max(abs(divBeforeDiagVec));
proj.maxAbsDivMassAfter = max(abs(divAfterDiagVec));
proj.maxAbsTargetDivMassRequested = max(abs(targetRequested));
proj.maxAbsTargetDivMass = max(abs(target));
proj.maxAbsTargetProjectionResidual = max(abs(target - targetRequested));
proj.maxAbsDivMassResidual = max(abs(residualVec));
proj.rmsDivMassBeforeFull = sqrt(mean(divBeforeVec.^2));
proj.rmsDivMassAfterFull = sqrt(mean(divAfterVecFull.^2));
proj.rmsDivMassResidualFull = sqrt(mean(residualVecFull.^2));
proj.maxAbsDivMassBeforeFull = max(abs(divBeforeVec));
proj.maxAbsDivMassAfterFull = max(abs(divAfterVecFull));
proj.maxAbsDivMassResidualFull = max(abs(residualVecFull));
proj.lowKCorrectionOnly = lowKCorrectionOnly;
proj.divMassReduction = proj.rmsDivMassResidual / max(proj.rmsDivMassBefore, eps);
proj.massDeltaBefore = sum(divBeforeVec);
proj.massDeltaAfter = sum(divAfterVecFull);
proj.massDeltaTargetRequested = sum(targetRequested);
proj.massDeltaTarget = sum(target);
proj.dt = dt;
proj.Lx = Lx;
proj.Ly = Ly;
proj.Nx = Nx;
proj.Ny = Ny;
proj.dx = dx;
proj.dy = dy;
proj.regularization = regularization;
proj.minCellCount = minCellCount;
proj.targetFilter = targetFilter;
proj.lowKMaxIndex = lowKMaxIndex;
proj.method = 'mass_flux_DDt_periodic_x_bounded_y';
end


function target = apply_mass_flux_target_filter(targetRaw, Nx, Ny, targetFilter, lowKMaxIndex)
%APPLY_MASS_FLUX_TARGET_FILTER Filter density-relaxation target before projection.
% The low-k mode is deliberately based on fft2, matching the diagnostic
% low-k density energy used elsewhere in this prototype. This is a spectral
% prototype filter; future bounded-y ports can replace it with a cosine or
% channel-aware basis without changing the high-level Q9 algorithm.
filterLower = lower(strrep(char(string(targetFilter)), '-', '_'));
targetRaw = double(targetRaw(:));
switch filterLower
    case {"none", "off", "identity", "raw"}
        target = targetRaw;
    case {"lowpass_fft", "lowk", "fft_lowk"}
        A = reshape(targetRaw, [Nx, Ny]);
        H = fft2(A);
        mask = make_low_k_mask(Nx, Ny, lowKMaxIndex);
        mask(1, 1) = false;
        H(~mask) = 0;
        target = real(ifft2(H));
        target = target(:);
    otherwise
        error('Unknown massFluxTargetFilter: %s', targetFilter);
end
target = target - mean(target);
end

function mask = make_low_k_mask(Nx, Ny, kmax)
ix = [0:floor(Nx/2), -ceil(Nx/2)+1:-1];
iy = [0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix) > Nx
    ix = ix(1:Nx);
end
if numel(iy) > Ny
    iy = iy(1:Ny);
end
[KX, KY] = ndgrid(ix, iy);
mask = sqrt(double(KX).^2 + double(KY).^2) <= kmax;
end

function [D, A, solver, A0] = cached_mass_flux_operators(Nx, Ny, dx, dy, regularization)
persistent cache
key = sprintf('Nx%d_Ny%d_dx%.17g_dy%.17g_reg%.17g', Nx, Ny, dx, dy, regularization);
if ~isempty(cache) && isfield(cache, 'key') && strcmp(cache.key, key)
    D = cache.D;
    A = cache.A;
    solver = cache.solver;
    A0 = cache.A0;
    return;
end

D = build_divergence_operator(Nx, Ny, dx, dy);
A0 = D * D';
A = A0;
if regularization > 0
    A = A + regularization * speye(size(A));
end
solver = [];
try
    if exist('decomposition', 'file') == 2 || exist('decomposition', 'builtin') == 5
        solver = decomposition(A, 'chol');
    end
catch
    solver = [];
end
cache = struct('key', key, 'D', D, 'A', A, 'A0', A0, 'solver', solver);
end

function target = project_onto_divergence_range(targetRequested, A0, A, solver)
%PROJECT_ONTO_DIVERGENCE_RANGE Return the attainable component in range(D).
% For A0 = D*D'', range(A0)=range(D). The small regularized solve avoids
% explicit rank decisions while removing incompatible nullspace modes.
if isempty(solver)
    mu = A \ targetRequested;
else
    mu = solver \ targetRequested;
end
target = A0 * mu;
target = target - mean(target);
end

function D = build_divergence_operator(Nx, Ny, dx, dy)
Nc = Nx * Ny;
maxEntries = Nc * 6;
rows = zeros(maxEntries, 1);
cols = zeros(maxEntries, 1);
vals = zeros(maxEntries, 1);
nzv = 0;

for iy = 1:Ny
    for ix = 1:Nx
        c = grid_index(ix, iy, Nx);

        ixp = ix + 1; if ixp > Nx, ixp = 1; end
        ixm = ix - 1; if ixm < 1, ixm = Nx; end
        cp = grid_index(ixp, iy, Nx);
        cm = grid_index(ixm, iy, Nx);
        [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, cp,  1/(2*dx));
        [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, cm, -1/(2*dx));

        if Ny == 1
            % No y derivative.
        elseif iy == 1
            c0 = Nc + grid_index(ix, iy, Nx);
            c1 = Nc + grid_index(ix, iy+1, Nx);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c1,  1/dy);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c0, -1/dy);
        elseif iy == Ny
            c0 = Nc + grid_index(ix, iy, Nx);
            c1 = Nc + grid_index(ix, iy-1, Nx);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c0,  1/dy);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c1, -1/dy);
        else
            cpY = Nc + grid_index(ix, iy+1, Nx);
            cmY = Nc + grid_index(ix, iy-1, Nx);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, cpY,  1/(2*dy));
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, cmY, -1/(2*dy));
        end
    end
end
D = sparse(rows(1:nzv), cols(1:nzv), vals(1:nzv), Nx*Ny, 2*Nx*Ny);
end

function id = grid_index(ix, iy, Nx)
id = ix + Nx * (iy - 1);
end

function [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, row, col, val)
nzv = nzv + 1;
rows(nzv) = row;
cols(nzv) = col;
vals(nzv) = val;
end

function val = get_param_default(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultValue;
end
end
