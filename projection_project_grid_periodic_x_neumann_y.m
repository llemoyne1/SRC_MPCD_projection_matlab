function proj = projection_project_grid_periodic_x_neumann_y(Ux, Uy, params, varargin)
%PROJECTION_PROJECT_GRID_PERIODIC_X_NEUMANN_Y Discrete pressure projection with walls.
%
%   proj = projection_project_grid_periodic_x_neumann_y(Ux, Uy, params)
%
% Projects a cell-centered velocity field on a domain periodic in x and
% bounded in y. This first wall-compatible prototype uses an algebraic
% Helmholtz-Hodge projection:
%
%       u_projected = u - D' lambda,
%       (D D') lambda = D u,
%
% where D is a sparse discrete divergence operator. This guarantees that the
% post-projection discrete divergence D*u is reduced to solver precision. The
% pressure-like variable is p = rho/dt * lambda.
%
% The y treatment is a simple bounded-domain finite-difference divergence:
% centered differences in the interior and one-sided differences on the first
% and last rows. It is intended as a robust MATLAB prototype before moving to
% a MAC-grid or GPU implementation.
%
% Required params fields:
%   Lx, Ly, dt
%
% Optional params/arguments:
%   params.rho0 or name-value 'rho'      default 1
%   'regularization'                     default 1e-12
%
% Output fields:
%   Ux, Uy, dUx, dUy, p, divBefore, divAfter, rmsDivBefore, rmsDivAfter.

if nargin < 3 || isempty(params)
    params = struct();
end
if ~isnumeric(Ux) || ~isnumeric(Uy) || ~isequal(size(Ux), size(Uy))
    error('Ux and Uy must be numeric arrays with identical size.');
end

rho = get_param_default(params, 'rho0', 1.0);
regularization = 1e-12;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "rho"
            rho = val;
        case "regularization"
            regularization = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Lx = get_param_default(params, 'Lx', size(Ux, 1));
Ly = get_param_default(params, 'Ly', size(Ux, 2));
dt = get_param_default(params, 'dt', 1.0);
if rho <= 0 || dt <= 0 || Lx <= 0 || Ly <= 0
    error('rho, dt, Lx and Ly must be strictly positive.');
end

[Nx, Ny] = size(Ux);
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

[D, A, solver] = cached_projection_operators(Nx, Ny, dx, dy, regularization);

u = [Ux(:); Uy(:)];
divBeforeVec = D * u;
rhs = divBeforeVec;

% Solve A*lambda = D*u. A tiny Tikhonov regularization, set when the
% cached operator is built, removes the pressure gauge without imposing a
% point constraint. The previous lambda(1)=0 constraint left a residual
% divergence for the bounded-y operator.
if isempty(solver)
    lambda = A \ rhs;
else
    lambda = solver \ rhs;
end
% lambda mean subtraction intentionally disabled: it changes D''*lambda for the bounded-y operator.

correction = -D' * lambda;
uProj = u + correction;
divAfterVec = D * uProj;

UxProj = reshape(uProj(1:Nc), [Nx, Ny]);
UyProj = reshape(uProj(Nc+1:end), [Nx, Ny]);
dUx = UxProj - Ux;
dUy = UyProj - Uy;

proj = struct();
proj.Ux = UxProj;
proj.Uy = UyProj;
proj.dUx = dUx;
proj.dUy = dUy;
proj.p = reshape((rho/dt) * lambda, [Nx, Ny]);
proj.lambda = reshape(lambda, [Nx, Ny]);
proj.divBefore = reshape(divBeforeVec, [Nx, Ny]);
proj.divAfter = reshape(divAfterVec, [Nx, Ny]);
proj.meanDivBefore = mean(divBeforeVec);
proj.rmsDivBefore = sqrt(mean(divBeforeVec.^2));
proj.rmsDivAfter = sqrt(mean(divAfterVec.^2));
proj.maxAbsDivBefore = max(abs(divBeforeVec));
proj.maxAbsDivAfter = max(abs(divAfterVec));
proj.divReduction = proj.rmsDivAfter / max(proj.rmsDivBefore, eps);
proj.rho = rho;
proj.dt = dt;
proj.Lx = Lx;
proj.Ly = Ly;
proj.Nx = Nx;
proj.Ny = Ny;
proj.dx = dx;
proj.dy = dy;
proj.regularization = regularization;
proj.method = 'algebraic_DDt_periodic_x_bounded_y';
end

function [D, A, solver] = cached_projection_operators(Nx, Ny, dx, dy, regularization)
persistent cache
key = sprintf('Nx%d_Ny%d_dx%.17g_dy%.17g_reg%.17g', Nx, Ny, dx, dy, regularization);
if ~isempty(cache) && isfield(cache, 'key') && strcmp(cache.key, key)
    D = cache.D;
    A = cache.A;
    solver = cache.solver;
    return;
end

D = build_divergence_operator(Nx, Ny, dx, dy);
A = D * D';
if regularization > 0
    A = A + regularization * speye(size(A));
end

% A = D*D' has a pressure-gauge null mode. The small regularization
% removes this gauge. Do not impose a point constraint here: replacing one
% row/column is not a true orthogonal projection for this bounded-y
% divergence operator and leaves a measurable residual divergence.

solver = [];
try
    if exist('decomposition', 'file') == 2 || exist('decomposition', 'builtin') == 5
        solver = decomposition(A, 'chol');
    end
catch
    solver = [];
end

cache = struct('key', key, 'D', D, 'A', A, 'solver', solver);
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

        % dUx/dx: centered periodic difference.
        ixp = ix + 1; if ixp > Nx, ixp = 1; end
        ixm = ix - 1; if ixm < 1, ixm = Nx; end
        cp = grid_index(ixp, iy, Nx);
        cm = grid_index(ixm, iy, Nx);
        [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, cp,  1/(2*dx));
        [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, cm, -1/(2*dx));

        % dUy/dy: bounded y. Use one-sided differences at the two walls and
        % centered differences in the interior. Uy unknowns live after all Ux
        % unknowns in the vector [Ux(:); Uy(:)].
        if Ny == 1
            % No y derivative possible.
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

D = sparse(rows(1:nzv), cols(1:nzv), vals(1:nzv), Nc, 2*Nc);
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
