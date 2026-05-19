function proj = projection_project_mass_flux_general_bc(N, Ux, Uy, params, varargin)
%PROJECTION_PROJECT_MASS_FLUX_GENERAL_BC General BC-aware Q9 flux projection.
%
%   proj = projection_project_mass_flux_general_bc(N, Ux, Uy, params)
%
% First validation implementation of a finite-volume elliptic operator for
% the Q9 mass-flux correction.  The correction is solved on cell centers but
% represented on faces, so wall/inlet/outlet/periodic boundary conditions act
% on the normal correction flux itself:
%
%       J_new = J + dJ,       dJ = -alpha grad(phi),
%       div(J_new) ~= target.
%
% Boundary convention for the correction dJ:
%   periodic : opposite faces are connected.
%   wall     : dJ_n = 0.
%   inlet    : dJ_n = 0, i.e. Q9 does not alter imposed inflow.
%   outlet   : phi = 0, i.e. pressure-release outlet that may absorb a
%              residual mass-flux incompatibility.
%
% The MATLAB version uses a sparse operator to validate the discretization.
% The operator is assembled from local face stencils, which is the intended
% OpenMP/GPU implementation path.

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
targetFilter = char(string(get_param_default(params, 'massFluxTargetFilter', 'elliptic_lowpass')));
lowKMaxIndex = get_param_default(params, 'massFluxLowKMaxIndex', get_param_default(params, 'lowKMaxIndex', 2));
alphaMode = char(string(get_param_default(params, 'massFluxEllipticAlphaMode', 'constant')));
solverMode = char(string(get_param_default(params, 'massFluxEllipticSolver', 'backslash')));
maxIter = get_param_default(params, 'massFluxEllipticMaxIter', 200);
tol = get_param_default(params, 'massFluxEllipticTol', 1e-8);
maxVelocityKick = get_param_default(params, 'massFluxMaxVelocityKick', Inf);
useSolveCache = get_param_default(params, 'massFluxEllipticUseSolveCache', true);
factorizationMode = char(string(get_param_default(params, 'massFluxEllipticFactorization', 'auto')));

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

N = double(N);
Ux = double(Ux);
Uy = double(Uy);
Mx = N .* Ux;
My = N .* Uy;

bc = parse_mass_flux_bc(params);
validate_periodic_pairs(bc);

alpha = build_face_alpha(N, gamma, alphaMode, minCellCount, bc);
[A0, A, solverInfo] = build_or_get_operator(Nx, Ny, dx, dy, bc, alpha, regularization, solverMode);
faceBase = build_base_face_fluxes(Mx, My, bc, params);
divBeforeVec = divergence_from_faces(faceBase.Jx, faceBase.Jy, dx, dy);

targetProjectionResidual = zeros(Nc, 1);
lowKCorrectionOnly = false;
switch modeLower
    case {"conservative", "divrho", "divrho_zero", "mass_conservative"}
        targetRaw = zeros(Nc, 1);
        targetRequested = targetRaw;
        target = targetRequested;
        canonicalMode = 'conservative';
    case {"relax_to_uniform", "relax", "density_relax", "rho_relax"}
        targetRaw = (beta / dt) * (N(:) - gamma);
        targetRaw = enforce_rhs_compatibility(targetRaw, bc, false);
        targetRequested = apply_general_target_filter(targetRaw, A0, params, targetFilter, lowKMaxIndex, dx, dy, bc);
        target = enforce_rhs_compatibility(targetRequested, bc, false);
        canonicalMode = 'relax_to_uniform';
    case {"relax_to_uniform_lowk", "relax_lowk", "density_relax_lowk", "rho_relax_lowk"}
        targetRaw = (beta / dt) * (N(:) - gamma);
        targetRaw = enforce_rhs_compatibility(targetRaw, bc, false);
        targetRequested = apply_general_target_filter(targetRaw, A0, params, targetFilter, lowKMaxIndex, dx, dy, bc);
        target = enforce_rhs_compatibility(targetRequested, bc, false);
        canonicalMode = 'relax_to_uniform_lowk';
        lowKCorrectionOnly = true;
    otherwise
        error('Unknown massFluxProjectionMode: %s', mode);
end

% With dJ = -alpha grad(phi), the positive finite-volume operator A0 gives
% div(dJ) = A0*phi.  Therefore A0*phi = target - div(J).
rhsFull = target - divBeforeVec;
if lowKCorrectionOnly
    rhsRequested = apply_general_target_filter(rhsFull, A0, params, targetFilter, lowKMaxIndex, dx, dy, bc);
else
    rhsRequested = rhsFull;
end
rhs = enforce_rhs_compatibility(rhsRequested, bc, has_only_neumann_like_bc(bc));
rhsSolve = apply_mode0_gauge_rhs(rhs, solverInfo);

[phi, solveInfo] = solve_projection_system(A, rhsSolve, solverMode, tol, maxIter, ...
    useSolveCache, factorizationMode, solverInfo.cacheKey);
if isfield(solverInfo, 'mode0GaugeActive') && solverInfo.mode0GaugeActive
    % The gauge fixes the additive constant of phi only.  Subtracting the
    % mean after the solve gives a reproducible zero-mean diagnostic
    % potential while preserving exactly the correction fluxes.
    phi = phi - mean(phi);
    physicalResidual = A0 * phi - rhs;
    solveInfo.relresPhysical = norm(physicalResidual) / max(norm(rhs), eps);
    solveInfo.maxAbsPhysicalResidual = max(abs(physicalResidual));
    solveInfo.mode0GaugeActive = true;
    solveInfo.mode0GaugeMode = solverInfo.mode0GaugeMode;
    solveInfo.mode0GaugeIndex = solverInfo.mode0GaugeIndex;
else
    solveInfo.relresPhysical = solveInfo.relres;
    solveInfo.maxAbsPhysicalResidual = NaN;
    solveInfo.mode0GaugeActive = false;
    solveInfo.mode0GaugeMode = 'none';
    solveInfo.mode0GaugeIndex = NaN;
end
faceCorr = correction_faces_from_phi(phi, Nx, Ny, dx, dy, bc, alpha);
faceProjected = faceBase;
faceProjected.Jx = faceProjected.Jx + faceCorr.Jx;
faceProjected.Jy = faceProjected.Jy + faceCorr.Jy;
divAfterVecFull = divergence_from_faces(faceProjected.Jx, faceProjected.Jy, dx, dy);
residualVecFull = divAfterVecFull - target;

if lowKCorrectionOnly
    divBeforeDiagVec = apply_general_target_filter(divBeforeVec, A0, params, targetFilter, lowKMaxIndex, dx, dy, bc);
    divAfterDiagVec = apply_general_target_filter(divAfterVecFull, A0, params, targetFilter, lowKMaxIndex, dx, dy, bc);
    residualVec = divAfterDiagVec - target;
else
    divBeforeDiagVec = divBeforeVec;
    divAfterDiagVec = divAfterVecFull;
    residualVec = residualVecFull;
end

[dMxCell, dMyCell] = cell_center_correction_from_faces(faceCorr.Jx, faceCorr.Jy, bc);
Nsafe = max(N, minCellCount);
dUx = dMxCell ./ Nsafe;
dUy = dMyCell ./ Nsafe;
dUx(N <= 0) = 0;
dUy(N <= 0) = 0;

limiterInfo = apply_velocity_kick_limiter(dUx, dUy, maxVelocityKick);
dUx = limiterInfo.dUx;
dUy = limiterInfo.dUy;
MxProj = Mx + N .* dUx;
MyProj = My + N .* dUy;
UxProj = Ux + dUx;
UyProj = Uy + dUy;
UxProj(N <= 0) = 0;
UyProj(N <= 0) = 0;

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
proj.lambda = reshape(phi, [Nx, Ny]);
proj.phi = proj.lambda;
proj.pi = reshape(phi / dt, [Nx, Ny]);
proj.targetDivMassRequested = reshape(targetRequested, [Nx, Ny]);
proj.targetDivMass = reshape(target, [Nx, Ny]);
proj.targetProjectionResidual = reshape(targetProjectionResidual, [Nx, Ny]);
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
proj.alphaMode = alphaMode;
proj.bc = bc;
proj.faceBase = faceBase;
proj.faceCorrection = faceCorr;
proj.faceProjected = faceProjected;
proj.solverInfo = solveInfo;
proj.operatorInfo = solverInfo;
proj.limiterInfo = limiterInfo.public;
proj.method = 'mass_flux_finite_volume_general_bc_v1';
end

function bc = parse_mass_flux_bc(params)
% Default is the Poiseuille/channel geometry used by the current runner.
bc = struct();
bc.left = make_bc_side('periodic');
bc.right = make_bc_side('periodic');
bc.bottom = make_bc_side('wall');
bc.top = make_bc_side('wall');

if isfield(params, 'massFluxBC') && ~isempty(params.massFluxBC)
    raw = params.massFluxBC;
elseif isfield(params, 'massFluxBoundary') && ~isempty(params.massFluxBoundary)
    raw = params.massFluxBoundary;
else
    raw = [];
end

if ~isempty(raw)
    sides = {'left','right','bottom','top'};
    for s = 1:numel(sides)
        side = sides{s};
        if isstruct(raw) && isfield(raw, side) && ~isempty(raw.(side))
            bc.(side) = parse_bc_side(raw.(side), bc.(side));
        end
    end
end

% Flat aliases are convenient from scripts.
bc.left = parse_flat_bc_alias(params, {'massFluxBCLeft','bcLeft','boundaryLeft','boundary_left'}, bc.left);
bc.right = parse_flat_bc_alias(params, {'massFluxBCRight','bcRight','boundaryRight','boundary_right'}, bc.right);
bc.bottom = parse_flat_bc_alias(params, {'massFluxBCBottom','bcBottom','boundaryBottom','boundary_bottom'}, bc.bottom);
bc.top = parse_flat_bc_alias(params, {'massFluxBCTop','bcTop','boundaryTop','boundary_top'}, bc.top);
end

function side = parse_flat_bc_alias(params, names, current)
side = current;
for k = 1:numel(names)
    nm = names{k};
    if isfield(params, nm) && ~isempty(params.(nm))
        side = parse_bc_side(params.(nm), side);
        return;
    end
end
end

function side = parse_bc_side(raw, current)
side = current;
if ischar(raw) || isstring(raw)
    side.type = canonical_bc_type(raw);
elseif isstruct(raw)
    if isfield(raw, 'type') && ~isempty(raw.type)
        side.type = canonical_bc_type(raw.type);
    end
    if isfield(raw, 'flux') && ~isempty(raw.flux)
        side.flux = raw.flux;
    elseif isfield(raw, 'normalFlux') && ~isempty(raw.normalFlux)
        side.flux = raw.normalFlux;
    elseif isfield(raw, 'massFlux') && ~isempty(raw.massFlux)
        side.flux = raw.massFlux;
    end
elseif isnumeric(raw)
    side.type = 'inlet';
    side.flux = raw;
else
    error('Unsupported boundary condition descriptor.');
end
end

function side = make_bc_side(type)
side = struct();
side.type = canonical_bc_type(type);
side.flux = [];
end

function t = canonical_bc_type(type)
t = lower(strrep(char(string(type)), '-', '_'));
switch t
    case {'periodic','periodic_bc'}
        t = 'periodic';
    case {'wall','solid','no_flux','noflux','neumann'}
        t = 'wall';
    case {'inlet','inflow','in_flow','imposed_inflow'}
        t = 'inlet';
    case {'outlet','outflow','out_flow','open','pressure_release','dirichlet'}
        t = 'outlet';
    otherwise
        error('Unknown mass-flux boundary type: %s', type);
end
end

function validate_periodic_pairs(bc)
if xor(strcmp(bc.left.type,'periodic'), strcmp(bc.right.type,'periodic'))
    error('x-periodic Q9 BC requires both left and right boundaries to be periodic.');
end
if xor(strcmp(bc.bottom.type,'periodic'), strcmp(bc.top.type,'periodic'))
    error('y-periodic Q9 BC requires both bottom and top boundaries to be periodic.');
end
end

function tf = has_only_neumann_like_bc(bc)
% In this first operator, outlets are Dirichlet releases.  Everything else
% is periodic or zero-correction-flux, hence Neumann-like.
tf = ~(strcmp(bc.left.type,'outlet') || strcmp(bc.right.type,'outlet') || ...
       strcmp(bc.bottom.type,'outlet') || strcmp(bc.top.type,'outlet'));
end

function alpha = build_face_alpha(N, gamma, alphaMode, minCellCount, bc)
[Nx, Ny] = size(N);
mode = lower(strrep(alphaMode, '-', '_'));
alpha = struct();
alpha.Jx = ones(Nx+1, Ny);
alpha.Jy = ones(Nx, Ny+1);
if any(strcmp(mode, {'constant','one','uniform'}))
    return;
end
if any(strcmp(mode, {'count','cell_count','density','occupancy'}))
    Nsafe = max(double(N), minCellCount);
    scale = max(double(gamma), eps);
    for j = 1:Ny
        for f = 2:Nx
            alpha.Jx(f,j) = 0.5*(Nsafe(f-1,j)+Nsafe(f,j))/scale;
        end
        if strcmp(bc.left.type,'periodic')
            a = 0.5*(Nsafe(Nx,j)+Nsafe(1,j))/scale;
            alpha.Jx(1,j) = a;
            alpha.Jx(Nx+1,j) = a;
        else
            alpha.Jx(1,j) = Nsafe(1,j)/scale;
            alpha.Jx(Nx+1,j) = Nsafe(Nx,j)/scale;
        end
    end
    for i = 1:Nx
        for f = 2:Ny
            alpha.Jy(i,f) = 0.5*(Nsafe(i,f-1)+Nsafe(i,f))/scale;
        end
        if strcmp(bc.bottom.type,'periodic')
            a = 0.5*(Nsafe(i,Ny)+Nsafe(i,1))/scale;
            alpha.Jy(i,1) = a;
            alpha.Jy(i,Ny+1) = a;
        else
            alpha.Jy(i,1) = Nsafe(i,1)/scale;
            alpha.Jy(i,Ny+1) = Nsafe(i,Ny)/scale;
        end
    end
else
    error('Unknown massFluxEllipticAlphaMode: %s', alphaMode);
end
end

function [A0, A, info] = build_or_get_operator(Nx, Ny, dx, dy, bc, alpha, regularization, solverMode)
% For alpha=constant, a small persistent cache avoids rebuilding in long runs.
% For variable alpha, rebuild: this is a validation implementation.
useCache = is_uniform_alpha(alpha);
isClosedNeumannLike = has_only_neumann_like_bc(bc);
gaugeMode = 'none';
gaugeIndex = NaN;
effectiveRegularization = regularization;
if isClosedNeumannLike
    % Closed periodic/Neumann-like operators have an exact constant null mode.
    % Do not remove it by adding epsilon*I: that changes the elliptic problem
    % and leaves the matrix numerically ill-conditioned.  Instead, fix the
    % additive constant of the potential by pinning one cell.
    gaugeMode = 'pin_first_cell';
    gaugeIndex = 1;
    effectiveRegularization = 0.0;
end
key = '';
persistent cache
if useCache
    key = sprintf('Nx%d_Ny%d_dx%.17g_dy%.17g_bc%s_regEff%.17g_gauge%s_solver%s', ...
        Nx, Ny, dx, dy, bc_signature(bc), effectiveRegularization, gaugeMode, solverMode);
    if ~isempty(cache) && isfield(cache, 'key') && strcmp(cache.key, key)
        A0 = cache.A0;
        A = cache.A;
        info = cache.info;
        return;
    end
end
A0 = build_fv_operator(Nx, Ny, dx, dy, bc, alpha);
A = A0;
if isClosedNeumannLike
    A = apply_mode0_gauge_operator(A, gaugeIndex);
elseif effectiveRegularization > 0
    A = A + effectiveRegularization * speye(Nx*Ny);
end
info = struct();
info.cached = false;
info.cacheable = useCache;
info.cacheKey = key;
info.solverMode = solverMode;
info.bcSignature = bc_signature(bc);
info.nnz = nnz(A);
info.regularization = regularization;
info.effectiveRegularization = effectiveRegularization;
info.hasOnlyNeumannLikeBC = isClosedNeumannLike;
info.mode0GaugeActive = isClosedNeumannLike;
info.mode0GaugeMode = gaugeMode;
info.mode0GaugeIndex = gaugeIndex;
if useCache
    info.cached = true;
    cache = struct('key', key, 'A0', A0, 'A', A, 'info', info);
end
end

function A = apply_mode0_gauge_operator(A, gaugeIndex)
%APPLY_MODE0_GAUGE_OPERATOR Pin one potential value for closed Neumann-like BC.
% This removes only the additive constant null mode of phi.  The correction
% flux dJ = -alpha grad(phi) is unchanged by the choice of gauge.
n = size(A, 1);
if isempty(gaugeIndex) || ~isfinite(gaugeIndex) || gaugeIndex < 1 || gaugeIndex > n
    error('Invalid mode-0 gauge index.');
end
gaugeIndex = round(gaugeIndex);
A(gaugeIndex, :) = sparse(1, n);
A(:, gaugeIndex) = sparse(n, 1);
A(gaugeIndex, gaugeIndex) = 1.0;
end

function rhsSolve = apply_mode0_gauge_rhs(rhs, solverInfo)
rhsSolve = rhs;
if isfield(solverInfo, 'mode0GaugeActive') && solverInfo.mode0GaugeActive
    idx = solverInfo.mode0GaugeIndex;
    if isempty(idx) || ~isfinite(idx)
        error('Missing mode-0 gauge index.');
    end
    rhsSolve(round(idx)) = 0.0;
end
end

function sig = bc_signature(bc)
sig = sprintf('L%s_R%s_B%s_T%s', bc.left.type, bc.right.type, bc.bottom.type, bc.top.type);
end

function tf = is_uniform_alpha(alpha)
tf = all(abs(alpha.Jx(:) - 1) < 10*eps) && all(abs(alpha.Jy(:) - 1) < 10*eps);
end

function A = build_fv_operator(Nx, Ny, dx, dy, bc, alpha)
Nc = Nx * Ny;
rows = zeros(8*Nc + 4*(Nx+Ny), 1);
cols = rows;
vals = rows;
nzv = 0;

% x faces.
if strcmp(bc.left.type, 'periodic')
    for j = 1:Ny
        cL = grid_index(Nx, j, Nx);
        cR = grid_index(1, j, Nx);
        coef = alpha.Jx(1,j) / (dx*dx);
        [rows, cols, vals, nzv] = add_pair_laplacian(rows, cols, vals, nzv, cL, cR, coef);
    end
else
    if strcmp(bc.left.type, 'outlet')
        for j = 1:Ny
            c = grid_index(1, j, Nx);
            coef = 2 * alpha.Jx(1,j) / (dx*dx);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c, coef);
        end
    end
    if strcmp(bc.right.type, 'outlet')
        for j = 1:Ny
            c = grid_index(Nx, j, Nx);
            coef = 2 * alpha.Jx(Nx+1,j) / (dx*dx);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c, coef);
        end
    end
end
for f = 2:Nx
    for j = 1:Ny
        cL = grid_index(f-1, j, Nx);
        cR = grid_index(f, j, Nx);
        coef = alpha.Jx(f,j) / (dx*dx);
        [rows, cols, vals, nzv] = add_pair_laplacian(rows, cols, vals, nzv, cL, cR, coef);
    end
end

% y faces.
if strcmp(bc.bottom.type, 'periodic')
    for i = 1:Nx
        cB = grid_index(i, Ny, Nx);
        cT = grid_index(i, 1, Nx);
        coef = alpha.Jy(i,1) / (dy*dy);
        [rows, cols, vals, nzv] = add_pair_laplacian(rows, cols, vals, nzv, cB, cT, coef);
    end
else
    if strcmp(bc.bottom.type, 'outlet')
        for i = 1:Nx
            c = grid_index(i, 1, Nx);
            coef = 2 * alpha.Jy(i,1) / (dy*dy);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c, coef);
        end
    end
    if strcmp(bc.top.type, 'outlet')
        for i = 1:Nx
            c = grid_index(i, Ny, Nx);
            coef = 2 * alpha.Jy(i,Ny+1) / (dy*dy);
            [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c, c, coef);
        end
    end
end
for f = 2:Ny
    for i = 1:Nx
        cB = grid_index(i, f-1, Nx);
        cT = grid_index(i, f, Nx);
        coef = alpha.Jy(i,f) / (dy*dy);
        [rows, cols, vals, nzv] = add_pair_laplacian(rows, cols, vals, nzv, cB, cT, coef);
    end
end

A = sparse(rows(1:nzv), cols(1:nzv), vals(1:nzv), Nc, Nc);
end

function face = build_base_face_fluxes(Mx, My, bc, params)
[Nx, Ny] = size(Mx);
face = struct();
face.Jx = zeros(Nx+1, Ny);
face.Jy = zeros(Nx, Ny+1);

% x interior faces.
for f = 2:Nx
    face.Jx(f,:) = 0.5 * (Mx(f-1,:) + Mx(f,:));
end
if strcmp(bc.left.type, 'periodic')
    jp = 0.5 * (Mx(Nx,:) + Mx(1,:));
    face.Jx(1,:) = jp;
    face.Jx(Nx+1,:) = jp;
else
    face.Jx(1,:) = boundary_flux_vector('left', bc.left, Mx(1,:), params, Ny);
    face.Jx(Nx+1,:) = boundary_flux_vector('right', bc.right, Mx(Nx,:), params, Ny);
end

% y interior faces.
for f = 2:Ny
    face.Jy(:,f) = 0.5 * (My(:,f-1) + My(:,f));
end
if strcmp(bc.bottom.type, 'periodic')
    jp = 0.5 * (My(:,Ny) + My(:,1));
    face.Jy(:,1) = jp;
    face.Jy(:,Ny+1) = jp;
else
    face.Jy(:,1) = boundary_flux_vector('bottom', bc.bottom, My(:,1), params, Nx);
    face.Jy(:,Ny+1) = boundary_flux_vector('top', bc.top, My(:,Ny), params, Nx);
end
end

function v = boundary_flux_vector(sideName, side, adjacentFlux, params, n)
if strcmp(side.type, 'wall')
    v = zeros(size(adjacentFlux));
elseif strcmp(side.type, 'inlet')
    if ~isempty(side.flux)
        v = expand_flux_value(side.flux, n, sideName);
    else
        v = adjacentFlux;
    end
elseif strcmp(side.type, 'outlet')
    if ~isempty(side.flux)
        v = expand_flux_value(side.flux, n, sideName);
    else
        v = adjacentFlux;
    end
else
    v = adjacentFlux;
end
v = reshape(v, size(adjacentFlux));
end

function v = expand_flux_value(val, n, sideName)
if isscalar(val)
    v = repmat(double(val), [n, 1]);
else
    v = double(val(:));
    if numel(v) ~= n
        error('Boundary flux for %s must be scalar or length %d.', sideName, n);
    end
end
end

function div = divergence_from_faces(Jx, Jy, dx, dy)
[Nxp1, Ny] = size(Jx);
Nx = Nxp1 - 1;
div = zeros(Nx, Ny);
for j = 1:Ny
    for i = 1:Nx
        div(i,j) = (Jx(i+1,j) - Jx(i,j))/dx + (Jy(i,j+1) - Jy(i,j))/dy;
    end
end
div = div(:);
end

function faceCorr = correction_faces_from_phi(phi, Nx, Ny, dx, dy, bc, alpha)
Phi = reshape(phi, [Nx, Ny]);
faceCorr = struct();
faceCorr.Jx = zeros(Nx+1, Ny);
faceCorr.Jy = zeros(Nx, Ny+1);

for f = 2:Nx
    faceCorr.Jx(f,:) = -alpha.Jx(f,:) .* (Phi(f,:) - Phi(f-1,:)) / dx;
end
if strcmp(bc.left.type, 'periodic')
    c = -alpha.Jx(1,:) .* (Phi(1,:) - Phi(Nx,:)) / dx;
    faceCorr.Jx(1,:) = c;
    faceCorr.Jx(Nx+1,:) = c;
else
    if strcmp(bc.left.type, 'outlet')
        faceCorr.Jx(1,:) = -alpha.Jx(1,:) .* (Phi(1,:) - 0) / (0.5*dx);
    end
    if strcmp(bc.right.type, 'outlet')
        faceCorr.Jx(Nx+1,:) = -alpha.Jx(Nx+1,:) .* (0 - Phi(Nx,:)) / (0.5*dx);
    end
end

for f = 2:Ny
    faceCorr.Jy(:,f) = -alpha.Jy(:,f) .* (Phi(:,f) - Phi(:,f-1)) / dy;
end
if strcmp(bc.bottom.type, 'periodic')
    c = -alpha.Jy(:,1) .* (Phi(:,1) - Phi(:,Ny)) / dy;
    faceCorr.Jy(:,1) = c;
    faceCorr.Jy(:,Ny+1) = c;
else
    if strcmp(bc.bottom.type, 'outlet')
        faceCorr.Jy(:,1) = -alpha.Jy(:,1) .* (Phi(:,1) - 0) / (0.5*dy);
    end
    if strcmp(bc.top.type, 'outlet')
        faceCorr.Jy(:,Ny+1) = -alpha.Jy(:,Ny+1) .* (0 - Phi(:,Ny)) / (0.5*dy);
    end
end
end

function [dMx, dMy] = cell_center_correction_from_faces(dJx, dJy, bc)
[Nxp1, Ny] = size(dJx);
Nx = Nxp1 - 1;
dMx = zeros(Nx, Ny);
dMy = zeros(Nx, Ny);
for j = 1:Ny
    for i = 1:Nx
        dMx(i,j) = 0.5 * (dJx(i,j) + dJx(i+1,j));
        dMy(i,j) = 0.5 * (dJy(i,j) + dJy(i,j+1));
    end
end
end

function target = apply_general_target_filter(targetRaw, A0, params, targetFilter, lowKMaxIndex, dx, dy, bc)
filterLower = lower(strrep(char(string(targetFilter)), '-', '_'));
targetRaw = double(targetRaw(:));
switch filterLower
    case {"none", "off", "identity", "raw"}
        target = targetRaw;
    case {"elliptic_lowpass", "operator_lowpass", "lowpass_operator", ...
          "lowpass_elliptic", "lowk_elliptic"}
        target = elliptic_lowpass_filter(targetRaw, A0, params, lowKMaxIndex, dx, dy, bc);
    case {"periodic_fft_lowk", "fft_lowk", "lowpass_fft", "lowk", ...
          "fft_lowk_exact", "periodic_lowpass_fft"}
        target = periodic_fft_lowk_filter_vec(targetRaw, params, lowKMaxIndex, bc);
    otherwise
        error('Unknown massFluxTargetFilter for general BC operator: %s', targetFilter);
end
target = enforce_rhs_compatibility(target, bc, has_only_neumann_like_bc(bc));
end


function filtered = periodic_fft_lowk_filter_vec(raw, params, lowKMaxIndex, bc)
%PERIODIC_FFT_LOWK_FILTER_VEC Exact idempotent low-k FFT mask for fully periodic checks.
%
% This diagnostic filter intentionally reproduces the low-k mask used by
% projection_project_mass_flux_periodic_fft while keeping the FV/general-BC
% elliptic solve for the correction.  It is useful for separating two issues:
%   (i) target subspace/filter differences,
%   (ii) elliptic inverse/discrete-gradient differences.
% It is only valid when both directions are periodic.

if ~is_fully_periodic_bc(bc)
    error('periodic_fft_lowk target filter requires fully periodic mass-flux BCs.');
end
raw = double(raw(:));
Nx = get_param_default(params, 'Nx', []);
Ny = get_param_default(params, 'Ny', []);
if isempty(Nx) || isempty(Ny) || Nx*Ny ~= numel(raw)
    n = round(sqrt(numel(raw)));
    if n*n ~= numel(raw)
        error('Cannot infer periodic FFT filter grid size. Provide params.Nx and params.Ny.');
    end
    Nx = n;
    Ny = n;
end
R = reshape(raw, [Nx, Ny]);
mask = make_periodic_low_k_mask(Nx, Ny, lowKMaxIndex);
mask(1,1) = false;
H = fft2(R);
H(~mask) = 0;
H(1,1) = 0;
filtered = real(ifft2(H));
filtered = filtered(:);
filtered = filtered - mean(filtered);
end

function tf = is_fully_periodic_bc(bc)
tf = strcmp(bc.left.type,'periodic') && strcmp(bc.right.type,'periodic') && ...
     strcmp(bc.bottom.type,'periodic') && strcmp(bc.top.type,'periodic');
end

function mask = make_periodic_low_k_mask(Nx, Ny, kmax)
ix = [0:floor(Nx/2), -ceil(Nx/2)+1:-1];
iy = [0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix) > Nx, ix = ix(1:Nx); end
if numel(iy) > Ny, iy = iy(1:Ny); end
[KX, KY] = ndgrid(ix, iy);
mask = sqrt(double(KX).^2 + double(KY).^2) <= double(kmax);
end

function filtered = elliptic_lowpass_filter(raw, A0, params, lowKMaxIndex, dx, dy, bc)
Nc = numel(raw);
if isfield(params, 'massFluxEllipticLowPassLengthCells') && ~isempty(params.massFluxEllipticLowPassLengthCells)
    lenCells = params.massFluxEllipticLowPassLengthCells;
else
    % Conservative default: a low-k index of 2 on 64 cells corresponds to
    % wavelengths O(20-30 cells).  The Helmholtz radius below damps high-k
    % content without attempting to reproduce the former fft2 mask exactly.
    nApprox = sqrt(double(Nc));
    lenCells = max(1.0, nApprox / max(2*double(lowKMaxIndex) + 1, 1));
end
passes = get_param_default(params, 'massFluxEllipticLowPassPasses', 1);
reg = get_param_default(params, 'massFluxEllipticLowPassRegularization', 0.0);
useSolveCache = get_param_default(params, 'massFluxEllipticUseSolveCache', true);
factorizationMode = char(string(get_param_default(params, 'massFluxEllipticFactorization', 'auto')));
ell = double(lenCells) * 0.5 * (dx + dy);
B = speye(Nc) + (ell^2) * A0;
if reg > 0
    B = B + reg * speye(Nc);
end
bKey = sprintf('LOWPASS_N%d_nnz%d_ell%.17g_reg%.17g_bc%s', ...
    Nc, nnz(A0), ell, reg, bc_signature(bc));
filtered = raw;
for p = 1:passes
    filtered = solve_cached_direct(B, filtered, bKey, useSolveCache, factorizationMode);
    filtered = enforce_rhs_compatibility(filtered, bc, has_only_neumann_like_bc(bc));
end
end

function rhs = enforce_rhs_compatibility(rhs, bc, forceMeanZero)
rhs = double(rhs(:));
if nargin < 3
    forceMeanZero = has_only_neumann_like_bc(bc);
end
if forceMeanZero || has_only_neumann_like_bc(bc)
    rhs = rhs - mean(rhs);
end
end

function [phi, info] = solve_projection_system(A, rhs, solverMode, tol, maxIter, useSolveCache, factorizationMode, cacheKey)
if nargin < 6 || isempty(useSolveCache)
    useSolveCache = true;
end
if nargin < 7 || isempty(factorizationMode)
    factorizationMode = 'auto';
end
if nargin < 8
    cacheKey = '';
end
mode = lower(strrep(char(string(solverMode)), '-', '_'));
info = struct();
info.mode = mode;
info.tol = tol;
info.maxIter = maxIter;
info.flag = NaN;
info.relres = NaN;
info.iter = NaN;
info.usedSolveCache = false;
info.factorizationMode = factorizationMode;
switch mode
    case {'backslash','direct','mldivide','cached','cached_backslash'}
        [phi, cacheInfo] = solve_cached_direct(A, rhs, cacheKey, useSolveCache, factorizationMode);
        info.flag = 0;
        info.relres = norm(A*phi-rhs)/max(norm(rhs), eps);
        info.iter = 1;
        info.usedSolveCache = cacheInfo.usedCache;
        info.factorizationCacheHit = cacheInfo.cacheHit;
    case {'pcg','cg'}
        [phi, flag, relres, iter] = pcg(A, rhs, tol, maxIter);
        info.flag = flag;
        info.relres = relres;
        info.iter = iter;
    otherwise
        error('Unknown massFluxEllipticSolver: %s', solverMode);
end
end

function [x, info] = solve_cached_direct(A, b, cacheKey, useSolveCache, factorizationMode)
%SOLVE_CACHED_DIRECT Cached sparse factorization for repeated MATLAB solves.
% The projection operator and the Helmholtz low-pass operator are usually
% constant throughout a run.  Reusing decomposition(A,...) avoids paying a
% sparse factorization at every MPCD step.  If decomposition is unavailable
% or fails, this helper falls back to A\b.
info = struct('usedCache', false, 'cacheHit', false, 'fallback', false);
if nargin < 3 || isempty(cacheKey)
    cacheKey = '';
end
if nargin < 4 || isempty(useSolveCache)
    useSolveCache = true;
end
if nargin < 5 || isempty(factorizationMode)
    factorizationMode = 'auto';
end
if ~useSolveCache || isempty(cacheKey) || exist('decomposition', 'file') ~= 2
    x = A \ b;
    return;
end
persistent directCache
if isempty(directCache)
    directCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
key = char(cacheKey);
try
    if isKey(directCache, key)
        D = directCache(key);
        info.cacheHit = true;
    else
        D = make_decomposition(A, factorizationMode);
        directCache(key) = D;
    end
    x = D \ b;
    info.usedCache = true;
catch
    % Robust fallback: correctness before speed for unanticipated BC combos.
    x = A \ b;
    info.fallback = true;
end
end

function D = make_decomposition(A, factorizationMode)
mode = lower(strrep(char(string(factorizationMode)), '-', '_'));
switch mode
    case {'chol','cholesky'}
        D = decomposition(A, 'chol');
    case {'ldl'}
        D = decomposition(A, 'ldl');
    case {'lu'}
        D = decomposition(A, 'lu');
    otherwise
        D = decomposition(A, 'auto');
end
end

function limiter = apply_velocity_kick_limiter(dUx, dUy, maxVelocityKick)
limiter = struct();
limiter.dUx = dUx;
limiter.dUy = dUy;
limiter.public = struct('enabled', false, 'maxVelocityKick', maxVelocityKick, ...
    'nLimitedCells', 0, 'maxRawKick', max(sqrt(dUx(:).^2 + dUy(:).^2)));
if isfinite(maxVelocityKick) && maxVelocityKick > 0
    mag = sqrt(dUx.^2 + dUy.^2);
    mask = mag > maxVelocityKick;
    scale = ones(size(mag));
    scale(mask) = maxVelocityKick ./ max(mag(mask), eps);
    limiter.dUx = dUx .* scale;
    limiter.dUy = dUy .* scale;
    limiter.public.enabled = true;
    limiter.public.nLimitedCells = nnz(mask);
end
end

function [rows, cols, vals, nzv] = add_pair_laplacian(rows, cols, vals, nzv, c1, c2, coef)
[rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c1, c1,  coef);
[rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c2, c2,  coef);
[rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c1, c2, -coef);
[rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, c2, c1, -coef);
end

function [rows, cols, vals, nzv] = add_entry(rows, cols, vals, nzv, row, col, val)
nzv = nzv + 1;
if nzv > numel(rows)
    rows = [rows; zeros(numel(rows),1)]; %#ok<AGROW>
    cols = [cols; zeros(numel(cols),1)]; %#ok<AGROW>
    vals = [vals; zeros(numel(vals),1)]; %#ok<AGROW>
end
rows(nzv) = row;
cols(nzv) = col;
vals(nzv) = val;
end

function id = grid_index(ix, iy, Nx)
id = ix + Nx * (iy - 1);
end

function val = get_param_default(params, name, defaultValue)
if nargin < 1 || isempty(params) || ~isstruct(params)
    val = defaultValue;
elseif isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultValue;
end
end
