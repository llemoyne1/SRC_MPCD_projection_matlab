function [stateOut, diag] = mpcd_step_classic_step_channel(state, params)
%MPCD_STEP_CLASSIC_STEP_CHANNEL Classic MPCD step with a bottom rectangular step.
%
% Geometry:
%   - periodic x,
%   - bounded y using mpcd_apply_wall_bc_y,
%   - rectangular bottom-attached solid step, grid-aligned by construction.

validate_state(state);

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dt = params.dt;
alphaDeg = get_param(params, 'alphaDeg', 90);
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);
useRandomGridShiftX = logical(get_param(params, 'useRandomGridShiftX', true));
useRandomGridShiftY = logical(get_param(params, 'useRandomGridShiftY', false));

x = state.x;
v = state.v;
Np = size(x, 1);
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

% Force kick and streaming.
v(:, 1) = v(:, 1) + dt * bodyForceX;
v(:, 2) = v(:, 2) + dt * bodyForceY;
x(:, 1) = mod(x(:, 1) + dt * v(:, 1), Lx);
x(:, 2) = x(:, 2) + dt * v(:, 2);
[x, v, wallInfoY] = mpcd_apply_wall_bc_y(x, v, params);
[x, v, stepInfo] = apply_step_bc(x, v, params);

% Collision-cell assignment.
if useRandomGridShiftX
    shiftX = (rand() - 0.5) * dx;
else
    shiftX = 0.0;
end
if useRandomGridShiftY
    shiftY = (rand() - 0.5) * dy;
else
    shiftY = 0.0;
end

xs = mod(x(:, 1) + shiftX, Lx);
ys = min(max(x(:, 2) + shiftY, 0), Ly - eps(Ly));
ix = floor(xs / dx) + 1;
iy = floor(ys / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
cellId = iy + Ny * (ix - 1);

Ncell = accumarray(cellId, 1, [Nc, 1], @sum, 0);
Px = accumarray(cellId, v(:, 1), [Nc, 1], @sum, 0);
Py = accumarray(cellId, v(:, 2), [Nc, 1], @sum, 0);

Ux = zeros(Nc, 1);
Uy = zeros(Nc, 1);
occ = Ncell > 0;
Ux(occ) = Px(occ) ./ Ncell(occ);
Uy(occ) = Py(occ) ./ Ncell(occ);

alpha = alphaDeg * pi / 180.0;
signs = 2.0 * (rand(Nc, 1) > 0.5) - 1.0;
angles = signs * alpha;
ca = cos(angles(cellId));
sa = sin(angles(cellId));

uxp = Ux(cellId);
uyp = Uy(cellId);
rvx = v(:, 1) - uxp;
rvy = v(:, 2) - uyp;

v(:, 1) = uxp + ca .* rvx - sa .* rvy;
v(:, 2) = uyp + sa .* rvx + ca .* rvy;

stateOut = state;
stateOut.x = x;
stateOut.v = v;

diag = struct();
diag.Np = Np;
diag.Nx = Nx;
diag.Ny = Ny;
diag.dx = dx;
diag.dy = dy;
diag.shiftX = shiftX;
diag.shiftY = shiftY;
diag.alphaDeg = alphaDeg;
diag.bodyForceX = bodyForceX;
diag.bodyForceY = bodyForceY;
diag.wallInfoY = wallInfoY;
diag.stepInfo = stepInfo;
diag.wallInfo = combine_wall_info(wallInfoY, stepInfo);
diag.NMean = mean(Ncell);
diag.NStd = std(double(Ncell));
diag.NMin = min(Ncell);
diag.NMax = max(Ncell);
diag.nEmptyCells = nnz(Ncell == 0);
diag.meanVx = mean(v(:, 1));
diag.meanVy = mean(v(:, 2));
diag.kBT = estimate_kBT(v);
diag.kineticEnergyMean = 0.5 * mean(sum(v.^2, 2));
end

function [x, v, info] = apply_step_bc(x, v, params)
Lx = params.Lx;
Ly = params.Ly;
x0 = get_param(params, 'stepX0', 0.15 * Lx);
x1 = get_param(params, 'stepX1', 0.35 * Lx);
h = get_param(params, 'stepHeight', 0.25 * Ly);
mode = lower(strrep(char(string(get_param(params, 'stepWallMode', 'bounceback'))), '-', '_'));
epsWall = get_param(params, 'stepWallEpsilon', 1e-10 * max(Lx, Ly));

oldV = v;
nHitTotal = 0;
maxPen = 0;
for pass = 1:6
    xp = mod(x(:, 1), Lx);
    inside = xp >= x0 & xp <= x1 & x(:, 2) >= 0 & x(:, 2) <= h;
    ids = find(inside);
    if isempty(ids)
        break;
    end
    nHitTotal = nHitTotal + numel(ids);
    dLeft = xp(ids) - x0;
    dRight = x1 - xp(ids);
    dTop = h - x(ids, 2);
    [dmin, face] = min([dLeft, dRight, dTop], [], 2);
    maxPen = max(maxPen, max(dmin));
    nx = zeros(numel(ids), 1);
    ny = zeros(numel(ids), 1);
    for k = 1:numel(ids)
        switch face(k)
            case 1
                x(ids(k), 1) = x0 - epsWall;
                nx(k) = -1;
            case 2
                x(ids(k), 1) = x1 + epsWall;
                nx(k) = 1;
            otherwise
                x(ids(k), 2) = h + epsWall;
                ny(k) = 1;
        end
    end
    x(ids, 1) = mod(x(ids, 1), Lx);
    x(ids, 2) = min(max(x(ids, 2), 0), Ly - eps(Ly));
    switch mode
        case {'specular','slip'}
            vn = v(ids, 1).*nx + v(ids, 2).*ny;
            v(ids, 1) = v(ids, 1) - 2 * vn .* nx;
            v(ids, 2) = v(ids, 2) - 2 * vn .* ny;
        case {'bounceback','no_slip','noslip'}
            v(ids, :) = -v(ids, :);
        otherwise
            error('Unknown stepWallMode: %s', mode);
    end
end

info = struct();
info.enabled = true;
info.mode = mode;
info.nStepHits = nHitTotal;
info.dPxStep = sum(v(:,1) - oldV(:,1));
info.dPyStep = sum(v(:,2) - oldV(:,2));
info.maxPenetration = maxPen;
end

function wallInfo = combine_wall_info(wallInfoY, stepInfo)
wallInfo = wallInfoY;
wallInfo.nStep = stepInfo.nStepHits;
wallInfo.dPxStep = stepInfo.dPxStep;
wallInfo.dPyStep = stepInfo.dPyStep;
wallInfo.stepMode = stepInfo.mode;
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v')
    error('state must be a struct with fields x and v.');
end
if size(state.x, 2) ~= 2 || size(state.v, 2) ~= 2 || size(state.x, 1) ~= size(state.v, 1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
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
