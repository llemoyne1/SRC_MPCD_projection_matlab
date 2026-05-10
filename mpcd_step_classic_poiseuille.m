function [stateOut, diag] = mpcd_step_classic_poiseuille(state, params)
%MPCD_STEP_CLASSIC_POISEUILLE Classic MPCD step for x-periodic y-wall channel.
%
%   [stateOut, diag] = mpcd_step_classic_poiseuille(state, params)
%
% Standalone kernel for the pressure-projection Poiseuille prototype:
% body-force kick, x-periodic streaming, y-wall collision, then SRD collision.

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

if Lx <= 0 || Ly <= 0 || Nx <= 0 || Ny <= 0 || dt <= 0
    error('Lx, Ly, Nx, Ny and dt must be strictly positive.');
end

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
[x, v, wallInfo] = mpcd_apply_wall_bc_y(x, v, params);

% Collision-cell assignment. The x shift is periodic. The y shift is off by
% default near walls; it can be enabled for experiments, with clamped cells.
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
diag.wallInfo = wallInfo;
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
