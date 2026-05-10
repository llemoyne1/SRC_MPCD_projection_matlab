function [stateOut, diag] = mpcd_step_classic_periodic(state, params)
%MPCD_STEP_CLASSIC_PERIODIC One classic 2D SRD/MPCD step in a periodic box.
%
%   [stateOut, diag] = mpcd_step_classic_periodic(state, params)
%
% This is a compact, standalone kernel for the pressure-projection prototype.
% It intentionally does not call the legacy redistribution/liquid-closure code.
%
% State fields:
%   state.x   Np-by-2 particle positions in [0,Lx) x [0,Ly)
%   state.v   Np-by-2 particle velocities
%
% Required params fields:
%   Lx, Ly, Nx, Ny, dt
%
% Optional params fields:
%   alphaDeg            SRD rotation angle in degrees, default 90
%   bodyForceX          acceleration in x, default 0
%   bodyForceY          acceleration in y, default 0
%   useRandomGridShift  random collision-grid shift, default true
%
% The order is:
%   force kick -> periodic streaming -> SRD collision.
%
% The collision conserves cell momentum exactly and uses one randomly signed
% rotation angle per occupied cell.

validate_state(state);

Lx = get_param(params, 'Lx', []);
Ly = get_param(params, 'Ly', []);
Nx = get_param(params, 'Nx', []);
Ny = get_param(params, 'Ny', []);
dt = get_param(params, 'dt', []);
alphaDeg = get_param(params, 'alphaDeg', 90);
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);
useRandomGridShift = logical(get_param(params, 'useRandomGridShift', true));

if isempty(Lx) || isempty(Ly) || isempty(Nx) || isempty(Ny) || isempty(dt)
    error('params must contain Lx, Ly, Nx, Ny and dt.');
end
if Lx <= 0 || Ly <= 0 || Nx <= 0 || Ny <= 0 || dt <= 0
    error('Lx, Ly, Nx, Ny and dt must be strictly positive.');
end

x = state.x;
v = state.v;
Np = size(x, 1);

dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

% Body force and periodic streaming.
v(:, 1) = v(:, 1) + dt * bodyForceX;
v(:, 2) = v(:, 2) + dt * bodyForceY;
x(:, 1) = mod(x(:, 1) + dt * v(:, 1), Lx);
x(:, 2) = mod(x(:, 2) + dt * v(:, 2), Ly);

% Collision-cell assignment, with optional random grid shift for Galilean
% invariance in the small-mean-free-path regime.
if useRandomGridShift
    shiftX = (rand() - 0.5) * dx;
    shiftY = (rand() - 0.5) * dy;
else
    shiftX = 0.0;
    shiftY = 0.0;
end

xs = mod(x(:, 1) + shiftX, Lx);
ys = mod(x(:, 2) + shiftY, Ly);
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
