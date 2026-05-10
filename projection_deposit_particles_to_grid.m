function G = projection_deposit_particles_to_grid(x, v, params, varargin)
%PROJECTION_DEPOSIT_PARTICLES_TO_GRID Deposit particle velocity on cell centers.
%
%   G = projection_deposit_particles_to_grid(x, v, params)
%
% Minimal particle-to-grid deposit for the pressure-projection prototype.
% The current version uses nearest-cell assignment, consistent with the MPCD
% cell representation. A future version can add bilinear/FLIP-style weights.
%
% Required params fields:
%   Lx, Ly, Nx, Ny
%
% Optional name-value arguments:
%   'periodicX'   default true
%   'periodicY'   default false
%   'minCount'    default 1
%
% Output fields:
%   N, Ux, Uy, Px, Py, rho, valid, dx, dy

if nargin < 3
    error('Usage: G = projection_deposit_particles_to_grid(x, v, params, ...)');
end
if size(x, 2) < 2 || size(v, 2) < 2 || size(x, 1) ~= size(v, 1)
    error('x and v must be N-by-2 arrays with matching particle count.');
end

periodicX = true;
periodicY = false;
minCount = 1;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "mincount"
            minCount = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;

xp = x(:, 1);
yp = x(:, 2);

if periodicX
    xp = mod(xp, Lx);
else
    xp = min(max(xp, 0), Lx - eps(Lx));
end

if periodicY
    yp = mod(yp, Ly);
else
    yp = min(max(yp, 0), Ly - eps(Ly));
end

ix = floor(xp / dx) + 1;
iy = floor(yp / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);

cellId = iy + Ny * (ix - 1);
Nc = Nx * Ny;

N = accumarray(cellId, 1, [Nc, 1], @sum, 0);
Px = accumarray(cellId, v(:, 1), [Nc, 1], @sum, 0);
Py = accumarray(cellId, v(:, 2), [Nc, 1], @sum, 0);

Ux = zeros(Nc, 1);
Uy = zeros(Nc, 1);
valid = N >= minCount;
Ux(valid) = Px(valid) ./ N(valid);
Uy(valid) = Py(valid) ./ N(valid);

G = struct();
G.N = reshape(N, [Ny, Nx]).';
G.Px = reshape(Px, [Ny, Nx]).';
G.Py = reshape(Py, [Ny, Nx]).';
G.Ux = reshape(Ux, [Ny, Nx]).';
G.Uy = reshape(Uy, [Ny, Nx]).';
G.valid = reshape(valid, [Ny, Nx]).';
G.rho = G.N / (dx * dy);
G.dx = dx;
G.dy = dy;
G.Nx = Nx;
G.Ny = Ny;
G.Lx = Lx;
G.Ly = Ly;
end
