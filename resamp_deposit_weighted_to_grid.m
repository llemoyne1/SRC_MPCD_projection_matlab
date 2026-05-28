function G = resamp_deposit_weighted_to_grid(x, v, m, params, varargin)
%RESAMP_DEPOSIT_WEIGHTED_TO_GRID Deposit weighted particles on cell centers.
%
%   G = resamp_deposit_weighted_to_grid(x, v, m, params)
%
% The historical deposit uses unit particle mass.  This weighted deposit
% separates population N from cell mass M and computes the hydrodynamic
% velocity as U = P/M, with P = sum_p m_p v_p.
%
% Optional name-value arguments:
%   'periodicX'   default true
%   'periodicY'   default false
%   'minMass'     default eps
%
% Output fields include:
%   N, M, Px, Py, Ux, Uy, rho, valid, cellId, dx, dy

if nargin < 4
    error('Usage: G = resamp_deposit_weighted_to_grid(x, v, m, params, ...)');
end
if size(x, 2) < 2 || size(v, 2) < 2 || size(x, 1) ~= size(v, 1)
    error('x and v must be Np-by-2 arrays with matching particle count.');
end
m = m(:);
if numel(m) ~= size(x, 1)
    error('m must have one entry per particle.');
end
if any(~isfinite(m)) || any(m < 0)
    error('Particle masses must be finite and non-negative.');
end

periodicX = true;
periodicY = false;
minMass = eps;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case {"minmass", "minm"}
            minMass = val;
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

cellId = resamp_cell_ids_periodic(x, params, 'periodicX', periodicX, 'periodicY', periodicY);
Nc = Nx * Ny;

N = accumarray(cellId, 1, [Nc, 1], @sum, 0);
M = accumarray(cellId, m, [Nc, 1], @sum, 0);
Px = accumarray(cellId, m .* v(:, 1), [Nc, 1], @sum, 0);
Py = accumarray(cellId, m .* v(:, 2), [Nc, 1], @sum, 0);

Ux = zeros(Nc, 1);
Uy = zeros(Nc, 1);
valid = M > minMass;
Ux(valid) = Px(valid) ./ M(valid);
Uy(valid) = Py(valid) ./ M(valid);

G = struct();
G.N = reshape(N, [Ny, Nx]).';
G.M = reshape(M, [Ny, Nx]).';
G.Px = reshape(Px, [Ny, Nx]).';
G.Py = reshape(Py, [Ny, Nx]).';
G.Ux = reshape(Ux, [Ny, Nx]).';
G.Uy = reshape(Uy, [Ny, Nx]).';
G.valid = reshape(valid, [Ny, Nx]).';
G.rho = G.M / (dx * dy);
G.cellId = cellId;
G.dx = dx;
G.dy = dy;
G.Nx = Nx;
G.Ny = Ny;
G.Lx = Lx;
G.Ly = Ly;
G.totalMass = sum(m, 'omitnan');
G.totalMomentum = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];
end
