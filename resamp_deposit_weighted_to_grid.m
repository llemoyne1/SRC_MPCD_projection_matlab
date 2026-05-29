function G = resamp_deposit_weighted_to_grid(x, v, m, params, varargin)
%RESAMP_DEPOSIT_WEIGHTED_TO_GRID Deposit weighted active particles on cell centers.
%
%   G = resamp_deposit_weighted_to_grid(x, v, m, params)
%
% Optional name-value arguments:
%   'periodicX'   default true
%   'periodicY'   default false
%   'minMass'     default eps
%   'activeMask'  default all true
%   'cellWetMask' default state/params/all-wet; stored in G.cellWetMask
%
% Output fields include N, M, Px, Py, Ux, Uy, rho, valid, cellId.  cellId has
% one entry per storage row; inactive rows are set to zero and do not enter
% accumarray deposits.

if nargin < 4
    error('Usage: G = resamp_deposit_weighted_to_grid(x, v, m, params, ...)');
end
if size(x, 2) < 2 || size(v, 2) < 2 || size(x, 1) ~= size(v, 1)
    error('x and v must be Np-by-2 arrays with matching particle count.');
end
m = m(:);
Np = size(x, 1);
if numel(m) ~= Np
    error('m must have one entry per particle slot.');
end
if any(~isfinite(m)) || any(m < 0)
    error('Particle masses must be finite and non-negative.');
end

periodicX = true;
periodicY = false;
minMass = eps;
activeMask = true(Np, 1);
cellWetMask = [];
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
        case {"activemask", "active"}
            activeMask = logical(val(:));
        case {"cellwetmask", "wetmask", "fluidmask"}
            cellWetMask = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end
if numel(activeMask) ~= Np
    error('activeMask must have one entry per particle slot.');
end
activeMask = activeMask & all(isfinite(x),2) & all(isfinite(v),2) & isfinite(m) & m >= 0;

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

cellId = zeros(Np, 1);
if any(activeMask)
    cellId(activeMask) = resamp_cell_ids_periodic(x(activeMask,:), params, ...
        'periodicX', periodicX, 'periodicY', periodicY);
    idA = cellId(activeMask);
    mA = m(activeMask);
    vA = v(activeMask, :);
    N = accumarray(idA, 1, [Nc, 1], @sum, 0);
    M = accumarray(idA, mA, [Nc, 1], @sum, 0);
    Px = accumarray(idA, mA .* vA(:, 1), [Nc, 1], @sum, 0);
    Py = accumarray(idA, mA .* vA(:, 2), [Nc, 1], @sum, 0);
else
    N = zeros(Nc, 1);
    M = zeros(Nc, 1);
    Px = zeros(Nc, 1);
    Py = zeros(Nc, 1);
end

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
G.activeMask = activeMask;
G.activeIndex = find(activeMask);
G.dx = dx;
G.dy = dy;
[wetMask, wetInfo] = resamp_cell_wet_mask(struct('cellWetMask', cellWetMask), params, ...
    'mode', ternary_empty(cellWetMask, 'auto', 'explicit'), 'cellWetMask', cellWetMask);
G.cellWetMask = wetMask;
G.nWetCells = wetInfo.nWetCells;
G.nDryCells = wetInfo.nDryCells;
G.wetFraction = wetInfo.wetFraction;
G.Nx = Nx;
G.Ny = Ny;
G.Lx = Lx;
G.Ly = Ly;
G.totalMass = sum(m(activeMask), 'omitnan');
if any(activeMask)
    G.totalMomentum = [sum(m(activeMask) .* v(activeMask,1), 'omitnan'), ...
                       sum(m(activeMask) .* v(activeMask,2), 'omitnan')];
else
    G.totalMomentum = [0 0];
end
G.Nactive = nnz(activeMask);
G.Ncapacity = Np;
G.Nfree = Np - nnz(activeMask);
end

function out = ternary_empty(x, a, b)
if isempty(x)
    out = a;
else
    out = b;
end
end
