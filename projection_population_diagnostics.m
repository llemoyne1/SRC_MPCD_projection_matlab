function pop = projection_population_diagnostics(x, params, varargin)
%PROJECTION_POPULATION_DIAGNOSTICS Cell-population diagnostics for projection runs.
%
%   pop = projection_population_diagnostics(x, params)
%
% Counts particles on the same cell grid used by the pressure-projection
% prototype. The function is deliberately independent from velocity deposit so
% it can be used before and after a projection step. Since pressure projection
% only changes velocities, a correct implementation must leave these diagnostics
% unchanged across the projection substep.
%
% Required params fields:
%   Lx, Ly, Nx, Ny
%
% Optional params fields:
%   gamma       target mean occupancy used for band diagnostics
%
% Optional name-value arguments:
%   'periodicX'      default true
%   'periodicY'      default false
%   'targetGamma'    default params.gamma if present, otherwise mean(N)
%   'bandFraction'   default 0.20
%
% Output fields include N, meanN, stdN, minN, maxN, nEmptyCells,
% fracEmptyCells, coefficientOfVariation, outBandFraction, xProfile,
% yProfile and entropyLike. N is an Nx-by-Ny array.

if nargin < 2
    error('Usage: pop = projection_population_diagnostics(x, params, ...)');
end
if size(x, 2) < 2
    error('x must be an Np-by-2 array.');
end

periodicX = true;
periodicY = false;
targetGamma = [];
bandFraction = 0.20;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "targetgamma"
            targetGamma = val;
        case "bandfraction"
            bandFraction = val;
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
Nvec = accumarray(cellId, 1, [Nc, 1], @sum, 0);
N = reshape(Nvec, [Ny, Nx]).';

if isempty(targetGamma)
    if isfield(params, 'gamma') && ~isempty(params.gamma)
        targetGamma = params.gamma;
    else
        targetGamma = mean(Nvec);
    end
end

lo = targetGamma * (1 - bandFraction);
hi = targetGamma * (1 + bandFraction);
meanN = mean(Nvec);
stdN = std(double(Nvec));

p = double(Nvec) / max(sum(Nvec), eps);
p = p(p > 0);
entropyLike = -sum(p .* log(p));

pop = struct();
pop.N = N;
pop.Nx = Nx;
pop.Ny = Ny;
pop.Np = size(x, 1);
pop.meanN = meanN;
pop.stdN = stdN;
pop.minN = min(Nvec);
pop.maxN = max(Nvec);
pop.nEmptyCells = nnz(Nvec == 0);
pop.fracEmptyCells = nnz(Nvec == 0) / max(numel(Nvec), 1);
pop.coefficientOfVariation = stdN / max(meanN, eps);
pop.targetGamma = targetGamma;
pop.bandFraction = bandFraction;
pop.outBandFraction = nnz(Nvec < lo | Nvec > hi) / max(numel(Nvec), 1);
pop.xProfile = mean(N, 2);
pop.yProfile = mean(N, 1).';
pop.entropyLike = entropyLike;
pop.periodicX = periodicX;
pop.periodicY = periodicY;
end
