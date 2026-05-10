function th = projection_thermal_diagnostics(x, v, params, varargin)
%PROJECTION_THERMAL_DIAGNOSTICS Separate hydrodynamic and thermal energies.
%
%   th = projection_thermal_diagnostics(x, v, params)
%
% Deposits particle velocities on the MPCD grid, interpolates the cell mean
% velocity back to particles with nearest-cell assignment, and separates:
%   hydrodynamic energy   1/2 <|U_cell|^2>
%   thermal energy        1/2 <|v_i - U_cell(i)|^2>
%
% The global kBT diagnostic based on v_i - mean(v) is also reported for
% comparison. In 2D, kBT = 1/2 <|c|^2> for unit particle mass.

if nargin < 3
    error('Usage: th = projection_thermal_diagnostics(x, v, params, ...)');
end
if size(x,2) ~= 2 || size(v,2) ~= 2 || size(x,1) ~= size(v,1)
    error('x and v must be Np-by-2 arrays with matching particle count.');
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

G = projection_deposit_particles_to_grid(x, v, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minCount', minCount);
ids = local_cell_ids(x, params, periodicX, periodicY);

Ux = G.Ux(:);
Uy = G.Uy(:);
Ucell = [Ux(ids), Uy(ids)];
vrel = v - Ucell;
umean = mean(v, 1);
cglobal = v - umean;

N = G.N(:);
validCells = N > 0;
localKBT = nan(numel(N), 1);
localHydroKE = nan(numel(N), 1);
for c = find(validCells).'
    pids = ids == c;
    if any(pids)
        cc = v(pids, :) - [Ux(c), Uy(c)];
        localKBT(c) = 0.5 * mean(sum(cc.^2, 2));
        localHydroKE(c) = 0.5 * (Ux(c)^2 + Uy(c)^2);
    end
end

th = struct();
th.Np = size(x,1);
th.G = G;
th.meanVx = umean(1);
th.meanVy = umean(2);
th.kBTGlobal = 0.5 * mean(sum(cglobal.^2, 2));
th.kBTCellRelative = 0.5 * mean(sum(vrel.^2, 2));
th.totalKineticEnergy = 0.5 * mean(sum(v.^2, 2));
th.hydroKineticEnergy = 0.5 * mean(sum(Ucell.^2, 2));
th.thermalKineticEnergy = th.kBTCellRelative;
th.energyClosureError = th.totalKineticEnergy - th.hydroKineticEnergy - th.thermalKineticEnergy;
th.localKBTMean = mean(localKBT(validCells), 'omitnan');
th.localKBTStd = std(localKBT(validCells), 0, 'omitnan');
th.localKBTMin = min(localKBT(validCells), [], 'omitnan');
th.localKBTMax = max(localKBT(validCells), [], 'omitnan');
th.localHydroKEMean = mean(localHydroKE(validCells), 'omitnan');
th.nEmptyCells = nnz(~validCells);
th.nValidCells = nnz(validCells);
th.populationStd = std(double(N));
th.populationMean = mean(double(N));
th.populationCV = th.populationStd / max(th.populationMean, eps);
end

function ids = local_cell_ids(x, params, periodicX, periodicY)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;

xp = x(:,1);
yp = x(:,2);
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
ids = iy + Ny * (ix - 1);
end
