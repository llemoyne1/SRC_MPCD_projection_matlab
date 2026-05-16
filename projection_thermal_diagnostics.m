function th = projection_thermal_diagnostics(x, v, params, varargin)
%PROJECTION_THERMAL_DIAGNOSTICS Separate hydrodynamic and thermal energies.
%
%   th = projection_thermal_diagnostics(x, v, params, ...)
%
% This diagnostic uses the same nearest-cell assignment convention as
% projection_apply_cell_thermostat.  Important: projection_deposit_particles_to_grid
% stores fields as Nx-by-Ny arrays after a transpose for plotting/field usage;
% therefore this function does NOT recover cell values using G.Ux(:) indexed by
% the raw particle cell ids.  It recomputes the cell means in the native 1D cell
% id ordering, so that the temperature diagnosed here is exactly comparable to
% the one imposed by projection_apply_cell_thermostat.
%
% In 2D, the reported cell-relative kBT is
%      kBT = 1/2 < |v_i - U_cell(i)|^2 >
% for unit particle mass.

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

% Keep the usual grid structure for downstream diagnostics/plots.
G = projection_deposit_particles_to_grid(x, v, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minCount', minCount);

ids = local_cell_ids(x, params, periodicX, periodicY);
Nc = params.Nx * params.Ny;
Np = size(x,1);

% Native 1D cell statistics, consistent with the thermostat.
nCell = accumarray(ids, 1, [Nc 1], @sum, 0);
sumVx = accumarray(ids, v(:,1), [Nc 1], @sum, 0);
sumVy = accumarray(ids, v(:,2), [Nc 1], @sum, 0);

populated = nCell > 0;
Ux1 = zeros(Nc,1);
Uy1 = zeros(Nc,1);
Ux1(populated) = sumVx(populated) ./ nCell(populated);
Uy1(populated) = sumVy(populated) ./ nCell(populated);

Ucell = [Ux1(ids), Uy1(ids)];
vrel = v - Ucell;
umean = mean(v, 1, 'omitnan');
cglobal = v - umean;

rel2 = sum(vrel.^2, 2);
sumRel2 = accumarray(ids, rel2, [Nc 1], @sum, 0);
localKBT = nan(Nc,1);
localKBT(populated) = 0.5 * sumRel2(populated) ./ nCell(populated);

localHydroKE = nan(Nc,1);
localHydroKE(populated) = 0.5 * (Ux1(populated).^2 + Uy1(populated).^2);

validCells = nCell >= minCount;

th = struct();
th.Np = Np;
th.G = G;
th.meanVx = umean(1);
th.meanVy = umean(2);
th.kBTGlobal = 0.5 * mean(sum(cglobal.^2, 2), 'omitnan');
th.kBTCellRelative = 0.5 * mean(rel2, 'omitnan');
th.totalKineticEnergy = 0.5 * mean(sum(v.^2, 2), 'omitnan');
th.hydroKineticEnergy = 0.5 * mean(sum(Ucell.^2, 2), 'omitnan');
th.thermalKineticEnergy = th.kBTCellRelative;
th.energyClosureError = th.totalKineticEnergy - th.hydroKineticEnergy - th.thermalKineticEnergy;
th.localKBTMean = mean(localKBT(validCells), 'omitnan');
th.localKBTStd = std(localKBT(validCells), 0, 'omitnan');
th.localKBTMin = min(localKBT(validCells), [], 'omitnan');
th.localKBTMax = max(localKBT(validCells), [], 'omitnan');
th.localHydroKEMean = mean(localHydroKE(validCells), 'omitnan');
th.nEmptyCells = nnz(~populated);
th.nValidCells = nnz(validCells);
th.populationStd = std(double(nCell));
th.populationMean = mean(double(nCell));
th.populationCV = th.populationStd / max(th.populationMean, eps);

% Extra consistency fields useful while validating the thermostat.
if any(validCells)
    th.kBTCellRelativeWeightedValid = sum(nCell(validCells) .* localKBT(validCells)) ./ sum(nCell(validCells));
else
    th.kBTCellRelativeWeightedValid = NaN;
end
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
