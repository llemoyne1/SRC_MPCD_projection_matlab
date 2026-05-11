function [stateOut, info] = projection_apply_virial_density_repair(stateIn, params, varargin)
%PROJECTION_APPLY_VIRIAL_DENSITY_REPAIR Local position repair plus velocity restoration.
%
%   [stateOut, info] = projection_apply_virial_density_repair(stateIn, params)
%
% Applies a liquid-like density repair after pressure projection by separating
% the two roles that were mixed in the velocity-only virial kick:
%
%   1. move particles slightly along -grad(Pvir) to repair the configuration,
%      where Pvir = Kvir * (N - gamma);
%   2. restore the cell-mean velocity field to the pressure-projected target
%      field, so the hydrodynamic flow is minimally perturbed by the density
%      repair.
%
% The operation is local on the grid and uses only nearest-cell quantities by
% default. This keeps the method compatible with future OpenMP/MPI/GPU ports.
%
% Main parameters:
%   useVirialDensityRepair              true/false gate
%   virialDensityRepairStrength         dimensionless displacement strength
%   virialDensityRepairK                virial pressure coefficient
%   virialDensityRepairSmoothPasses     smoothing passes applied to Pvir
%   virialDensityRepairMaxDisplacementFraction  limiter in cell-size units
%   virialDensityRepairMinCellCount     denominator floor for N
%   virialDensityRepairRestoreVelocity  restore target U field after repair
%   virialDensityRepairInterpolationMethod       usually 'nearest'
%
% Optional name-value arguments:
%   'periodicX'   default true
%   'periodicY'   default false
%   'targetGrid'  grid deposited from the pressure-projected state. If omitted,
%                 it is computed from stateIn.x/stateIn.v.

if nargin < 2
    error('Usage: [stateOut, info] = projection_apply_virial_density_repair(stateIn, params)');
end
if ~isstruct(stateIn) || ~isfield(stateIn, 'x') || ~isfield(stateIn, 'v')
    error('stateIn must be a struct with fields x and v.');
end
if size(stateIn.x, 1) ~= size(stateIn.v, 1) || size(stateIn.x, 2) < 2 || size(stateIn.v, 2) < 2
    error('stateIn.x and stateIn.v must be Np-by-2 arrays with the same particle count.');
end

periodicX = true;
periodicY = false;
targetGrid = [];
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "targetgrid"
            targetGrid = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

useRepair = logical(get_param(params, 'useVirialDensityRepair', false));
strength = get_param(params, 'virialDensityRepairStrength', 0.0);
Kvir = get_param(params, 'virialDensityRepairK', get_param(params, 'virialK', 1.0));
smoothPasses = round(get_param(params, 'virialDensityRepairSmoothPasses', get_param(params, 'virialSmoothPasses', 0)));
minCellCount = get_param(params, 'virialDensityRepairMinCellCount', get_param(params, 'virialMinCellCount', 1.0));
maxDispFrac = get_param(params, 'virialDensityRepairMaxDisplacementFraction', 0.05);
restoreVelocity = logical(get_param(params, 'virialDensityRepairRestoreVelocity', true));
method = char(string(get_param(params, 'virialDensityRepairInterpolationMethod', 'nearest')));

stateOut = stateIn;
info = empty_info(useRepair, strength, Kvir, smoothPasses, minCellCount, maxDispFrac, restoreVelocity, method);
if ~useRepair || strength == 0 || Kvir == 0
    return;
end

if isempty(targetGrid)
    targetGrid = projection_deposit_particles_to_grid(stateIn.x, stateIn.v, params, ...
        'periodicX', periodicX, 'periodicY', periodicY, 'minCount', 1);
end

Nbefore = double(targetGrid.N);
gamma = get_param(params, 'gamma', mean(Nbefore(:)));
if gamma <= 0
    error('params.gamma must be positive.');
end

Pvir = Kvir * (Nbefore - gamma);
for ipass = 1:max(0, smoothPasses)
    Pvir = smooth_cell_field(Pvir, periodicX, periodicY);
end

[gradPx, gradPy] = gradient_cell_centered(Pvir, params, periodicX, periodicY);
Nden = max(Nbefore, minCellCount);

% Convert the pressure gradient into a bounded configurational displacement.
% Multiplication by dx^2/dy^2 gives a length scale: grad(P) ~ count/length,
% so dx^2*grad(P)/N is O(dx * relative density contrast).
dxCell = params.Lx / params.Nx;
dyCell = params.Ly / params.Ny;
dXGrid = -strength * dxCell^2 * gradPx ./ Nden;
dYGrid = -strength * dyCell^2 * gradPy ./ Nden;

disp = projection_interpolate_grid_delta_to_particles(stateIn.x, dXGrid, dYGrid, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'method', method);

maxDisp = maxDispFrac * min(dxCell, dyCell);
if isfinite(maxDisp) && maxDisp > 0
    mag = sqrt(sum(disp.^2, 2));
    tooLarge = mag > maxDisp;
    if any(tooLarge)
        scale = maxDisp ./ max(mag(tooLarge), eps);
        disp(tooLarge, :) = disp(tooLarge, :) .* scale;
    end
elseif maxDisp <= 0
    disp(:) = 0;
end

xRepaired = stateIn.x + disp;
if periodicX
    xRepaired(:, 1) = mod(xRepaired(:, 1), params.Lx);
else
    xRepaired(:, 1) = min(max(xRepaired(:, 1), 0), params.Lx - eps(params.Lx));
end
if periodicY
    xRepaired(:, 2) = mod(xRepaired(:, 2), params.Ly);
    nWallClamped = 0;
else
    yRaw = xRepaired(:, 2);
    nWallClamped = nnz(yRaw < 0 | yRaw >= params.Ly);
    xRepaired(:, 2) = min(max(yRaw, 0), params.Ly - eps(params.Ly));
end

vRestored = stateIn.v;
GafterMoveRaw = projection_deposit_particles_to_grid(xRepaired, vRestored, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minCount', 1);

UxTarget = targetGrid.Ux;
UyTarget = targetGrid.Uy;
velocityRestoreDeltaRms = 0.0;
velocityRestoreDeltaMaxAbs = 0.0;
velocityRestoreResidualBeforeRms = velocity_field_difference_rms(GafterMoveRaw, targetGrid);
velocityRestoreResidualAfterRms = velocityRestoreResidualBeforeRms;

if restoreVelocity
    dUxRestore = zeros(params.Nx, params.Ny);
    dUyRestore = zeros(params.Nx, params.Ny);
    valid = GafterMoveRaw.N > 0;
    dUxRestore(valid) = UxTarget(valid) - GafterMoveRaw.Ux(valid);
    dUyRestore(valid) = UyTarget(valid) - GafterMoveRaw.Uy(valid);

    dvRestore = projection_interpolate_grid_delta_to_particles(xRepaired, dUxRestore, dUyRestore, params, ...
        'periodicX', periodicX, 'periodicY', periodicY, 'method', 'nearest');
    vRestored = vRestored + dvRestore;
    velocityRestoreDeltaRms = sqrt(mean(sum(dvRestore.^2, 2)));
    velocityRestoreDeltaMaxAbs = max(sqrt(sum(dvRestore.^2, 2)));
else
    dvRestore = zeros(size(vRestored));
end

GafterRestore = projection_deposit_particles_to_grid(xRepaired, vRestored, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minCount', 1);
if restoreVelocity
    velocityRestoreResidualAfterRms = velocity_field_difference_rms(GafterRestore, targetGrid);
end

stateOut.x = xRepaired;
stateOut.v = vRestored;

Nafter = double(GafterRestore.N);
outBandBefore = mean(abs(Nbefore(:) - gamma) > 0.2*gamma);
outBandAfter = mean(abs(Nafter(:) - gamma) > 0.2*gamma);

info.enabled = true;
info.method = method;
info.meanNBefore = mean(Nbefore(:));
info.stdNBefore = std(Nbefore(:));
info.minNBefore = min(Nbefore(:));
info.maxNBefore = max(Nbefore(:));
info.outBand20Before = outBandBefore;
info.meanNAfter = mean(Nafter(:));
info.stdNAfter = std(Nafter(:));
info.minNAfter = min(Nafter(:));
info.maxNAfter = max(Nafter(:));
info.outBand20After = outBandAfter;
info.deltaStdN = info.stdNAfter - info.stdNBefore;
info.deltaOutBand20 = info.outBand20After - info.outBand20Before;
info.pVirMean = mean(Pvir(:));
info.pVirRms = sqrt(mean(Pvir(:).^2));
info.pVirMaxAbs = max(abs(Pvir(:)));
info.gradPvirRms = sqrt(mean(gradPx(:).^2 + gradPy(:).^2));
info.gradPvirMaxAbs = max(sqrt(gradPx(:).^2 + gradPy(:).^2));
info.gridDisplacementRms = sqrt(mean(dXGrid(:).^2 + dYGrid(:).^2));
info.gridDisplacementMaxAbs = max(sqrt(dXGrid(:).^2 + dYGrid(:).^2));
info.particleDisplacementRms = sqrt(mean(sum(disp.^2, 2)));
info.particleDisplacementMaxAbs = max(sqrt(sum(disp.^2, 2)));
info.particleDisplacementMeanX = mean(disp(:, 1));
info.particleDisplacementMeanY = mean(disp(:, 2));
info.nLimitedParticles = nnz(sqrt(sum(disp.^2, 2)) >= maxDisp & isfinite(maxDisp) & maxDisp > 0);
info.nWallClamped = nWallClamped;
info.velocityRestored = restoreVelocity;
info.velocityRestoreDeltaRms = velocityRestoreDeltaRms;
info.velocityRestoreDeltaMaxAbs = velocityRestoreDeltaMaxAbs;
info.velocityRestoreResidualBeforeRms = velocityRestoreResidualBeforeRms;
info.velocityRestoreResidualAfterRms = velocityRestoreResidualAfterRms;
info.NBefore = Nbefore;
info.NAfter = Nafter;
info.Pvir = Pvir;
info.gradPvirX = gradPx;
info.gradPvirY = gradPy;
info.dXGrid = dXGrid;
info.dYGrid = dYGrid;
info.displacement = disp;
info.dvRestore = dvRestore;
info.Gtarget = targetGrid;
info.GafterMoveRaw = GafterMoveRaw;
info.GafterRestore = GafterRestore;
end

function rmsDiff = velocity_field_difference_rms(G, Gtarget)
valid = G.N > 0;
if ~any(valid(:))
    rmsDiff = NaN;
    return;
end
dUx = G.Ux(valid) - Gtarget.Ux(valid);
dUy = G.Uy(valid) - Gtarget.Uy(valid);
rmsDiff = sqrt(mean(dUx(:).^2 + dUy(:).^2));
end

function info = empty_info(useRepair, strength, Kvir, smoothPasses, minCellCount, maxDispFrac, restoreVelocity, method)
info = struct();
info.enabled = false;
info.requested = logical(useRepair);
info.strength = strength;
info.Kvir = Kvir;
info.smoothPasses = smoothPasses;
info.minCellCount = minCellCount;
info.maxDisplacementFraction = maxDispFrac;
info.restoreVelocity = logical(restoreVelocity);
info.method = method;
info.meanNBefore = NaN;
info.stdNBefore = NaN;
info.minNBefore = NaN;
info.maxNBefore = NaN;
info.outBand20Before = NaN;
info.meanNAfter = NaN;
info.stdNAfter = NaN;
info.minNAfter = NaN;
info.maxNAfter = NaN;
info.outBand20After = NaN;
info.deltaStdN = NaN;
info.deltaOutBand20 = NaN;
info.pVirMean = 0.0;
info.pVirRms = 0.0;
info.pVirMaxAbs = 0.0;
info.gradPvirRms = 0.0;
info.gradPvirMaxAbs = 0.0;
info.gridDisplacementRms = 0.0;
info.gridDisplacementMaxAbs = 0.0;
info.particleDisplacementRms = 0.0;
info.particleDisplacementMaxAbs = 0.0;
info.particleDisplacementMeanX = 0.0;
info.particleDisplacementMeanY = 0.0;
info.nLimitedParticles = 0;
info.nWallClamped = 0;
info.velocityRestored = logical(restoreVelocity);
info.velocityRestoreDeltaRms = 0.0;
info.velocityRestoreDeltaMaxAbs = 0.0;
info.velocityRestoreResidualBeforeRms = NaN;
info.velocityRestoreResidualAfterRms = NaN;
info.NBefore = [];
info.NAfter = [];
info.Pvir = [];
info.gradPvirX = [];
info.gradPvirY = [];
info.dXGrid = [];
info.dYGrid = [];
info.displacement = [];
info.dvRestore = [];
info.Gtarget = [];
info.GafterMoveRaw = [];
info.GafterRestore = [];
end

function F = smooth_cell_field(F, periodicX, periodicY)
[Nx, Ny] = size(F);
ixp = [2:Nx, 1];
ixm = [Nx, 1:Nx-1];
if ~periodicX
    ixp(end) = Nx;
    ixm(1) = 1;
end

Fxp = F(ixp, :);
Fxm = F(ixm, :);

Fyp = F;
Fym = F;
if Ny > 1
    Fyp(:, 1:Ny-1) = F(:, 2:Ny);
    Fym(:, 2:Ny) = F(:, 1:Ny-1);
    if periodicY
        Fyp(:, Ny) = F(:, 1);
        Fym(:, 1) = F(:, Ny);
    end
end

F = 0.5*F + 0.125*(Fxp + Fxm + Fyp + Fym);
end

function [dFx, dFy] = gradient_cell_centered(F, params, periodicX, periodicY)
[Nx, Ny] = size(F);
dx = params.Lx / Nx;
dy = params.Ly / Ny;

ixp = [2:Nx, 1];
ixm = [Nx, 1:Nx-1];
if periodicX
    dFx = (F(ixp, :) - F(ixm, :)) / (2*dx);
else
    dFx = zeros(Nx, Ny);
    if Nx > 1
        dFx(1, :) = (F(2, :) - F(1, :)) / dx;
        dFx(Nx, :) = (F(Nx, :) - F(Nx-1, :)) / dx;
    end
    if Nx > 2
        dFx(2:Nx-1, :) = (F(3:Nx, :) - F(1:Nx-2, :)) / (2*dx);
    end
end

dFy = zeros(Nx, Ny);
if Ny > 1
    if periodicY
        iyp = [2:Ny, 1];
        iym = [Ny, 1:Ny-1];
        dFy = (F(:, iyp) - F(:, iym)) / (2*dy);
    else
        dFy(:, 1) = (F(:, 2) - F(:, 1)) / dy;
        dFy(:, Ny) = (F(:, Ny) - F(:, Ny-1)) / dy;
        if Ny > 2
            dFy(:, 2:Ny-1) = (F(:, 3:Ny) - F(:, 1:Ny-2)) / (2*dy);
        end
    end
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
