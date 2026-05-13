function stepDiag = projection_step_diagnostics(G, params)
%PROJECTION_STEP_DIAGNOSTICS Vorticity and recirculation diagnostics for a rectangular step.
%
%   stepDiag = projection_step_diagnostics(G, params)
%
% The solid step is assumed attached to the lower wall.  Diagnostics focus on
% the downstream shear layer and near-wall recirculation zone behind the step.

Ux = double(G.Ux);
Uy = double(G.Uy);
[Nx, Ny] = size(Ux);
dx = G.dx;
dy = G.dy;
Lx = get_param(params, 'Lx', G.Lx);
Ly = get_param(params, 'Ly', G.Ly);
solidMask = false(Nx, Ny);
if isfield(params, 'solidMask') && ~isempty(params.solidMask)
    solidMask = logical(params.solidMask);
elseif isfield(params, 'solidMaskForProjection') && ~isempty(params.solidMaskForProjection)
    solidMask = logical(params.solidMaskForProjection);
else
    solidMask = mpcd_step_mask(params);
end
if ~isequal(size(solidMask), size(Ux))
    solidMask = false(Nx, Ny);
end
fluidMask = ~solidMask;
Ux(solidMask) = 0;
Uy(solidMask) = 0;

% Periodic x derivative.
dUyDx = (Uy([2:Nx 1], :) - Uy([Nx 1:Nx-1], :)) / (2 * dx);

% Bounded y derivative.
dUxDy = zeros(Nx, Ny);
if Ny >= 3
    dUxDy(:, 2:Ny-1) = (Ux(:, 3:Ny) - Ux(:, 1:Ny-2)) / (2 * dy);
    dUxDy(:, 1) = (Ux(:, 2) - Ux(:, 1)) / dy;
    dUxDy(:, Ny) = (Ux(:, Ny) - Ux(:, Ny-1)) / dy;
elseif Ny == 2
    dUxDy(:, 1) = (Ux(:, 2) - Ux(:, 1)) / dy;
    dUxDy(:, 2) = dUxDy(:, 1);
end
omega = dUyDx - dUxDy;
wf = omega(fluidMask);

x1 = get_param(params, 'stepX1', 0.35 * Lx);
h = get_param(params, 'stepHeight', 0.25 * Ly);
xCenters = ((1:Nx) - 0.5) * dx;
yCenters = ((1:Ny) - 0.5) * dy;
[Xc, Yc] = ndgrid(xCenters, yCenters);
downDist = mod(Xc - x1, Lx);

recircMaxX = get_param(params, 'stepRecirculationMaxDistance', 0.50 * Lx);
recircYMax = get_param(params, 'stepRecirculationYMax', min(2.0*h, 0.55*Ly));
recircMask = fluidMask & downDist > 0 & downDist <= recircMaxX & Yc <= recircYMax;
recircCells = recircMask & Ux < 0;
if any(recircCells(:))
    recirculationLength = max(downDist(recircCells));
    recirculationArea = nnz(recircCells) * dx * dy;
    recirculationMinUx = min(Ux(recircCells));
else
    recirculationLength = 0;
    recirculationArea = 0;
    recirculationMinUx = NaN;
end

shearXMax = get_param(params, 'stepShearLayerMaxDistance', 0.40 * Lx);
shearBand = get_param(params, 'stepShearLayerHalfHeight', max(2*dy, 0.08 * Ly));
shearMask = fluidMask & downDist > 0 & downDist <= shearXMax & abs(Yc - h) <= shearBand;
omegaShear = omega(shearMask);
UxRecirc = Ux(recircMask);
UyRecirc = Uy(recircMask);

probeDist = get_param(params, 'stepProbeDistanceHeights', 3.0) * h;
probeX = mod(x1 + probeDist, Lx);
probeY = min(max(0.5*h, dy/2), Ly - dy/2);
probe = sample_probe(probeX, probeY, Ux, Uy, omega, solidMask, G);

stepDiag = struct();
stepDiag.omega = omega;
stepDiag.solidMask = solidMask;
stepDiag.fluidMask = fluidMask;
stepDiag.meanOmega = mean(wf, 'omitnan');
stepDiag.rmsOmega = sqrt(mean(wf.^2, 'omitnan'));
stepDiag.meanAbsOmega = mean(abs(wf), 'omitnan');
stepDiag.maxAbsOmega = max(abs(wf));
stepDiag.meanEnstrophy = 0.5 * mean(wf.^2, 'omitnan');
stepDiag.totalEnstrophy = 0.5 * sum(wf.^2, 'omitnan') * dx * dy;
stepDiag.shearOmegaRms = sqrt(mean(omegaShear.^2, 'omitnan'));
stepDiag.shearMeanAbsOmega = mean(abs(omegaShear), 'omitnan');
stepDiag.shearCellFraction = nnz(shearMask) / max(nnz(fluidMask), 1);
stepDiag.recirculationLength = recirculationLength;
stepDiag.recirculationArea = recirculationArea;
stepDiag.recirculationCellFraction = nnz(recircCells) / max(nnz(recircMask), 1);
stepDiag.recirculationMinUx = recirculationMinUx;
stepDiag.recirculationMeanUx = mean(UxRecirc, 'omitnan');
stepDiag.recirculationRmsUy = sqrt(mean(UyRecirc.^2, 'omitnan'));
stepDiag.probe = probe;
stepDiag.probeOmega = probe.omega;
stepDiag.probeUx = probe.Ux;
stepDiag.probeUy = probe.Uy;
end

function probe = sample_probe(xq, yq, Ux, Uy, omega, solidMask, G)
Lx = G.Lx;
Ly = G.Ly;
dx = G.dx;
dy = G.dy;
[Nx, Ny] = size(Ux);
xq = mod(xq, Lx);
yq = min(max(yq, 0), Ly - eps(Ly));
ix = floor(xq / dx) + 1;
iy = floor(yq / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
if solidMask(ix, iy)
    [ix, iy] = nearest_fluid_cell(ix, iy, solidMask);
end
probe = struct();
probe.x = (ix - 0.5) * dx;
probe.y = (iy - 0.5) * dy;
probe.ix = ix;
probe.iy = iy;
probe.omega = omega(ix, iy);
probe.Ux = Ux(ix, iy);
probe.Uy = Uy(ix, iy);
probe.isSolid = solidMask(ix, iy);
end

function [ixBest, iyBest] = nearest_fluid_cell(ix0, iy0, solidMask)
[Nx, Ny] = size(solidMask);
[ixGrid, iyGrid] = ndgrid(1:Nx, 1:Ny);
fluid = ~solidMask;
if ~any(fluid(:))
    ixBest = ix0;
    iyBest = iy0;
    return;
end
dxPeriodic = abs(ixGrid - ix0);
dxPeriodic = min(dxPeriodic, Nx - dxPeriodic);
dy = abs(iyGrid - iy0);
d2 = dxPeriodic.^2 + dy.^2;
d2(~fluid) = Inf;
[~, id] = min(d2(:));
[ixBest, iyBest] = ind2sub(size(solidMask), id);
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
