function vort = projection_vorticity_diagnostics(G, params)
%PROJECTION_VORTICITY_DIAGNOSTICS Vorticity/enstrophy diagnostics on grid fields.
%
%   vort = projection_vorticity_diagnostics(G, params)
%
% Uses omega_z = dUy/dx - dUx/dy. x is periodic and y is bounded. If a
% solid mask is provided in params.solidMask or params.solidMaskForProjection,
% solid cells are excluded from scalar summaries but kept as zero-velocity
% cells for finite differences. This deliberately retains wall vorticity near
% the cylinder boundary.
%
% Wake probes:
%   - center probe on the wake centerline,
%   - upper/lower off-center probes for shedding-frequency estimation.
% The off-center probes are important because centerline signals can contain a
% doubled frequency or a weak alternating component depending on symmetry.

Ux = double(G.Ux);
Uy = double(G.Uy);
[Nx, Ny] = size(Ux);
dx = G.dx;
dy = G.dy;

solidMask = false(Nx, Ny);
if isfield(params, 'solidMask') && ~isempty(params.solidMask)
    solidMask = logical(params.solidMask);
elseif isfield(params, 'solidMaskForProjection') && ~isempty(params.solidMaskForProjection)
    solidMask = logical(params.solidMaskForProjection);
end
if ~isequal(size(solidMask), size(Ux))
    solidMask = false(Nx, Ny);
end

Ux(solidMask) = 0;
Uy(solidMask) = 0;

% Periodic x derivative.
dUyDx = (Uy([2:Nx 1], :) - Uy([Nx 1:Nx-1], :)) / (2 * dx);

% Bounded y derivative: central in the interior, one-sided at walls.
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
fluidMask = ~solidMask;
wf = omega(fluidMask);

vort = struct();
vort.omega = omega;
vort.solidMask = solidMask;
vort.fluidMask = fluidMask;
vort.meanOmega = mean(wf, 'omitnan');
vort.rmsOmega = sqrt(mean(wf.^2, 'omitnan'));
vort.meanAbsOmega = mean(abs(wf), 'omitnan');
vort.maxAbsOmega = max(abs(wf));
vort.meanEnstrophy = 0.5 * mean(wf.^2, 'omitnan');
vort.totalEnstrophy = 0.5 * sum(wf.^2, 'omitnan') * dx * dy;

% Wake probes behind the cylinder. Distance and offset are specified in
% cylinder diameters, so the same settings remain meaningful when R changes.
Lx = get_param(params, 'Lx', G.Lx);
Ly = get_param(params, 'Ly', G.Ly);
R = get_param(params, 'cylinderRadius', 0.10 * Ly);
D = 2 * R;
cx = get_param(params, 'cylinderCenterX', 0.35 * Lx);
cy = get_param(params, 'cylinderCenterY', 0.50 * Ly);
probeDistance = get_param(params, 'wakeProbeDistanceDiameters', 2.0) * D;
probeOffset = get_param(params, 'wakeProbeYOffsetDiameters', 0.35) * D;

probeX = mod(cx + probeDistance, Lx);
center = sample_probe(probeX, cy, Ux, Uy, omega, solidMask, G);
upper = sample_probe(probeX, cy + probeOffset, Ux, Uy, omega, solidMask, G);
lower = sample_probe(probeX, cy - probeOffset, Ux, Uy, omega, solidMask, G);

vort.probeDistanceDiameters = get_param(params, 'wakeProbeDistanceDiameters', 2.0);
vort.probeYOffsetDiameters = get_param(params, 'wakeProbeYOffsetDiameters', 0.35);
vort.centerProbe = center;
vort.upperProbe = upper;
vort.lowerProbe = lower;

% Backward-compatible names: keep the centerline probe as wakeOmega/wakeUx/wakeUy.
vort.probeX = center.x;
vort.probeY = center.y;
vort.probeIx = center.ix;
vort.probeIy = center.iy;
vort.wakeOmega = center.omega;
vort.wakeUx = center.Ux;
vort.wakeUy = center.Uy;

% Preferred off-center shedding signal.
vort.wakeOmegaUpper = upper.omega;
vort.wakeUxUpper = upper.Ux;
vort.wakeUyUpper = upper.Uy;
vort.wakeOmegaLower = lower.omega;
vort.wakeUxLower = lower.Ux;
vort.wakeUyLower = lower.Uy;
vort.wakeOmegaAntiSym = 0.5 * (upper.omega - lower.omega);
vort.wakeUyAntiSym = 0.5 * (upper.Uy - lower.Uy);
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
