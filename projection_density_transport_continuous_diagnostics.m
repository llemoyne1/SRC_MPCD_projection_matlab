function tr = projection_density_transport_continuous_diagnostics(Gclassic, Gprojected, params)
%PROJECTION_DENSITY_TRANSPORT_CONTINUOUS_DIAGNOSTICS Continuous density transport proxy.
%
%   tr = projection_density_transport_continuous_diagnostics(Gclassic, Gprojected, params)
%
% Integer one-step population forecasts can be insensitive when dt is small,
% because few particles cross cell boundaries. This diagnostic evaluates the
% continuum proxy
%      deltaN = -dt div(N U)
% on the grid for the classic and projected fields.

if nargin < 3
    error('Usage: tr = projection_density_transport_continuous_diagnostics(Gclassic, Gprojected, params)');
end

N = double(Gclassic.N);
fluxClassicX = N .* Gclassic.Ux;
fluxClassicY = N .* Gclassic.Uy;
fluxProjectedX = N .* Gprojected.Ux;
fluxProjectedY = N .* Gprojected.Uy;

divClassic = divergence_periodic_x_bounded_y(fluxClassicX, fluxClassicY, params);
divProjected = divergence_periodic_x_bounded_y(fluxProjectedX, fluxProjectedY, params);

deltaClassic = -params.dt * divClassic;
deltaProjected = -params.dt * divProjected;
deltaDiff = deltaProjected - deltaClassic;

tr = struct();
tr.deltaClassic = deltaClassic;
tr.deltaProjected = deltaProjected;
tr.deltaProjectedMinusClassic = deltaDiff;
tr.divFluxClassic = divClassic;
tr.divFluxProjected = divProjected;
tr.classic = summarize(deltaClassic, N);
tr.projected = summarize(deltaProjected, N);
tr.projectedMinusClassic = summarize(deltaDiff, N);
tr.stdN = std(N(:));
tr.meanN = mean(N(:));
tr.nEmptyCells = nnz(N(:) == 0);
tr.massDeltaClassic = sum(deltaClassic(:));
tr.massDeltaProjected = sum(deltaProjected(:));
tr.massDeltaProjectedMinusClassic = sum(deltaDiff(:));
end

function div = divergence_periodic_x_bounded_y(Fx, Fy, params)
[Nx, Ny] = size(Fx);
dx = params.Lx / Nx;
dy = params.Ly / Ny;

ixp = [2:Nx, 1];
ixm = [Nx, 1:Nx-1];
dFx = (Fx(ixp, :) - Fx(ixm, :)) / (2*dx);

dFy = zeros(Nx, Ny);
if Ny > 1
    dFy(:, 1) = (Fy(:, 2) - Fy(:, 1)) / dy;
    dFy(:, Ny) = (Fy(:, Ny) - Fy(:, Ny-1)) / dy;
end
if Ny > 2
    dFy(:, 2:Ny-1) = (Fy(:, 3:Ny) - Fy(:, 1:Ny-2)) / (2*dy);
end

div = dFx + dFy;
end

function s = summarize(delta, N)
d = delta(:);
s = struct();
s.rms = sqrt(mean(d.^2));
s.maxAbs = max(abs(d));
s.meanAbs = mean(abs(d));
s.std = std(d);
s.sum = sum(d);
s.rmsOverGamma = s.rms / max(mean(N(:)), eps);
s.maxAbsOverGamma = s.maxAbs / max(mean(N(:)), eps);
s.predictedStdAfter = std(N(:) + d);
s.predictedStdDelta = s.predictedStdAfter - std(N(:));
end
