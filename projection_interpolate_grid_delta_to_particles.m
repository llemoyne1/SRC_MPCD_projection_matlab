function dv = projection_interpolate_grid_delta_to_particles(x, dUx, dUy, params, varargin)
%PROJECTION_INTERPOLATE_GRID_DELTA_TO_PARTICLES Interpolate grid correction to particles.
%
%   dv = projection_interpolate_grid_delta_to_particles(x, dUx, dUy, params)
%
% Interpolation from cell-centered grid values to particles. This function is
% intentionally independent from the legacy simulator so it can be tested
% before coupling the pressure projection to MPCD.
%
% Optional name-value arguments:
%   'periodicX'   default true
%   'periodicY'   default false
%   'method'      'bilinear' or 'nearest', default 'bilinear'
%
% Grid convention:
%   dUx and dUy are Nx-by-Ny arrays.
%   Cell center (ix,iy) is located at:
%       x = (ix - 0.5) * dx, y = (iy - 0.5) * dy.

if ~isequal(size(dUx), size(dUy))
    error('dUx and dUy must have the same size.');
end
[Nx, Ny] = size(dUx);
if Nx ~= params.Nx || Ny ~= params.Ny
    error('Grid size mismatch: got %dx%d, params has %dx%d.', Nx, Ny, params.Nx, params.Ny);
end
if size(x, 2) < 2
    error('x must be an Np-by-2 array.');
end

periodicX = true;
periodicY = false;
method = "bilinear";
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "method"
            method = lower(string(val));
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Lx = params.Lx;
Ly = params.Ly;
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

switch method
    case "nearest"
        ix = floor(xp / dx) + 1;
        iy = floor(yp / dy) + 1;
        ix = min(max(ix, 1), Nx);
        iy = min(max(iy, 1), Ny);
        ind = sub2ind([Nx, Ny], ix, iy);
        dv = zeros(size(x, 1), 2);
        dv(:, 1) = dUx(ind);
        dv(:, 2) = dUy(ind);
        return;

    case "bilinear"
        % Continuous 1-based cell-center coordinate. A particle located at
        % the center of cell 1 has c=1. Values near the left periodic edge
        % have c<1 and correctly interpolate between cells Nx and 1.
        cx = xp / dx + 0.5;
        cy = yp / dy + 0.5;

        ix0 = floor(cx);
        iy0 = floor(cy);
        wx = cx - ix0;
        wy = cy - iy0;

        ix1 = ix0 + 1;
        iy1 = iy0 + 1;

        ix0 = fix_index(ix0, Nx, periodicX);
        ix1 = fix_index(ix1, Nx, periodicX);
        iy0 = fix_index(iy0, Ny, periodicY);
        iy1 = fix_index(iy1, Ny, periodicY);

    otherwise
        error('Unknown interpolation method: %s', method);
end

w00 = (1 - wx) .* (1 - wy);
w10 = wx .* (1 - wy);
w01 = (1 - wx) .* wy;
w11 = wx .* wy;

ind00 = sub2ind([Nx, Ny], ix0, iy0);
ind10 = sub2ind([Nx, Ny], ix1, iy0);
ind01 = sub2ind([Nx, Ny], ix0, iy1);
ind11 = sub2ind([Nx, Ny], ix1, iy1);

dv = zeros(size(x, 1), 2);
dv(:, 1) = w00.*dUx(ind00) + w10.*dUx(ind10) + w01.*dUx(ind01) + w11.*dUx(ind11);
dv(:, 2) = w00.*dUy(ind00) + w10.*dUy(ind10) + w01.*dUy(ind01) + w11.*dUy(ind11);
end

function idx = fix_index(idx, N, periodic)
if periodic
    idx = mod(idx - 1, N) + 1;
else
    idx = min(max(idx, 1), N);
end
end
