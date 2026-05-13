function solidMask = mpcd_step_mask(params)
%MPCD_STEP_MASK Rectangular bottom-attached step mask on the Eulerian grid.
%
%   solidMask = mpcd_step_mask(params)
%
% The mask is true in solid cells.  Cell centers are used, so choosing
% step dimensions aligned with dx/dy gives a strictly grid-aligned obstacle.

Lx = get_param(params, 'Lx', 2.0);
Ly = get_param(params, 'Ly', 1.0);
Nx = get_param(params, 'Nx', 48);
Ny = get_param(params, 'Ny', 24);
x0 = get_param(params, 'stepX0', 0.15 * Lx);
x1 = get_param(params, 'stepX1', 0.35 * Lx);
h  = get_param(params, 'stepHeight', 0.25 * Ly);

if x1 <= x0
    error('stepX1 must be greater than stepX0.');
end
if h <= 0 || h >= Ly
    error('stepHeight must be in (0,Ly).');
end

dx = Lx / Nx;
dy = Ly / Ny;
xc = ((1:Nx) - 0.5) * dx;
yc = ((1:Ny) - 0.5) * dy;
[Xc, Yc] = ndgrid(xc, yc);
solidMask = (Xc >= x0) & (Xc <= x1) & (Yc <= h);
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
