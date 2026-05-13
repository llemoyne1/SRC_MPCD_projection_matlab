function mask = mpcd_cylinder_mask(params)
%MPCD_CYLINDER_MASK Cell-center solid mask for a circular obstacle.
%
%   mask = mpcd_cylinder_mask(params)
%
% Returns an Nx-by-Ny logical mask.  mask(ix,iy)=true means the cell center
% lies inside the fixed cylinder.  The x direction is periodic, matching the
% channel projection convention used in the Q9 prototypes.

required = {'Lx','Ly','Nx','Ny'};
for k = 1:numel(required)
    if ~isfield(params, required{k})
        error('params.%s is required.', required{k});
    end
end

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
R = get_param(params, 'cylinderRadius', 0.10 * Ly);
cx = get_param(params, 'cylinderCenterX', 0.35 * Lx);
cy = get_param(params, 'cylinderCenterY', 0.50 * Ly);

xCenters = ((0:Nx-1).' + 0.5) * Lx / Nx;
yCenters = ((0:Ny-1) + 0.5) * Ly / Ny;
[X, Y] = ndgrid(xCenters, yCenters);
DX = X - cx;
DX = DX - Lx * round(DX / Lx);
DY = Y - cy;
mask = (DX.^2 + DY.^2) <= R^2;
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
