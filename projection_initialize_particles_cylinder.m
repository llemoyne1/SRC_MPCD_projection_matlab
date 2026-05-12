function [state, info] = projection_initialize_particles_cylinder(params)
%PROJECTION_INITIALIZE_PARTICLES_CYLINDER Initialize particles outside a cylinder.
%
%   [state, info] = projection_initialize_particles_cylinder(params)
%
% This initializer mirrors initialPopulationMode='exact_per_cell', but only
% for fluid cells.  Cell centers inside the circular obstacle are skipped;
% boundary cells are filled by rejection sampling so particles are not created
% inside the true circle.

if nargin < 1 || isempty(params)
    error('params is required.');
end

mode = lower(strrep(char(string(get_param(params, 'initialPopulationMode', 'exact_per_cell'))), '-', '_'));
if ~ismember(mode, {'exact_per_cell','uniform_by_cell','homogeneous_per_cell','cell_exact'})
    error('projection_initialize_particles_cylinder currently supports exact_per_cell only. Received %s.', mode);
end

gamma = get_param(params, 'gamma', 20);
gammaInt = round(gamma);
if abs(gamma - gammaInt) > 1e-12 || gammaInt < 0
    error('exact_per_cell cylinder initialization requires non-negative integer gamma. Received %.16g.', gamma);
end

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;
solidMask = mpcd_cylinder_mask(params);
fluidCells = find(~solidMask);
nFluidCells = numel(fluidCells);
Np = gammaInt * nFluidCells;

x = zeros(Np, 2);
cursor = 0;
for ix = 1:Nx
    for iy = 1:Ny
        if solidMask(ix, iy)
            continue;
        end
        for p = 1:gammaInt
            cursor = cursor + 1;
            accepted = false;
            for attempt = 1:200
                xp = (ix - 1 + rand()) * dx;
                yp = (iy - 1 + rand()) * dy;
                if ~point_inside_cylinder(xp, yp, params)
                    x(cursor, :) = [xp, yp];
                    accepted = true;
                    break;
                end
            end
            if ~accepted
                % Fallback for highly cut boundary cells: place near the cell
                % center and push out along the cylinder normal if needed.
                xp = (ix - 0.5) * dx;
                yp = (iy - 0.5) * dy;
                [xp, yp] = push_point_outside_cylinder(xp, yp, params);
                x(cursor, :) = [xp, yp];
            end
        end
    end
end

kBT = get_param(params, 'kBT', 1.0);
v = sqrt(kBT) * randn(Np, 2);
if logical(get_param(params, 'initialVelocityZeroGlobalMean', true)) && Np > 0
    v(:, 1) = v(:, 1) - mean(v(:, 1));
    v(:, 2) = v(:, 2) - mean(v(:, 2));
end
v(:, 1) = v(:, 1) + get_param(params, 'initialMeanVelocityX', 0.0);
v(:, 2) = v(:, 2) + get_param(params, 'initialMeanVelocityY', 0.0);

state = struct();
state.x = x;
state.v = v;

pop = projection_population_diagnostics(state.x, params, 'periodicX', true, 'periodicY', false);
info = struct();
info.mode = 'exact_per_fluid_cell';
info.Np = Np;
info.gamma = gamma;
info.Nx = Nx;
info.Ny = Ny;
info.nFluidCells = nFluidCells;
info.nSolidCells = nnz(solidMask);
info.solidMask = solidMask;
info.population = pop;
info.initialMeanNFluid = mean(double(pop.N(~solidMask)), 'omitnan');
info.initialStdNFluid = std(double(pop.N(~solidMask)), 'omitnan');
info.initialEmptyFluidCells = nnz(pop.N(~solidMask) == 0);
info.initialSolidParticles = nnz(pop.N(solidMask));
end

function inside = point_inside_cylinder(xp, yp, params)
Lx = params.Lx;
R = get_param(params, 'cylinderRadius', 0.10 * params.Ly);
cx = get_param(params, 'cylinderCenterX', 0.35 * params.Lx);
cy = get_param(params, 'cylinderCenterY', 0.50 * params.Ly);
dx = xp - cx;
dx = dx - Lx * round(dx / Lx);
dy = yp - cy;
inside = (dx*dx + dy*dy) <= R*R;
end

function [xp, yp] = push_point_outside_cylinder(xp, yp, params)
Lx = params.Lx;
Ly = params.Ly;
R = get_param(params, 'cylinderRadius', 0.10 * Ly);
cx = get_param(params, 'cylinderCenterX', 0.35 * params.Lx);
cy = get_param(params, 'cylinderCenterY', 0.50 * Ly);
epsWall = 1e-10 * max(Ly, R);
dx = xp - cx;
dx = dx - Lx * round(dx / Lx);
dy = yp - cy;
r = hypot(dx, dy);
if r < R + epsWall
    if r < epsWall
        nx = 1.0;
        ny = 0.0;
    else
        nx = dx / r;
        ny = dy / r;
    end
    cxImage = xp - dx;
    xp = mod(cxImage + (R + epsWall) * nx, Lx);
    yp = cy + (R + epsWall) * ny;
    yp = min(max(yp, 0), Ly - eps(Ly));
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
