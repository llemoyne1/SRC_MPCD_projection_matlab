function [state, info] = projection_initialize_particles_step_channel(params)
%PROJECTION_INITIALIZE_PARTICLES_STEP_CHANNEL Initialize particles outside a rectangular step.
%
%   [state, info] = projection_initialize_particles_step_channel(params)
%
% Supports exact_per_cell-style initialization on the fluid cells only.

mode = lower(strrep(char(string(get_param(params, 'initialPopulationMode', 'exact_per_fluid_cell'))), '-', '_'));
if ~ismember(mode, {'exact_per_fluid_cell','exact_per_cell','uniform_by_cell','homogeneous_per_cell','cell_exact'})
    error('projection_initialize_particles_step_channel supports exact_per_fluid_cell/exact_per_cell only. Received %s.', mode);
end

gamma = get_param(params, 'gamma', 20);
gammaInt = round(gamma);
if abs(gamma - gammaInt) > 1e-12 || gammaInt < 0
    error('Step exact initialization requires non-negative integer gamma. Received %.16g.', gamma);
end

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;
solidMask = mpcd_step_mask(params);
nFluidCells = nnz(~solidMask);
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
            xp = (ix - 1 + rand()) * dx;
            yp = (iy - 1 + rand()) * dy;
            if point_inside_step(xp, yp, params)
                [xp, yp] = push_point_outside_step(xp, yp, params);
            end
            x(cursor, :) = [xp, yp];
        end
    end
end

kBT = get_param(params, 'kBT', 0.02);
v = sqrt(kBT) * randn(Np, 2);
if logical(get_param(params, 'initialVelocityZeroGlobalMean', true)) && Np > 0
    v(:, 1) = v(:, 1) - mean(v(:, 1), 'omitnan');
    v(:, 2) = v(:, 2) - mean(v(:, 2), 'omitnan');
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

function inside = point_inside_step(xp, yp, params)
Lx = get_param(params, 'Lx', 2.0);
Ly = get_param(params, 'Ly', 1.0);
x0 = get_param(params, 'stepX0', 0.15 * Lx);
x1 = get_param(params, 'stepX1', 0.35 * Lx);
h = get_param(params, 'stepHeight', 0.25 * Ly);
xp = mod(xp, Lx);
inside = xp >= x0 && xp <= x1 && yp >= 0 && yp <= h;
end

function [xp, yp] = push_point_outside_step(xp, yp, params)
Lx = get_param(params, 'Lx', 2.0);
Ly = get_param(params, 'Ly', 1.0);
x0 = get_param(params, 'stepX0', 0.15 * Lx);
x1 = get_param(params, 'stepX1', 0.35 * Lx);
h = get_param(params, 'stepHeight', 0.25 * Ly);
epsWall = get_param(params, 'stepWallEpsilon', 1e-10 * max(Lx, Ly));
xp = mod(xp, Lx);
if ~(xp >= x0 && xp <= x1 && yp >= 0 && yp <= h)
    return;
end
[dmin, id] = min([xp - x0, x1 - xp, h - yp]); %#ok<ASGLU>
switch id
    case 1
        xp = x0 - epsWall;
    case 2
        xp = x1 + epsWall;
    otherwise
        yp = h + epsWall;
end
xp = mod(xp, Lx);
yp = min(max(yp, 0), Ly - eps(Ly));
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
