function [stateOut, diag] = resamp_step_classic_periodic_weighted(state, params)
%RESAMP_STEP_CLASSIC_PERIODIC_WEIGHTED Weighted 2D SRD/MPCD step in a periodic box.
%
% Supports both compact states and preallocated pool states.  Inactive pool
% slots are ignored by streaming, forcing, collision and diagnostics.

validate_state(state);
Lx = get_param(params, 'Lx', []);
Ly = get_param(params, 'Ly', []);
Nx = get_param(params, 'Nx', []);
Ny = get_param(params, 'Ny', []);
dt = get_param(params, 'dt', []);
alphaDeg = get_param(params, 'alphaDeg', 90);
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);
useRandomGridShift = logical(get_param(params, 'useRandomGridShift', true));
if isempty(Lx) || isempty(Ly) || isempty(Nx) || isempty(Ny) || isempty(dt)
    error('params must contain Lx, Ly, Nx, Ny and dt.');
end

activeMask = resamp_active_mask(state);
idx = find(activeMask);
stateOut = state;
if isempty(idx)
    diag = empty_diag(state, params);
    return;
end

x = state.x(idx, :);
v = state.v(idx, :);
m = state.m(idx);
Np = numel(idx);
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

% Acceleration kicks are independent of particle mass.
v(:,1) = v(:,1) + dt * bodyForceX;
v(:,2) = v(:,2) + dt * bodyForceY;
[v, forceInfo] = projection_apply_taylor_green_forcing(x, v, params);

% Periodic streaming.
x(:,1) = mod(x(:,1) + dt * v(:,1), Lx);
x(:,2) = mod(x(:,2) + dt * v(:,2), Ly);

% Random shifted SRD collision grid.
if useRandomGridShift
    shiftX = (rand() - 0.5) * dx;
    shiftY = (rand() - 0.5) * dy;
else
    shiftX = 0.0;
    shiftY = 0.0;
end
xs = mod(x(:,1) + shiftX, Lx);
ys = mod(x(:,2) + shiftY, Ly);
ix = floor(xs / dx) + 1;
iy = floor(ys / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
cellId = iy + Ny * (ix - 1);

Ncell = accumarray(cellId, 1, [Nc, 1], @sum, 0);
Mcell = accumarray(cellId, m, [Nc, 1], @sum, 0);
Px = accumarray(cellId, m .* v(:,1), [Nc, 1], @sum, 0);
Py = accumarray(cellId, m .* v(:,2), [Nc, 1], @sum, 0);
Ux = zeros(Nc, 1);
Uy = zeros(Nc, 1);
occ = Mcell > eps;
Ux(occ) = Px(occ) ./ Mcell(occ);
Uy(occ) = Py(occ) ./ Mcell(occ);

momentumBeforeCollision = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];

alpha = alphaDeg * pi / 180.0;
signs = 2.0 * (rand(Nc, 1) > 0.5) - 1.0;
angles = signs * alpha;
ca = cos(angles(cellId));
sa = sin(angles(cellId));
uxp = Ux(cellId);
uyp = Uy(cellId);
rvx = v(:,1) - uxp;
rvy = v(:,2) - uyp;
v(:,1) = uxp + ca .* rvx - sa .* rvy;
v(:,2) = uyp + sa .* rvx + ca .* rvy;

momentumAfterCollision = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];

stateOut.x(idx,:) = x;
stateOut.v(idx,:) = v;
stateOut.m(idx) = m;
stateOut.Nactive = nnz(resamp_active_mask(stateOut));
stateOut.Ncapacity = size(stateOut.x, 1);

thermostatAfterStep = logical(get_param(params, 'thermostatAfterStep', ...
    get_param(params, 'thermostatAfterProjection', false)));
if thermostatAfterStep
    [stateOut.v, thermostatInfo] = resamp_apply_cell_thermostat_weighted(stateOut.x, stateOut.v, stateOut.m, params, ...
        'periodicX', true, 'periodicY', true, 'activeMask', resamp_active_mask(stateOut));
else
    thermostatInfo = empty_thermostat_info();
end

md = resamp_population_mass_diagnostics(stateOut, params, 'periodicX', true, 'periodicY', true);
pool = resamp_particle_pool_info(stateOut);

diag = struct();
diag.Np = Np;
diag.NpActive = pool.Nactive;
diag.Ncapacity = pool.Ncapacity;
diag.Nfree = pool.Nfree;
diag.Nx = Nx;
diag.Ny = Ny;
diag.dx = dx;
diag.dy = dy;
diag.shiftX = shiftX;
diag.shiftY = shiftY;
diag.alphaDeg = alphaDeg;
diag.bodyForceX = bodyForceX;
diag.bodyForceY = bodyForceY;
diag.taylorGreenForce = forceInfo;
diag.thermostatAfterStep = thermostatAfterStep;
diag.thermostat = thermostatInfo;
diag.NMean = md.NMean;
diag.NStd = md.NStd;
diag.NMin = md.NMin;
diag.NMax = md.NMax;
diag.MMean = md.MMean;
diag.MStd = md.MStd;
diag.MRelRms = md.MRelRms;
diag.mParticleMean = md.mParticleMean;
diag.mParticleStd = md.mParticleStd;
diag.nEmptyCells = md.nEmptyCells;
diag.meanVxWeighted = md.meanVxWeighted;
diag.meanVyWeighted = md.meanVyWeighted;
diag.kBTWeighted = md.kBTWeighted;
diag.kineticEnergyMeanWeighted = md.kineticEnergyMeanWeighted;
diag.totalMass = md.totalMass;
diag.totalMomentum = md.totalMomentum;
diag.momentumBeforeCollision = momentumBeforeCollision;
diag.momentumAfterCollision = momentumAfterCollision;
diag.collisionDeltaPNorm = norm(momentumAfterCollision - momentumBeforeCollision);
diag.massDiagnostics = md;
end

function diag = empty_diag(state, params)
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);
pool = resamp_particle_pool_info(state);
diag = struct('Np', 0, 'NpActive', pool.Nactive, 'Ncapacity', pool.Ncapacity, ...
    'Nfree', pool.Nfree, 'Nx', params.Nx, 'Ny', params.Ny, 'dx', params.Lx/params.Nx, ...
    'dy', params.Ly/params.Ny, 'shiftX', 0, 'shiftY', 0, 'alphaDeg', get_param(params,'alphaDeg',90), ...
    'bodyForceX', get_param(params,'bodyForceX',0), 'bodyForceY', get_param(params,'bodyForceY',0), ...
    'taylorGreenForce', struct(), 'thermostatAfterStep', false, 'thermostat', empty_thermostat_info(), ...
    'NMean', md.NMean, 'NStd', md.NStd, 'NMin', md.NMin, 'NMax', md.NMax, ...
    'MMean', md.MMean, 'MStd', md.MStd, 'MRelRms', md.MRelRms, ...
    'mParticleMean', md.mParticleMean, 'mParticleStd', md.mParticleStd, ...
    'nEmptyCells', md.nEmptyCells, 'meanVxWeighted', md.meanVxWeighted, 'meanVyWeighted', md.meanVyWeighted, ...
    'kBTWeighted', md.kBTWeighted, 'kineticEnergyMeanWeighted', md.kineticEnergyMeanWeighted, ...
    'totalMass', md.totalMass, 'totalMomentum', md.totalMomentum, ...
    'momentumBeforeCollision', [0 0], 'momentumAfterCollision', [0 0], 'collisionDeltaPNorm', 0, ...
    'massDiagnostics', md);
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must be a struct with fields x, v and m.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
if numel(state.m) ~= size(state.x,1)
    error('state.m must have one entry per particle slot.');
end
if any(~isfinite(state.m(:))) || any(state.m(:) < 0)
    error('state.m must contain finite non-negative masses.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end

function info = empty_thermostat_info()
info = struct();
info.enabled = false;
info.kBTTarget = NaN;
info.kBTBefore = NaN;
info.kBTAfter = NaN;
info.meanScale = NaN;
info.maxScale = NaN;
end
