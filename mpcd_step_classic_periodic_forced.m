function [stateOut, diag] = mpcd_step_classic_periodic_forced(state, params)
%MPCD_STEP_CLASSIC_PERIODIC_FORCED Classic 2D SRD/MPCD step in a periodic box with optional TG forcing.
%
% Order:
%   uniform body-force kick + TG forcing kick -> periodic streaming -> SRD collision
%   -> optional common post-step cell thermostat.

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

x = state.x;
v = state.v;
Np = size(x, 1);
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

% Uniform and modal body-force kicks.
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
Px = accumarray(cellId, v(:,1), [Nc, 1], @sum, 0);
Py = accumarray(cellId, v(:,2), [Nc, 1], @sum, 0);
Ux = zeros(Nc, 1);
Uy = zeros(Nc, 1);
occ = Ncell > 0;
Ux(occ) = Px(occ) ./ Ncell(occ);
Uy(occ) = Py(occ) ./ Ncell(occ);

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

% Optional common post-step thermostat.  This is intentionally available
% for the classic method too, so Taylor--Green classic/Q6/Q9 comparisons
% use the same thermal closure.  The cell thermostat rescales only relative
% velocities inside each cell and therefore preserves the cell mean velocity.
thermostatAfterStep = logical(get_param(params, 'thermostatAfterStep', ...
    get_param(params, 'thermostatAfterProjection', false)));
if thermostatAfterStep
    [v, thermostatInfo] = projection_apply_cell_thermostat(x, v, params, ...
        'periodicX', true, 'periodicY', true);
    thermCheck = projection_thermal_diagnostics(x, v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);

% fprintf(['THERMO classic: target=%.6g, strength=%.3g, cells=%d, ', ...
%          'info before=%.6g, info after=%.6g, diag cellRel=%.6g, diag localMean=%.6g, ', ...
%          'scale mean/min/max=%.6g/%.6g/%.6g\n'], ...
%     thermostatInfo.targetKBT, thermostatInfo.strength, ...
%     thermostatInfo.nThermostattedCells, ...
%     thermostatInfo.meanKBTBefore, thermostatInfo.meanKBTAfter, ...
%     thermCheck.kBTCellRelative, thermCheck.localKBTMean, ...
%     thermostatInfo.meanScale, thermostatInfo.minScale, thermostatInfo.maxScaleApplied);
else
    thermostatInfo = empty_thermostat_info();
end

stateOut = struct('x', x, 'v', v);

diag = struct();
diag.Np = Np;
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
diag.NMean = mean(Ncell, 'omitnan');
diag.NStd = std(double(Ncell), 0, 'omitnan');
diag.NMin = min(Ncell);
diag.NMax = max(Ncell);
diag.nEmptyCells = nnz(Ncell == 0);
diag.meanVx = mean(v(:,1), 'omitnan');
diag.meanVy = mean(v(:,2), 'omitnan');
diag.kBT = estimate_kBT(v);
diag.kineticEnergyMean = 0.5 * mean(sum(v.^2, 2), 'omitnan');
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v')
    error('state must be a struct with fields x and v.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
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

function kBT = estimate_kBT(v)
u = mean(v, 1, 'omitnan');
c = v - u;
kBT = 0.5 * mean(sum(c.^2, 2), 'omitnan');
end
