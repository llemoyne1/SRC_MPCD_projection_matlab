function [stateOut, diag] = mpcd_step_classic_poiseuille(state, params)
%MPCD_STEP_CLASSIC_POISEUILLE Classic SRC/MPCD step in an x-periodic channel.
%
%   [stateOut, diag] = mpcd_step_classic_poiseuille(state, params)
%
% Order:
%   body-force kick -> x-periodic streaming -> y-wall collision -> SRD collision
%   -> optional common post-step cell thermostat.
%
% The post-step thermostat is the corrected/vectorized cell thermostat used
% in the current Taylor--Green validation.  It rescales only velocity
% fluctuations around the local cell mean and therefore preserves local cell
% momentum up to roundoff.

validate_state(state);

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dt = params.dt;
alphaDeg = get_param(params, 'alphaDeg', 90);
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);
useRandomGridShiftX = logical(get_param(params, 'useRandomGridShiftX', true));
useRandomGridShiftY = logical(get_param(params, 'useRandomGridShiftY', false));

if Lx <= 0 || Ly <= 0 || Nx <= 0 || Ny <= 0 || dt <= 0
    error('Lx, Ly, Nx, Ny and dt must be strictly positive.');
end

x = state.x;
v = state.v;
Np = size(x, 1);
dx = Lx / Nx;
dy = Ly / Ny;
Nc = Nx * Ny;

% Force kick and streaming.
v(:, 1) = v(:, 1) + dt * bodyForceX;
v(:, 2) = v(:, 2) + dt * bodyForceY;
x(:, 1) = mod(x(:, 1) + dt * v(:, 1), Lx);
x(:, 2) = x(:, 2) + dt * v(:, 2);
[x, v, wallInfo] = mpcd_apply_wall_bc_y(x, v, params);

% SRD collision. Optional virtual wall particles contribute aggregate
% collision mass/momentum in the near-wall cells and are discarded after the
% collision; only real particles are updated.
[v, collisionInfo] = mpcd_srd_collision_channel_virtual_walls(x, v, params);

thermostatAfterStep = logical(get_param(params, 'thermostatAfterStep', false));
if thermostatAfterStep
    [v, thermostatInfo] = projection_apply_cell_thermostat(x, v, params, ...
        'periodicX', true, 'periodicY', false);
else
    thermostatInfo = empty_thermostat_info();
end

stateOut = state;
stateOut.x = x;
stateOut.v = v;

thermal = projection_thermal_diagnostics(x, v, params, 'periodicX', true, 'periodicY', false, 'minCount', 1);

diag = struct();
diag.Np = Np;
diag.Nx = Nx;
diag.Ny = Ny;
diag.dx = dx;
diag.dy = dy;
diag.shiftX = collisionInfo.shiftX;
diag.shiftY = collisionInfo.shiftY;
diag.alphaDeg = alphaDeg;
diag.bodyForceX = bodyForceX;
diag.bodyForceY = bodyForceY;
diag.wallInfo = wallInfo;
diag.thermostatAfterStep = thermostatAfterStep;
diag.thermostat = thermostatInfo;
diag.thermal = thermal;
diag.NMean = collisionInfo.NMeanReal;
diag.NStd = collisionInfo.NStdReal;
diag.NMin = collisionInfo.NMinReal;
diag.NMax = collisionInfo.NMaxReal;
diag.nEmptyCells = collisionInfo.nEmptyRealCells;
diag.NMeanTotalCollision = collisionInfo.NMeanTotal;
diag.NStdTotalCollision = collisionInfo.NStdTotal;
diag.wallVirtualParticles = collisionInfo.wallVirtualParticles;
diag.nVirtualWallCells = collisionInfo.nVirtualCells;
diag.nVirtualWallParticles = collisionInfo.nVirtualParticlesTotal;
diag.meanVx = mean(v(:, 1), 'omitnan');
diag.meanVy = mean(v(:, 2), 'omitnan');
diag.kBT = estimate_kBT(v);
diag.kBTCell = thermal.kBTCellRelative;
diag.kineticEnergyMean = 0.5 * mean(sum(v.^2, 2), 'omitnan');
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v')
    error('state must be a struct with fields x and v.');
end
if size(state.x, 2) ~= 2 || size(state.v, 2) ~= 2 || size(state.x, 1) ~= size(state.v, 1)
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
info = struct('enabled', false, 'targetKBT', NaN, 'strength', NaN, ...
    'minParticlesPerCell', NaN, 'maxScale', NaN, 'nThermostattedCells', 0, ...
    'meanScale', NaN, 'minScale', NaN, 'maxScaleApplied', NaN, ...
    'meanKBTBefore', NaN, 'meanKBTAfter', NaN, 'rmsVelocityChange', 0.0);
end

function kBT = estimate_kBT(v)
u = mean(v, 1, 'omitnan');
c = v - u;
kBT = 0.5 * mean(sum(c.^2, 2), 'omitnan');
end
