function [stateOut, diag] = mpcd_step_classic_piston(state, params)
%MPCD_STEP_CLASSIC_PISTON Classic SRC/MPCD step in a moving-piston channel.
%
%   [stateOut, diag] = mpcd_step_classic_piston(state, params)
%
% Clean piston kernel for the wallVP-v2 branch.  The active domain is
% [0,Lx] x [0,yTop], with fixed Nx-by-Ny grid and dy=yTop/Ny.  The top wall
% is allowed to move normally with velocity Up.  Virtual wall particles are
% handled by mpcd_srd_collision_channel_virtual_walls using the active height
% and the wall velocities [Ubottom,Vbottom] and [Utop,Vtop].
%
% Order:
%   optional affine compression -> force kick -> streaming -> moving-wall BC
%   -> SRC collision with optional wall virtual particles -> optional cell
%   thermostat.
%
% Only real particles are stored/advected.  Virtual wall particles contribute
% to the collision mean and then disappear.

validate_state(state);

Lx = params.Lx;
Ly0 = get_param(params, 'Ly0', get_param(params, 'LyReference', params.Ly));
Nx = params.Nx;
Ny = params.Ny;
dt = params.dt;
alphaDeg = get_param(params, 'alphaDeg', 90);
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);

if Lx <= 0 || Ly0 <= 0 || Nx <= 0 || Ny <= 0 || dt <= 0
    error('Lx, Ly0, Nx, Ny and dt must be strictly positive.');
end

yPrev = get_param(params, 'pistonYPrev', Ly0);
yTop = get_param(params, 'pistonYCurrent', Ly0);
Up = get_param(params, 'pistonVyCurrent', get_param(params, 'pistonVy', 0.0));
yTop = min(max(yTop, eps(Ly0)), Ly0);
yPrev = min(max(yPrev, eps(Ly0)), Ly0);

x = state.x;
v = state.v;
Np = size(x, 1);

% Optional quasi-static geometric compression.  This preserves a homogeneous
% initial density field during slow piston motion and avoids dumping the whole
% displacement into the first layer below the moving wall.
if logical(get_param(params, 'pistonAffineReposition', true)) && abs(yTop - yPrev) > eps(Ly0)
    scale = yTop / max(yPrev, eps(Ly0));
    x(:, 2) = min(max(x(:, 2), 0), max(yPrev - eps(yPrev), 0));
    x(:, 2) = min(max(x(:, 2) * scale, 0), max(yTop - eps(yTop), 0));
end

% Force kick and streaming in the active domain.
v(:, 1) = v(:, 1) + dt * bodyForceX;
v(:, 2) = v(:, 2) + dt * bodyForceY;
x(:, 1) = mod(x(:, 1) + dt * v(:, 1), Lx);
x(:, 2) = x(:, 2) + dt * v(:, 2);
[x, v, wallInfo] = apply_piston_wall_bc_y(x, v, params, yTop, Up);

% Collision with active-height parameters.  The top virtual wall velocity is
% the piston velocity Up in the normal direction.
collParams = params;
collParams.Ly = yTop;
collParams.Ly0 = Ly0;
collParams.Ubottom = get_param(params, 'Ubottom', 0.0);
collParams.Utop = get_param(params, 'Utop', 0.0);
collParams.Vbottom = get_param(params, 'Vbottom', 0.0);
collParams.Vtop = Up;
[v, collisionInfo] = mpcd_srd_collision_channel_virtual_walls(x, v, collParams);

% Mechanical diagnostics: impact impulses from the geometric wall BC plus
% momentum exchanged through virtual wall particles during SRC collision.
% pressureTopWall is positive when the fluid pushes upward on the top wall;
% with a downward piston velocity Up<0, powerOnFluid = -F_wall*Up is positive.
wallInfo.pressureTopWallImpact = wallInfo.pressureTopWall;
wallInfo.pressureBottomWallImpact = wallInfo.pressureBottomWall;
wallInfo.pressureTopWallVP = getf(collisionInfo, 'pressureTopWallVP', 0.0);
wallInfo.pressureBottomWallVP = getf(collisionInfo, 'pressureBottomWallVP', 0.0);
wallInfo.pressureTopWallTotal = wallInfo.pressureTopWallImpact + wallInfo.pressureTopWallVP;
wallInfo.pressureBottomWallTotal = wallInfo.pressureBottomWallImpact + wallInfo.pressureBottomWallVP;
wallInfo.impulseOnTopWallVPY = getf(collisionInfo, 'impulseTopWallVPy', 0.0);
wallInfo.impulseOnBottomWallVPY = getf(collisionInfo, 'impulseBottomWallVPy', 0.0);
wallInfo.impulseOnTopWallTotalY = wallInfo.impulseOnTopWallY + wallInfo.impulseOnTopWallVPY;
wallInfo.impulseOnBottomWallTotalY = wallInfo.impulseOnBottomWallY + wallInfo.impulseOnBottomWallVPY;
wallInfo.topWallForceOnFluidY = -wallInfo.pressureTopWallTotal * Lx;
wallInfo.bottomWallForceOnFluidY = -wallInfo.pressureBottomWallTotal * Lx;
wallInfo.pistonPowerOnFluid = wallInfo.topWallForceOnFluidY * Up;
prevWork = 0.0;
if isfield(state, 'piston') && isfield(state.piston, 'workOnFluidCumulative') && ~isempty(state.piston.workOnFluidCumulative)
    prevWork = state.piston.workOnFluidCumulative;
end
wallInfo.pistonWorkIncrement = wallInfo.pistonPowerOnFluid * dt;
wallInfo.pistonWorkOnFluidCumulative = prevWork + wallInfo.pistonWorkIncrement;

thermostatAfterStep = logical(get_param(params, 'thermostatAfterStep', false));
if thermostatAfterStep
    thermoParams = collParams;
    [v, thermostatInfo] = projection_apply_cell_thermostat(x, v, thermoParams, ...
        'periodicX', true, 'periodicY', false);
else
    thermostatInfo = empty_thermostat_info();
end

stateOut = state;
stateOut.x = x;
stateOut.v = v;
stateOut.piston = struct('yTop', yTop, 'yPrev', yPrev, 'Up', Up, ...
    'Ly0', Ly0, 'compression', 1.0 - yTop / Ly0, 'activeHeight', yTop, ...
    'workOnFluidCumulative', wallInfo.pistonWorkOnFluidCumulative, ...
    'lastWorkIncrement', wallInfo.pistonWorkIncrement, ...
    'lastPowerOnFluid', wallInfo.pistonPowerOnFluid);

thermal = projection_thermal_diagnostics(x, v, collParams, ...
    'periodicX', true, 'periodicY', false, 'minCount', 1);

pop = projection_population_diagnostics(x, collParams, ...
    'periodicX', true, 'periodicY', false, 'targetGamma', get_param(params, 'gamma', NaN));

diag = struct();
diag.Np = Np;
diag.Nx = Nx;
diag.Ny = Ny;
diag.Lx = Lx;
diag.Ly = yTop;
diag.Ly0 = Ly0;
diag.dx = Lx / Nx;
diag.dy = yTop / Ny;
diag.alphaDeg = alphaDeg;
diag.bodyForceX = bodyForceX;
diag.bodyForceY = bodyForceY;
diag.wallInfo = wallInfo;
diag.piston = stateOut.piston;
diag.collisionInfo = collisionInfo;
diag.pressureTopWallImpact = wallInfo.pressureTopWallImpact;
diag.pressureTopWallVP = wallInfo.pressureTopWallVP;
diag.pressureTopWallTotal = wallInfo.pressureTopWallTotal;
diag.pressureBottomWallImpact = wallInfo.pressureBottomWallImpact;
diag.pressureBottomWallVP = wallInfo.pressureBottomWallVP;
diag.pressureBottomWallTotal = wallInfo.pressureBottomWallTotal;
diag.pistonPowerOnFluid = wallInfo.pistonPowerOnFluid;
diag.pistonWorkIncrement = wallInfo.pistonWorkIncrement;
diag.pistonWorkOnFluidCumulative = wallInfo.pistonWorkOnFluidCumulative;
diag.thermostatAfterStep = thermostatAfterStep;
diag.thermostat = thermostatInfo;
diag.thermal = thermal;
diag.population = pop;
diag.NMean = pop.meanN;
diag.NStd = pop.stdN;
diag.NMin = pop.minN;
diag.NMax = pop.maxN;
diag.nEmptyCells = pop.nEmptyCells;
diag.shiftX = collisionInfo.shiftX;
diag.shiftY = collisionInfo.shiftY;
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

function [x, v, info] = apply_piston_wall_bc_y(x, v, params, yTop, Up)
modeBottom = lower(char(string(get_param(params, 'wallModeY', 'bounceback'))));
modeTopTangential = lower(char(string(get_param(params, 'pistonTangentialMode', get_param(params, 'wallModeY', 'bounceback')))));
Ubottom = get_param(params, 'Ubottom', 0.0);
Utop = get_param(params, 'Utop', 0.0);
wallSigma = get_param(params, 'wallSigma', sqrt(max(get_param(params, 'kBT', 1.0), 0)));

info = struct('nBot', 0, 'nTop', 0, 'dEwall', 0, ...
              'dPxBot', 0, 'dPxTop', 0, 'dPyBot', 0, 'dPyTop', 0, ...
              'dPxBotSigned', 0, 'dPxTopSigned', 0, ...
              'dPyBotSigned', 0, 'dPyTopSigned', 0, ...
              'impulseOnBottomWallY', 0, 'impulseOnTopWallY', 0, ...
              'pressureBottomWall', NaN, 'pressureTopWall', NaN, ...
              'yTop', yTop, 'Up', Up);

for pass = 1:10
    hitB = x(:, 2) < 0;
    if any(hitB)
        ids = find(hitB);
        vOld = v(ids, :);
        x(ids, 2) = -x(ids, 2);
        switch modeBottom
            case 'thermalize'
                v(ids, 1) = Ubottom + wallSigma * randn(numel(ids), 1);
                v(ids, 2) = abs(wallSigma * randn(numel(ids), 1));
            case 'bounceback'
                v(ids, 1) = 2*Ubottom - v(ids, 1);
                v(ids, 2) = abs(v(ids, 2));
            case 'specular'
                v(ids, 2) = abs(v(ids, 2));
            otherwise
                error('Unknown wallModeY for bottom wall: %s', modeBottom);
        end
        dv = v(ids, :) - vOld;
        info.nBot = info.nBot + numel(ids);
        info.dEwall = info.dEwall + 0.5 * sum(sum(v(ids, :).^2 - vOld.^2, 2));
        info.dPxBot = info.dPxBot + sum(abs(dv(:, 1)));
        info.dPyBot = info.dPyBot + sum(abs(dv(:, 2)));
        info.dPxBotSigned = info.dPxBotSigned + sum(dv(:, 1));
        info.dPyBotSigned = info.dPyBotSigned + sum(dv(:, 2));
        info.impulseOnBottomWallY = info.impulseOnBottomWallY - sum(dv(:, 2));
    end

    hitT = x(:, 2) > yTop;
    if any(hitT)
        ids = find(hitT);
        vOld = v(ids, :);
        x(ids, 2) = 2*yTop - x(ids, 2);
        switch modeTopTangential
            case 'thermalize'
                v(ids, 1) = Utop + wallSigma * randn(numel(ids), 1);
            case 'bounceback'
                v(ids, 1) = 2*Utop - v(ids, 1);
            case 'specular'
                % tangential velocity unchanged
            otherwise
                error('Unknown pistonTangentialMode: %s', modeTopTangential);
        end
        % Moving-wall normal reflection in the piston frame.
        v(ids, 2) = 2*Up - v(ids, 2);
        dv = v(ids, :) - vOld;
        info.nTop = info.nTop + numel(ids);
        info.dEwall = info.dEwall + 0.5 * sum(sum(v(ids, :).^2 - vOld.^2, 2));
        info.dPxTop = info.dPxTop + sum(abs(dv(:, 1)));
        info.dPyTop = info.dPyTop + sum(abs(dv(:, 2)));
        info.dPxTopSigned = info.dPxTopSigned + sum(dv(:, 1));
        info.dPyTopSigned = info.dPyTopSigned + sum(dv(:, 2));
        info.impulseOnTopWallY = info.impulseOnTopWallY - sum(dv(:, 2));
    end

    if ~any(hitB) && ~any(hitT)
        break;
    end
end

x(:, 2) = min(max(x(:, 2), 0), max(yTop - eps(max(yTop, 1)), 0));
LxLocal = get_param(params, 'Lx', 1.0);
dtLocal = get_param(params, 'dt', 1.0);
info.pressureTopWall = info.impulseOnTopWallY / max(dtLocal * LxLocal, eps);
info.pressureBottomWall = info.impulseOnBottomWallY / max(dtLocal * LxLocal, eps);
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


function value = getf(s, name, defaultValue)
if nargin < 3
    defaultValue = NaN;
end
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
