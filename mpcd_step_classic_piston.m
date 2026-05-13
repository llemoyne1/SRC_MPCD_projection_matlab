function [stateOut, diag] = mpcd_step_classic_piston(state, params)
%MPCD_STEP_CLASSIC_PISTON Classic MPCD step in an x-periodic moving-piston channel.
%
%   [stateOut, diag] = mpcd_step_classic_piston(state, params)
%
% This is the Q9-native piston kernel.  It reuses the historical piston
% kinematics in minimal form:
%   - top wall position yTop = pistonYCurrent,
%   - optional affine compression from pistonYPrev to pistonYCurrent,
%   - top wall normal velocity pistonVyCurrent,
%   - x periodicity and a bounded lower wall.
%
% The computational grid is the active domain [0,Lx] x [0,yTop].  The number
% of cells Nx-by-Ny is kept fixed, so dy changes with the piston height.  This
% avoids projecting through an empty inactive layer above the piston.

validate_state(state);

Lx = params.Lx;
Ly0 = get_param(params, 'Ly0', get_param(params, 'LyReference', params.Ly));
Nx = params.Nx;
Ny = params.Ny;
dt = params.dt;
alphaDeg = get_param(params, 'alphaDeg', 90);
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);
useRandomGridShiftX = logical(get_param(params, 'useRandomGridShiftX', true));
useRandomGridShiftY = logical(get_param(params, 'useRandomGridShiftY', false));

if Lx <= 0 || Ly0 <= 0 || Nx <= 0 || Ny <= 0 || dt <= 0
    error('Lx, Ly, Nx, Ny and dt must be strictly positive.');
end

yPrev = get_param(params, 'pistonYPrev', Ly0);
yTop = get_param(params, 'pistonYCurrent', Ly0);
Up = get_param(params, 'pistonVyCurrent', get_param(params, 'pistonVy', 0.0));
yTop = min(max(yTop, eps(Ly0)), Ly0);
yPrev = min(max(yPrev, eps(Ly0)), Ly0);

x = state.x;
v = state.v;
Np = size(x, 1);

% Quasi-static geometric compression inherited from the historical piston
% implementation.  It spreads the imposed volume change through the fluid
% instead of letting only the upper layer absorb the piston displacement.
if logical(get_param(params, 'pistonAffineReposition', true)) && abs(yTop - yPrev) > eps(Ly0)
    scale = yTop / max(yPrev, eps(Ly0));
    x(:, 2) = min(max(x(:, 2), 0), max(yPrev - eps(yPrev), 0));
    x(:, 2) = min(max(x(:, 2) * scale, 0), max(yTop - eps(yTop), 0));
end

dx = Lx / Nx;
dy = yTop / Ny;
Nc = Nx * Ny;

% Force kick and streaming.
v(:, 1) = v(:, 1) + dt * bodyForceX;
v(:, 2) = v(:, 2) + dt * bodyForceY;
x(:, 1) = mod(x(:, 1) + dt * v(:, 1), Lx);
x(:, 2) = x(:, 2) + dt * v(:, 2);
[x, v, wallInfo] = apply_piston_wall_bc_y(x, v, params, yTop, Up);

% Collision-cell assignment on the active piston height.
if useRandomGridShiftX
    shiftX = (rand() - 0.5) * dx;
else
    shiftX = 0.0;
end
if useRandomGridShiftY
    shiftY = (rand() - 0.5) * dy;
else
    shiftY = 0.0;
end

xs = mod(x(:, 1) + shiftX, Lx);
ys = min(max(x(:, 2) + shiftY, 0), yTop - eps(yTop));
ix = floor(xs / dx) + 1;
iy = floor(ys / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
cellId = iy + Ny * (ix - 1);

Ncell = accumarray(cellId, 1, [Nc, 1], @sum, 0);
Px = accumarray(cellId, v(:, 1), [Nc, 1], @sum, 0);
Py = accumarray(cellId, v(:, 2), [Nc, 1], @sum, 0);

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
rvx = v(:, 1) - uxp;
rvy = v(:, 2) - uyp;

v(:, 1) = uxp + ca .* rvx - sa .* rvy;
v(:, 2) = uyp + sa .* rvx + ca .* rvy;

stateOut = state;
stateOut.x = x;
stateOut.v = v;
stateOut.piston = struct('yTop', yTop, 'yPrev', yPrev, 'Up', Up, ...
    'Ly0', Ly0, 'compression', 1.0 - yTop / Ly0, 'activeHeight', yTop);

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
diag.wallInfo = wallInfo;
diag.piston = stateOut.piston;
diag.NMean = mean(Ncell);
diag.NStd = std(double(Ncell));
diag.NMin = min(Ncell);
diag.NMax = max(Ncell);
diag.nEmptyCells = nnz(Ncell == 0);
diag.meanVx = mean(v(:, 1));
diag.meanVy = mean(v(:, 2));
diag.kBT = estimate_kBT(v);
diag.kineticEnergyMean = 0.5 * mean(sum(v.^2, 2));
end

function [x, v, info] = apply_piston_wall_bc_y(x, v, params, yTop, Up)
modeBottom = lower(char(string(get_param(params, 'wallModeY', 'bounceback'))));
modeTopTangential = lower(char(string(get_param(params, 'pistonTangentialMode', 'specular'))));
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

for pass = 1:8
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
        % Moving piston normal reflection.
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

x(:, 2) = min(max(x(:, 2), 0), max(yTop - eps(yTop), 0));

% Instantaneous wall pressures averaged over this MPCD step.
% pressureTopWall is positive when the fluid exerts an upward force on the
% piston. The impulse convention is the impulse received by the wall, i.e.
% minus the particle momentum change. In 2D the wall area is the box length Lx.
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

function kBT = estimate_kBT(v)
u = mean(v, 1);
c = v - u;
kBT = 0.5 * mean(sum(c.^2, 2));
end
