function [xOut, vOut, closureInfo, params] = mpcd_apply_liquid_closure_step(params, xRef, vRef, xRed, vRed)
%MPCD_APPLY_LIQUID_CLOSURE_STEP
% Fermeture liquide appliquee apres redistribution.
%
% La reference est l'etat physique pre-redistribution issu du step MPCD.
% La redistribution impose la densite cible. La fermeture liquide procede
% ensuite en deux etapes :
%   1) reparation du champ de vitesse detruit par la redistribution ;
%   2) kick viriel base sur Pdrive = Pkin + Pvir avec
%      Pvir = Kvirial * (rho_mpcd - rho_target).
%
% La pression motrice a controler est donc Pdrive construite sur la
% reference pre-redistribution, tandis que l'etat sur lequel on applique
% correction et kick est l'etat redistribue.

statsRef = closure_cell_stats(xRef, vRef, params);
statsRed = closure_cell_stats(xRed, vRed, params);

[targetOcc, rhoTarget, activeFrac] = closure_target_density_map(params);

[duRepairXRaw, duRepairYRaw, repairDefect] = closure_repair_field(statsRef, statsRed, targetOcc);
betaRepairOpt = optimal_scalar_beta(statsRef.ux, statsRef.uy, statsRed.ux, statsRed.uy, duRepairXRaw, duRepairYRaw, activeFrac > 0);
if params.useOptimalBetaRepair
    betaRepair = betaRepairOpt;
else
    betaRepair = params.betaRepair;
end

vRepair = apply_cell_kick_particles(vRed, statsRed.cellId, betaRepair .* duRepairXRaw, betaRepair .* duRepairYRaw);
xRepair = xRed;
statsRepair = closure_cell_stats(xRepair, vRepair, params);

[PkinDrive, PvirDrive, Pdrive] = closure_drive_pressure(statsRef, rhoTarget, params);
[dPdx, dPdy] = gradient_boundary_aware(Pdrive, params);

rhoKick = rhoTarget;
maskActive = rhoKick > 0;
duEOSX = zeros(size(Pdrive));
duEOSY = zeros(size(Pdrive));
duEOSX(maskActive) = -(params.dt ./ rhoKick(maskActive)) .* dPdx(maskActive);
duEOSY(maskActive) = -(params.dt ./ rhoKick(maskActive)) .* dPdy(maskActive);

vEOS = apply_cell_kick_particles(vRepair, statsRepair.cellId, params.betaEOS .* duEOSX, params.betaEOS .* duEOSY);
xEOS = xRepair;

xOut = xEOS;
vOut = vEOS;

closureInfo = struct();
closureInfo.betaRepairOpt = betaRepairOpt;
closureInfo.betaRepairApplied = betaRepair;
closureInfo.betaEOSApplied = params.betaEOS;
closureInfo.targetOcc = targetOcc;
closureInfo.rhoTarget = rhoTarget;
closureInfo.repairDefect = repairDefect;
closureInfo.PkinDrive = PkinDrive;
closureInfo.PvirDrive = PvirDrive;
closureInfo.Pdrive = Pdrive;
closureInfo.duEOSX = duEOSX;
closureInfo.duEOSY = duEOSY;

if params.liquidClosureDumpStages
    closureInfo.dumps.reference = struct('x', xRef, 'v', vRef);
    closureInfo.dumps.redistributed = struct('x', xRed, 'v', vRed);
    closureInfo.dumps.repaired = struct('x', xRepair, 'v', vRepair);
    closureInfo.dumps.postKick = struct('x', xEOS, 'v', vEOS);
    params.lastLiquidClosureDump = closureInfo.dumps;
end
end

function stats = closure_cell_stats(x, v, params)
Nx = params.Nx;
Ny = params.Ny;
aX = params.Lx / params.Nx;
aY = params.Ly / params.Ny;
Vc = aX * aY;

cellId = closure_cell_id_from_pos(x, params);

Nvec = accumarray(cellId, 1, [Nx*Ny, 1], @sum, 0);
Pxvec = accumarray(cellId, v(:,1), [Nx*Ny, 1], @sum, 0);
Pyvec = accumarray(cellId, v(:,2), [Nx*Ny, 1], @sum, 0);

uxvec = zeros(Nx*Ny, 1);
uyvec = zeros(Nx*Ny, 1);
maskN = Nvec > 0;
uxvec(maskN) = Pxvec(maskN) ./ Nvec(maskN);
uyvec(maskN) = Pyvec(maskN) ./ Nvec(maskN);

uxPart = uxvec(cellId);
uyPart = uyvec(cellId);
rel2 = (v(:,1) - uxPart).^2 + (v(:,2) - uyPart).^2;
sumRel2 = accumarray(cellId, rel2, [Nx*Ny, 1], @sum, 0);

dof = 2 * max(Nvec - 1, 0);
kBTloc = nan(Nx*Ny, 1);
maskT = dof > 0;
kBTloc(maskT) = sumRel2(maskT) ./ dof(maskT);

rhoVec = Nvec / Vc;
pkinVec = nan(Nx*Ny, 1);
pkinVec(maskT) = rhoVec(maskT) .* kBTloc(maskT);

stats.cellId = cellId;
stats.Nvec = Nvec;
stats.Pxvec = Pxvec;
stats.Pyvec = Pyvec;
stats.uxvec = uxvec;
stats.uyvec = uyvec;
stats.kBTvec = kBTloc;
stats.rhovec = rhoVec;
stats.pkinvec = pkinVec;

stats.N = reshape(Nvec, [Ny, Nx]);
stats.ux = reshape(uxvec, [Ny, Nx]);
stats.uy = reshape(uyvec, [Ny, Nx]);
stats.kBT = reshape(kBTloc, [Ny, Nx]);
stats.rho = reshape(rhoVec, [Ny, Nx]);
stats.pkin = reshape(pkinVec, [Ny, Nx]);
end

function cellId = closure_cell_id_from_pos(x, params)
Nx = params.Nx;
Ny = params.Ny;
aX = params.Lx / params.Nx;
aY = params.Ly / params.Ny;

ix = floor(x(:,1) / aX) + 1;
iy = floor(x(:,2) / aY) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
cellId = sub2ind([Ny, Nx], iy, ix);
end

function [targetOcc, rhoTarget, activeFrac] = closure_target_density_map(params)
Nx = params.Nx;
Ny = params.Ny;
aX = params.Lx / params.Nx;
aY = params.Ly / params.Ny;
Vc = aX * aY;

gammaEff = closure_effective_gamma(params);
activeFrac = closure_active_cell_fraction_map(params, [Ny, Nx]);

targetOcc = gammaEff .* activeFrac;
rhoTarget = targetOcc ./ Vc;
end

function gammaEff = closure_effective_gamma(params)
if isfield(params, 'gammaCurrent') && ~isempty(params.gammaCurrent)
    gammaEff = params.gammaCurrent;
else
    gammaEff = params.gamma;
end
end

function fracMap = closure_active_cell_fraction_map(params, sz)
Ny = sz(1);
Nx = sz(2);

if isfield(params, 'activeCellFracCurrent') && ~isempty(params.activeCellFracCurrent)
    fracMap = params.activeCellFracCurrent;
    return;
end

fracMap = ones(Ny, Nx);
if isfield(params, 'useMovingPiston') && params.useMovingPiston
    Ly = params.Ly;
    yTop = params.pistonYCurrent;
    aY = Ly / Ny;
    fracRows = zeros(Ny, 1);
    for iy = 1:Ny
        y0 = (iy - 1) * aY;
        y1 = iy * aY;
        fracRows(iy) = max(0, min(yTop, y1) - y0) / aY;
    end
    fracMap = repmat(fracRows, 1, Nx);
end
end

function [duX, duY, defect] = closure_repair_field(statsRef, statsRed, targetOcc)
Ny = size(statsRef.ux, 1);
Nx = size(statsRef.ux, 2);

targetOccSafe = max(targetOcc, 1);
dPx = reshape(statsRef.Pxvec - statsRed.Pxvec, [Ny, Nx]);
dPy = reshape(statsRef.Pyvec - statsRed.Pyvec, [Ny, Nx]);

duX = dPx ./ targetOccSafe;
duY = dPy ./ targetOccSafe;

defect = struct();
defect.dPx = dPx;
defect.dPy = dPy;
defect.duX = duX;
defect.duY = duY;
end

function beta = optimal_scalar_beta(uxRef, uyRef, uxCur, uyCur, duX, duY, mask)
ex = uxCur(mask) - uxRef(mask);
ey = uyCur(mask) - uyRef(mask);
ax = duX(mask);
ay = duY(mask);

den = sum(ax.^2 + ay.^2);
if den <= 0
    beta = 0.0;
    return;
end
num = -sum(ex .* ax + ey .* ay);
beta = max(num / den, 0.0);
end

function [PkinDrive, PvirDrive, Pdrive] = closure_drive_pressure(statsRef, rhoTarget, params)
PkinDrive = statsRef.pkin;
PvirDrive = params.Kvirial .* (statsRef.rho - rhoTarget);
Pdrive = PkinDrive + PvirDrive;
Pdrive(~isfinite(Pdrive)) = 0;
end

function vOut = apply_cell_kick_particles(vIn, cellId, dUx, dUy)
dUxVec = dUx(:);
dUyVec = dUy(:);
vOut = vIn;
vOut(:,1) = vOut(:,1) + dUxVec(cellId);
vOut(:,2) = vOut(:,2) + dUyVec(cellId);
end

function [dFdx, dFdy] = gradient_boundary_aware(F, params)
dx = params.Lx / params.Nx;
dy = params.Ly / params.Ny;
periodicX = closure_is_periodic_pair(params, 'x');
periodicY = closure_is_periodic_pair(params, 'y');

dFdx = ddx_2d(F, dx, periodicX);
dFdy = ddy_2d(F, dy, periodicY);
end

function dFdx = ddx_2d(F, dx, periodicX)
[Ny, Nx] = size(F);
dFdx = zeros(Ny, Nx);
if periodicX
    dFdx = (circshift(F, [0, -1]) - circshift(F, [0, 1])) / (2 * dx);
    return;
end
if Nx == 1
    return;
end
dFdx(:,1) = (F(:,2) - F(:,1)) / dx;
dFdx(:,Nx) = (F(:,Nx) - F(:,Nx-1)) / dx;
if Nx > 2
    dFdx(:,2:Nx-1) = (F(:,3:Nx) - F(:,1:Nx-2)) / (2 * dx);
end
end

function dFdy = ddy_2d(F, dy, periodicY)
[Ny, Nx] = size(F);
dFdy = zeros(Ny, Nx);
if periodicY
    dFdy = (circshift(F, [-1, 0]) - circshift(F, [1, 0])) / (2 * dy);
    return;
end
if Ny == 1
    return;
end
dFdy(1,:) = (F(2,:) - F(1,:)) / dy;
dFdy(Ny,:) = (F(Ny,:) - F(Ny-1,:)) / dy;
if Ny > 2
    dFdy(2:Ny-1,:) = (F(3:Ny,:) - F(1:Ny-2,:)) / (2 * dy);
end
end

function tf = closure_is_periodic_pair(params, axisName)
switch lower(axisName)
    case 'x'
        tf = strcmpi(closure_boundary_mode(params, 'left'), 'periodic') && ...
             strcmpi(closure_boundary_mode(params, 'right'), 'periodic');
    case 'y'
        tf = strcmpi(closure_boundary_mode(params, 'bottom'), 'periodic') && ...
             strcmpi(closure_boundary_mode(params, 'top'), 'periodic');
    otherwise
        error('Axe inconnu: %s', axisName);
end
end

function mode = closure_boundary_mode(params, side)
fieldName = ['boundary_' lower(side)];
if isfield(params, fieldName) && ~isempty(params.(fieldName))
    mode = lower(char(params.(fieldName)));
    return;
end
switch lower(side)
    case {'left','right'}
        mode = lower(char(params.xBoundary));
    case 'bottom'
        mode = lower(char(params.yWallMode));
    case 'top'
        if isfield(params, 'useMovingPiston') && params.useMovingPiston
            mode = 'piston';
        else
            mode = lower(char(params.yWallMode));
        end
    otherwise
        error('Cote inconnu: %s', side);
end
end
