function [x, v, type, r0, stepInfo, params] = mpcd_incompressible_step_liquidclosure(params, x, v, type, r0, it)
%MPCD_INCOMPRESSIBLE_STEP_LIQUIDCLOSURE
% Coeur MPCD sur un seul pas de temps avec fermeture liquide externe.
%
% Sequence retenue :
%   MPCD natif -> reference pre-redistribution -> redistribution ->
%   correction + kick viriel -> reorientation interfaciale.
% La reorientation ne fait pas partie de la reference.

bodyForceX = params.bodyForceX;
g = params.g;
dt = params.dt;
keepMeanFlow = params.keepMeanFlow;
useIncompressibleRedistribution = params.useIncompressibleRedistribution;
redistribAfterCollision = params.redistribAfterCollision;
useLiquidClosure = params.useLiquidClosure;

auditDumpStages = params.liquidClosureDumpStages;

params = update_runtime_moving_piston(params, it);
[x, r0] = apply_affine_piston_compression(x, r0, params);
gammaCell = get_effective_gamma(params);

% 1) Streaming + forcage volumique
v(:,1) = v(:,1) + bodyForceX*dt;
v(:,2) = v(:,2) + g*dt;
xTrial = x + dt*v;

% 2) Conditions limites
[x, v, wallInfo] = apply_bc_general(xTrial, v, params);
pBotInst = wallInfo.dPyBot / max(params.Lx*dt, eps);
pTopInst = wallInfo.dPyTop / max(params.Lx*dt, eps);
pMeanInst = 0.5*(pBotInst + pTopInst);

% 3) Collision SRD : l'etat obtenu ici devient la reference physique
% pre-redistribution pour la fermeture liquide.
cid = srd_cell_id_with_random_shift_general(x, params);
v = srd_collision_step_fast(v, cid, params);
if keepMeanFlow
    v(:,1) = v(:,1) - mean(v(:,1));
    v(:,2) = v(:,2) - mean(v(:,2));
end
xRef = x;
vRef = v;

% 4) Redistribution incompressible
if useIncompressibleRedistribution && redistribAfterCollision
    [x, v, stats, Ncell, diag] = enforce_cell_occupancy_local_gradient_fast(x, v, false, params);
else
    Ncell = cell_occupancy(x, params);
    diag = default_redist_diag(Ncell, params);
    stats = zeros(1,6);
end
xRed = x;
vRed = v;

% 5) Reorientation interfaciale 
if isfield(params,'useInterfaceVelocityReorientation') && params.useInterfaceVelocityReorientation
    if auditDumpStages
        xBeforeReorient = x;
        vBeforeReorient = v;
    end
    v = apply_interface_velocity_reorientation_step(v, x, params);
    if auditDumpStages
        closureInfo.dumps.preReorientation = struct('x', xBeforeReorient, 'v', vBeforeReorient);
        closureInfo.dumps.postReorientation = struct('x', x, 'v', v);
    end
end


% 6) Fermeture liquide externe : reparation puis kick viriel construit
% depuis la reference pre-redistribution.
if useLiquidClosure
    [x, v, closureInfo, params] = mpcd_apply_liquid_closure_step(params, xRef, vRef, xRed, vRed);
else
    closureInfo = struct();
    if auditDumpStages
        closureInfo.dumps.reference = struct('x', xRef, 'v', vRef);
        closureInfo.dumps.redistributed = struct('x', xRed, 'v', vRed);
        closureInfo.dumps.repaired = struct('x', xRed, 'v', vRed);
        closureInfo.dumps.postKick = struct('x', xRed, 'v', vRed);
    end
end


stepInfo = struct();
stepInfo.wallInfo = wallInfo;
stepInfo.Ncell = Ncell;
stepInfo.redistDiag = diag;
stepInfo.redistribStats = stats;
stepInfo.gammaCell = gammaCell;
stepInfo.pBotInst = pBotInst;
stepInfo.pTopInst = pTopInst;
stepInfo.pMeanInst = pMeanInst;
stepInfo.pRow = [wallInfo.dPyBot, wallInfo.dPyTop, pBotInst, pTopInst, pMeanInst];
stepInfo.liquidClosure = closureInfo;
stepInfo.referenceCaptured = true;

if auditDumpStages
    params.lastLiquidClosureDump = closureInfo.dumps;
end

end

function params = prepare_step_params(params, x)
if nargin < 2, x = []; end
if ~isfield(params,'n') || isempty(params.n)
    if ~isempty(x)
        params.n = size(x,1);
    else
        params.n = 0;
    end
else
    if ~isempty(x), params.n = size(x,1); end
end
if ~isfield(params, 'Lx') || ~isfield(params, 'Ly') || ~isfield(params, 'Nx') || ~isfield(params, 'Ny')
    error('params doit contenir au minimum Lx, Ly, Nx et Ny.');
end
if ~isfield(params, 'a0') || isempty(params.a0), params.a0 = params.Lx / params.Nx; end
if ~isfield(params, 'Nc') || isempty(params.Nc), params.Nc = params.Nx * params.Ny; end
if ~isfield(params, 'dt') || isempty(params.dt), params.dt = 5e-3; end
if ~isfield(params, 'alphaDeg') || isempty(params.alphaDeg), params.alphaDeg = 170; end
if ~isfield(params, 'alpha') || isempty(params.alpha), params.alpha = deg2rad(params.alphaDeg); end
if ~isfield(params, 'kBT') || isempty(params.kBT), params.kBT = 1.0; end
if ~isfield(params, 'g') || isempty(params.g), params.g = 0.0; end
if ~isfield(params, 'bodyForceX') || isempty(params.bodyForceX), params.bodyForceX = 0.0; end
if ~isfield(params, 'useThermostat') || isempty(params.useThermostat), params.useThermostat = true; end
if ~isfield(params, 'keepMeanFlow') || isempty(params.keepMeanFlow), params.keepMeanFlow = false; end
if ~isfield(params, 'xBoundary') || isempty(params.xBoundary), params.xBoundary = 'specular'; end
if ~isfield(params, 'yWallMode') || isempty(params.yWallMode), params.yWallMode = 'specular'; end
if ~isfield(params, 'boundary_left') || isempty(params.boundary_left), params.boundary_left = params.xBoundary; end
if ~isfield(params, 'boundary_right') || isempty(params.boundary_right), params.boundary_right = params.xBoundary; end
if ~isfield(params, 'boundary_bottom') || isempty(params.boundary_bottom), params.boundary_bottom = params.yWallMode; end
if ~isfield(params, 'boundary_top') || isempty(params.boundary_top)
    if isfield(params,'useMovingPiston') && params.useMovingPiston
        params.boundary_top = 'piston';
    else
        params.boundary_top = params.yWallMode;
    end
end
if ~isfield(params, 'Utop') || isempty(params.Utop), params.Utop = 0; end
if ~isfield(params, 'Ubottom') || isempty(params.Ubottom), params.Ubottom = 0; end
if ~isfield(params, 'wallSigma') || isempty(params.wallSigma), params.wallSigma = sqrt(max(params.kBT,0)); end
if ~isfield(params, 'useMovingPiston') || isempty(params.useMovingPiston), params.useMovingPiston = false; end
if isfield(params,'boundary_top') && ~isempty(params.boundary_top) && strcmpi(params.boundary_top,'piston')
    params.useMovingPiston = true;
end
if ~isfield(params, 'pistonY0') || isempty(params.pistonY0), params.pistonY0 = params.Ly; end
if ~isfield(params, 'pistonVy') || isempty(params.pistonVy), params.pistonVy = 0; end
if ~isfield(params, 'pistonYMin') || isempty(params.pistonYMin), params.pistonYMin = 0.5*params.Ly; end
if ~isfield(params, 'pistonStopOnMin') || isempty(params.pistonStopOnMin), params.pistonStopOnMin = true; end
if ~isfield(params, 'occAbsFloor') || isempty(params.occAbsFloor), params.occAbsFloor = 1; end
if ~isfield(params, 'pistonActiveMargin') || isempty(params.pistonActiveMargin), params.pistonActiveMargin = 0; end
if ~isfield(params, 'pistonFracCut') || isempty(params.pistonFracCut), params.pistonFracCut = 0.20; end
if ~isfield(params, 'pistonAffineReposition') || isempty(params.pistonAffineReposition), params.pistonAffineReposition = true; end
if ~isfield(params, 'gamma') || isempty(params.gamma)
    params.gamma = params.n / max(params.Nx*params.Ny,1);
end
if ~isfield(params, 'useIncompressibleRedistribution') || isempty(params.useIncompressibleRedistribution), params.useIncompressibleRedistribution = true; end
if ~isfield(params, 'redistribAfterCollision') || isempty(params.redistribAfterCollision), params.redistribAfterCollision = true; end
if ~isfield(params, 'coef') || isempty(params.coef), params.coef = 0.20; end
if ~isfield(params, 'highMode') || isempty(params.highMode), params.highMode = 'coef'; end
if ~isfield(params, 'lowMode') || isempty(params.lowMode), params.lowMode = 'coef'; end
if ~isfield(params, 'maxRedistribPasses') || isempty(params.maxRedistribPasses), params.maxRedistribPasses = 2; end
if ~isfield(params, 'enableMomentumCorrectionPostRedistribution') || isempty(params.enableMomentumCorrectionPostRedistribution)
    params.enableMomentumCorrectionPostRedistribution = true;
end
if ~isfield(params, 'useLocalFluidFractionThresholds') || isempty(params.useLocalFluidFractionThresholds), params.useLocalFluidFractionThresholds = true; end
if ~isfield(params, 'fluidFracGain') || isempty(params.fluidFracGain), params.fluidFracGain = 0.8; end
if ~isfield(params, 'lowThrFloor') || isempty(params.lowThrFloor), params.lowThrFloor = 20; end
if ~isfield(params, 'highThrFloor') || isempty(params.highThrFloor), params.highThrFloor = 25; end
if ~isfield(params, 'lowThrBulkOverride') || isempty(params.lowThrBulkOverride), params.lowThrBulkOverride = NaN; end
if ~isfield(params, 'lowThrInterfaceOverride') || isempty(params.lowThrInterfaceOverride), params.lowThrInterfaceOverride = 10; end
if ~isfield(params, 'useBottomWallWetting') || isempty(params.useBottomWallWetting), params.useBottomWallWetting = false; end
if ~isfield(params, 'bottomWettingTargetOverride') || isempty(params.bottomWettingTargetOverride), params.bottomWettingTargetOverride = NaN; end
if ~isfield(params, 'bottomWettingMaxRange') || isempty(params.bottomWettingMaxRange), params.bottomWettingMaxRange = inf; end
if ~isfield(params, 'useInterfaceVelocityReorientation') || isempty(params.useInterfaceVelocityReorientation), params.useInterfaceVelocityReorientation = false; end
if ~isfield(params, 'interfaceReorientBeta') || isempty(params.interfaceReorientBeta), params.interfaceReorientBeta = 0.36; end
if ~isfield(params, 'interfaceReorientFrac') || isempty(params.interfaceReorientFrac), params.interfaceReorientFrac = 0.35; end
if ~isfield(params, 'useCurvatureWeightedReorientation') || isempty(params.useCurvatureWeightedReorientation), params.useCurvatureWeightedReorientation = true; end
if ~isfield(params, 'curvatureReorientBetaGain') || isempty(params.curvatureReorientBetaGain), params.curvatureReorientBetaGain = 5; end
if ~isfield(params, 'curvatureReorientFracGain') || isempty(params.curvatureReorientFracGain), params.curvatureReorientFracGain = 5; end
if ~isfield(params, 'interfaceReorientBetaMax') || isempty(params.interfaceReorientBetaMax), params.interfaceReorientBetaMax = 1; end
if ~isfield(params, 'interfaceReorientFracMax') || isempty(params.interfaceReorientFracMax), params.interfaceReorientFracMax = 1; end
if ~isfield(params, 'curvatureAngleMax') || isempty(params.curvatureAngleMax), params.curvatureAngleMax = pi/6; end

legacySurfaceTopology = true;
if isfield(params,'isSurface') && ~isempty(params.isSurface)
    legacySurfaceTopology = logical(params.isSurface);
elseif isfield(params,'caseType') && ~isempty(params.caseType)
    caseTypeLegacy = lower(char(string(params.caseType)));
    legacySurfaceTopology = strcmp(caseTypeLegacy,'ellipse');
end
if ~isfield(params, 'noSurfaceCase') || isempty(params.noSurfaceCase)
    params.noSurfaceCase = ~legacySurfaceTopology;
end
if ~isfield(params, 'useLocalTargetsForPiston') || isempty(params.useLocalTargetsForPiston)
    params.useLocalTargetsForPiston = params.useMovingPiston;
end
if ~isfield(params, 'redistributionEnableSurfaceTopology') || isempty(params.redistributionEnableSurfaceTopology)
    params.redistributionEnableSurfaceTopology = legacySurfaceTopology;
end
if ~isfield(params, 'redistributionUseActiveVolumeTargets') || isempty(params.redistributionUseActiveVolumeTargets)
    params.redistributionUseActiveVolumeTargets = params.useLocalTargetsForPiston;
end
if ~isfield(params, 'redistributionBulkOverflowPolicy') || isempty(params.redistributionBulkOverflowPolicy)
    if params.redistributionEnableSurfaceTopology
        params.redistributionBulkOverflowPolicy = 'nearest_interface';
    else
        params.redistributionBulkOverflowPolicy = 'none';
    end
end
if ~isfield(params, 'redistributionWallWettingEnabled') || isempty(params.redistributionWallWettingEnabled)
    params.redistributionWallWettingEnabled = ...
        strcmpi(get_boundary_mode_general(params,'left'),   'wall_wetting') || ...
        strcmpi(get_boundary_mode_general(params,'right'),  'wall_wetting') || ...
        strcmpi(get_boundary_mode_general(params,'bottom'), 'wall_wetting') || ...
        strcmpi(get_boundary_mode_general(params,'top'),    'wall_wetting');
    if isfield(params,'useBottomWallWetting') && ~isempty(params.useBottomWallWetting)
        params.redistributionWallWettingEnabled = params.redistributionWallWettingEnabled || params.useBottomWallWetting;
    end
end
if ~isfield(params, 'redistributionPreserveWettingLayer') || isempty(params.redistributionPreserveWettingLayer)
    params.redistributionPreserveWettingLayer = params.redistributionWallWettingEnabled;
end
if ~isfield(params, 'redistributionWettingTargetOverride') || isempty(params.redistributionWettingTargetOverride)
    if isfield(params,'bottomWettingTargetOverride') && ~isempty(params.bottomWettingTargetOverride)
        params.redistributionWettingTargetOverride = params.bottomWettingTargetOverride;
    else
        params.redistributionWettingTargetOverride = NaN;
    end
end
if ~isfield(params, 'redistributionWettingMaxRange') || isempty(params.redistributionWettingMaxRange)
    if isfield(params,'bottomWettingMaxRange') && ~isempty(params.bottomWettingMaxRange)
        params.redistributionWettingMaxRange = params.bottomWettingMaxRange;
    else
        params.redistributionWettingMaxRange = inf;
    end
end
if ~isfield(params, 'redistributionWettingSides') || isempty(params.redistributionWettingSides)
    params.redistributionWettingSides = resolve_redistribution_wetting_sides(params);
end
if ~isfield(params, 'neighborsCore') || isempty(params.neighborsCore)
    params.neighborsCore = precompute_neighbors_general(params);
end
end

function [x, v, info] = apply_bc_general(x, v, params)
Lx = params.Lx;
Ly = params.Ly;
leftMode   = get_boundary_mode_general(params, 'left');
rightMode  = get_boundary_mode_general(params, 'right');
bottomMode = get_boundary_mode_general(params, 'bottom');
topMode    = get_boundary_mode_general(params, 'top');
Utop = params.Utop;
Ubottom = params.Ubottom;
wallSigma = params.wallSigma;
useMovingPiston = isfield(params,'useMovingPiston') && params.useMovingPiston;

yTop = resolve_top_boundary_position_general(params, Ly, topMode);
Up   = resolve_top_boundary_velocity_general(params, Utop, topMode);

info = struct('nBot',0,'nTop',0,'nLeft',0,'nRight',0, ...
              'dEwall',0,'dPxBot',0,'dPxTop',0,'dPyBot',0,'dPyTop',0, ...
              'dPxLeft',0,'dPxRight',0,'dPyLeft',0,'dPyRight',0);

for k = 1:4
    % --- Limite gauche ---
    hitL = x(:,1) < 0;
    if any(hitL)
        ids = find(hitL);
        vOld = v(ids,:);
        switch leftMode
            case 'periodic'
                x(ids,1) = x(ids,1) + Lx;
            case {'thermalize','bounceback','specular','wall_wetting','piston'}
                x(ids,1) = -x(ids,1);
                switch leftMode
                    case 'thermalize'
                        v(ids,1) = abs(wallSigma*randn(numel(ids),1));
                        v(ids,2) = wallSigma*randn(numel(ids),1);
                    case 'bounceback'
                        v(ids,1) = abs(v(ids,1));
                        v(ids,2) = -v(ids,2);
                    otherwise % specular / wall_wetting / piston
                        v(ids,1) = abs(v(ids,1));
                end
                dv = v(ids,:) - vOld;
                info.nLeft = info.nLeft + numel(ids);
                info.dEwall = info.dEwall + 0.5*sum(sum(v(ids,:).^2 - vOld.^2,2));
                info.dPxLeft = info.dPxLeft + sum(abs(dv(:,1)));
                info.dPyLeft = info.dPyLeft + sum(dv(:,2));
            otherwise
                error('Condition limite gauche non reconnue: %s', leftMode);
        end
    end

    % --- Limite droite ---
    hitR = x(:,1) > Lx;
    if any(hitR)
        ids = find(hitR);
        vOld = v(ids,:);
        switch rightMode
            case 'periodic'
                x(ids,1) = x(ids,1) - Lx;
            case {'thermalize','bounceback','specular','wall_wetting','piston'}
                x(ids,1) = 2*Lx - x(ids,1);
                switch rightMode
                    case 'thermalize'
                        v(ids,1) = -abs(wallSigma*randn(numel(ids),1));
                        v(ids,2) = wallSigma*randn(numel(ids),1);
                    case 'bounceback'
                        v(ids,1) = -abs(v(ids,1));
                        v(ids,2) = -v(ids,2);
                    otherwise % specular / wall_wetting / piston
                        v(ids,1) = -abs(v(ids,1));
                end
                dv = v(ids,:) - vOld;
                info.nRight = info.nRight + numel(ids);
                info.dEwall = info.dEwall + 0.5*sum(sum(v(ids,:).^2 - vOld.^2,2));
                info.dPxRight = info.dPxRight + sum(abs(dv(:,1)));
                info.dPyRight = info.dPyRight + sum(dv(:,2));
            otherwise
                error('Condition limite droite non reconnue: %s', rightMode);
        end
    end

    % --- Limite basse ---
    hitB = x(:,2) < 0;
    if any(hitB)
        ids = find(hitB);
        vOld = v(ids,:);
        switch bottomMode
            case 'periodic'
                x(ids,2) = x(ids,2) + yTop;
            case 'thermalize'
                x(ids,2) = -x(ids,2);
                v(ids,1) = Ubottom + wallSigma*randn(numel(ids),1);
                v(ids,2) = abs(wallSigma*randn(numel(ids),1));
            case 'bounceback'
                x(ids,2) = -x(ids,2);
                v(ids,1) = 2*Ubottom - v(ids,1);
                v(ids,2) = abs(v(ids,2));
            case {'specular','wall_wetting','piston'}
                x(ids,2) = -x(ids,2);
                v(ids,2) = abs(v(ids,2));
            otherwise
                error('Condition limite basse non reconnue: %s', bottomMode);
        end
        if ~strcmp(bottomMode,'periodic')
            dv = v(ids,:) - vOld;
            info.nBot = info.nBot + numel(ids);
            info.dEwall = info.dEwall + 0.5*sum(sum(v(ids,:).^2 - vOld.^2,2));
            info.dPxBot = info.dPxBot + sum(dv(:,1));
            info.dPyBot = info.dPyBot + sum(abs(dv(:,2)));
        end
    end

    % --- Limite haute ---
    hitT = x(:,2) > yTop;
    if any(hitT)
        ids = find(hitT);
        vOld = v(ids,:);
        switch topMode
            case 'periodic'
                x(ids,2) = x(ids,2) - yTop;
            case 'thermalize'
                x(ids,2) = 2*yTop - x(ids,2);
                v(ids,1) = Utop + wallSigma*randn(numel(ids),1);
                v(ids,2) = -abs(wallSigma*randn(numel(ids),1));
            case 'bounceback'
                x(ids,2) = 2*yTop - x(ids,2);
                v(ids,1) = 2*Utop - v(ids,1);
                v(ids,2) = -abs(v(ids,2));
            case {'specular','wall_wetting'}
                x(ids,2) = 2*yTop - x(ids,2);
                v(ids,2) = -abs(v(ids,2));
            case 'piston'
                x(ids,2) = 2*yTop - x(ids,2);
                v(ids,2) = 2*Up - v(ids,2);
            otherwise
                error('Condition limite haute non reconnue: %s', topMode);
        end
        if ~strcmp(topMode,'periodic')
            dv = v(ids,:) - vOld;
            info.nTop = info.nTop + numel(ids);
            info.dEwall = info.dEwall + 0.5*sum(sum(v(ids,:).^2 - vOld.^2,2));
            info.dPxTop = info.dPxTop + sum(dv(:,1));
            info.dPyTop = info.dPyTop + sum(abs(dv(:,2)));
        end
    end
end

if is_periodic_pair_general(params, 'x')
    x(:,1) = x(:,1) - floor(x(:,1)/Lx)*Lx;
else
    x(:,1) = min(max(x(:,1), 0), Lx);
end

if is_periodic_pair_general(params, 'y')
    x(:,2) = x(:,2) - floor(x(:,2)/max(yTop,eps))*max(yTop,eps);
elseif strcmpi(topMode,'piston') || useMovingPiston
    x(:,2) = min(max(x(:,2), 0), max(yTop - eps, 0));
else
    x(:,2) = min(max(x(:,2), 0), yTop);
end
end

function mode = get_boundary_mode_general(params, side)
fieldName = ['boundary_' lower(side)];
if isfield(params, fieldName) && ~isempty(params.(fieldName))
    mode = lower(char(params.(fieldName)));
    return;
end

switch lower(side)
    case {'left','right'}
        if isfield(params,'xBoundary') && ~isempty(params.xBoundary)
            mode = lower(params.xBoundary);
        else
            mode = 'specular';
        end
    case 'bottom'
        if isfield(params,'yWallMode') && ~isempty(params.yWallMode)
            mode = lower(params.yWallMode);
        else
            mode = 'specular';
        end
    case 'top'
        if isfield(params,'useMovingPiston') && params.useMovingPiston
            mode = 'piston';
        elseif isfield(params,'yWallMode') && ~isempty(params.yWallMode)
            mode = lower(params.yWallMode);
        else
            mode = 'specular';
        end
    otherwise
        error('Cote inconnu: %s', side);
end
end

function tf = is_periodic_pair_general(params, axisName)
switch lower(axisName)
    case 'x'
        tf = strcmpi(get_boundary_mode_general(params,'left'), 'periodic') && ...
             strcmpi(get_boundary_mode_general(params,'right'), 'periodic');
    case 'y'
        tf = strcmpi(get_boundary_mode_general(params,'bottom'), 'periodic') && ...
             strcmpi(get_boundary_mode_general(params,'top'), 'periodic');
    otherwise
        error('Axe inconnu: %s', axisName);
end
end

function yTop = resolve_top_boundary_position_general(params, Ly, topMode)
if nargin < 3 || isempty(topMode)
    topMode = get_boundary_mode_general(params, 'top');
end
if strcmpi(topMode,'piston') || (isfield(params,'useMovingPiston') && params.useMovingPiston)
    if isfield(params,'pistonYCurrent') && ~isempty(params.pistonYCurrent)
        yTop = params.pistonYCurrent;
    elseif isfield(params,'pistonY0') && ~isempty(params.pistonY0)
        yTop = params.pistonY0;
    else
        yTop = Ly;
    end
else
    yTop = Ly;
end
yTop = min(max(yTop, 0), Ly);
end

function Up = resolve_top_boundary_velocity_general(params, Utop, topMode)
if nargin < 3 || isempty(topMode)
    topMode = get_boundary_mode_general(params, 'top');
end
if strcmpi(topMode,'piston') || (isfield(params,'useMovingPiston') && params.useMovingPiston)
    if isfield(params,'pistonVyCurrent') && ~isempty(params.pistonVyCurrent)
        Up = params.pistonVyCurrent;
    elseif isfield(params,'pistonVy') && ~isempty(params.pistonVy)
        Up = params.pistonVy;
    else
        Up = 0;
    end
else
    Up = Utop;
end
end

function cid = srd_cell_id_with_random_shift_general(x, params)
Lx = params.Lx;
Ly = params.Ly;
a0 = params.a0;
Nx = params.Nx;
Ny = params.Ny;
shift = (rand(1,2) - 0.5)*a0;
xs = x + shift;
if is_periodic_pair_general(params,'x')
    xs(:,1) = xs(:,1) - floor(xs(:,1)/Lx)*Lx;
else
    xs(:,1) = min(max(xs(:,1), 0), Lx - eps);
end
topMode = get_boundary_mode_general(params, 'top');
yTop = resolve_top_boundary_position_general(params, Ly, topMode);
if is_periodic_pair_general(params,'y')
    xs(:,2) = xs(:,2) - floor(xs(:,2)/max(yTop,eps))*max(yTop,eps);
else
    xs(:,2) = min(max(xs(:,2), 0), yTop - eps);
end
ix = floor(xs(:,1)/a0) + 1;
iy = floor(xs(:,2)/a0) + 1;
ix = max(1, min(Nx, ix));
iy = max(1, min(Ny, iy));
cid = sub2ind([Ny, Nx], iy, ix);
end

function v = srd_collision_step_fast(v, cid, params)
Nc = params.Nc;
alpha = params.alpha;
kBT = params.kBT;
useThermostat = params.useThermostat;
cnt  = accumarray(cid, 1,      [Nc,1], @sum, 0);
sumx = accumarray(cid, v(:,1), [Nc,1], @sum, 0);
sumy = accumarray(cid, v(:,2), [Nc,1], @sum, 0);
ux = sumx ./ max(cnt,1);
uy = sumy ./ max(cnt,1);

vxr = v(:,1) - ux(cid);
vyr = v(:,2) - uy(cid);

sgn = ones(Nc,1);
sgn(rand(Nc,1) < 0.5) = -1;
theta = alpha * sgn;
theta(cnt < 2) = 0;

c = cos(theta); s = sin(theta);
ci = c(cid); si = s(cid);

vxr2 =  ci.*vxr - si.*vyr;
vyr2 =  si.*vxr + ci.*vyr;

if useThermostat
    rel2 = vxr2.^2 + vyr2.^2;
    sumRel2 = accumarray(cid, rel2, [Nc,1], @sum, 0);
    dof = 2*max(cnt - 1, 0);
    target = dof * kBT;
    lambda = ones(Nc,1);
    ok = (cnt > 1) & (sumRel2 > 0);
    lambda(ok) = sqrt(target(ok) ./ sumRel2(ok));
    lam = lambda(cid);
    vxr2 = lam .* vxr2;
    vyr2 = lam .* vyr2;
end

v(:,1) = ux(cid) + vxr2;
v(:,2) = uy(cid) + vyr2;
end

function [x, v, stats, Ncell, diag] = enforce_cell_occupancy_local_gradient_fast( ...
    x, v, doDivDiag, params)

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
gamma = get_effective_gamma(params);
activeFracMap = get_active_cell_fraction_map(params, [Ny, Nx]);
activeMask = get_active_cell_mask(params, [Ny, Nx]);
activeMaskVec = activeMask(:);
neighborsCore = params.neighborsCore;
flags = resolve_redistribution_flags(params);

stats = zeros(1,6);
diag  = zeros(1,10);

[cntVec, partsInCell] = rebuild_parts_in_cell(x, params);
[PbeforeX, PbeforeY] = cell_total_momentum(partsInCell, v, params);
recvMask = false(Nx*Ny,1);
momCorrBulkBulkMask = false(Nx*Ny,1);
momCorrSkipMask = false(Nx*Ny,1);

PglobBefore = sum(v,1);
Ebefore = 0.5 * sum(v(:,1).^2 + v(:,2).^2);
[highThr, lowThrBulk, ~] = occupancy_thresholds(params);

enableMomentumCorrection = params.enableMomentumCorrectionPostRedistribution;
pass = 0;
keepRedistributing = true; %#ok<NASGU>

while pass < params.maxRedistribPasses
    pass = pass + 1;
    nHigh = 0; nLow = 0; movedOut = 0; movedIn = 0;

    Nmat = reshape(cntVec, [Ny, Nx]);
    Nsm  = smooth_occupancy_field(Nmat);
    [~, ~, lowThrLocMap, highThrLocMap] = ...
        build_local_fluid_fraction_threshold_maps(Nsm, params);

    bulkMask = build_topological_bulk_mask(cntVec, params) & activeMaskVec;
    if ~flags.enableSurfaceTopology
        bulkMask = (cntVec > 0) & activeMaskVec;
    end
    lowThrActiveMap  = lowThrLocMap;
    highThrActiveMap = highThrLocMap;
    if ~flags.useActiveVolumeTargets
        lowThrActiveMap(bulkMask)  = lowThrBulk;
        highThrActiveMap(bulkMask) = highThr;
    end
    bulkMaskDenseBefore = bulkMask;

    % ================================================================
    % 1) Cellules trop denses
    % ================================================================
    highThrVec = reshape(highThrActiveMap, [], 1);
    denseList = find(activeMaskVec & (cntVec > highThrVec));
    if ~isempty(denseList)
        denseList = denseList(randperm(numel(denseList)));
    end

    for ii = 1:numel(denseList)
        c = denseList(ii);
        nc = cntVec(c);
        if nc <= 0, continue; end

        [iyc, ixc] = ind2sub([Ny, Nx], c);
        highThrLocal = highThrActiveMap(iyc, ixc);
        if nc <= highThrLocal
            continue;
        end

        nHigh = nHigh + 1;
        wasBulkDense = (~flags.enableSurfaceTopology) || is_bulk_like_cell(c, cntVec, params);

        if ~flags.enableSurfaceTopology
            nHat = [0, 0];
            inwardWeight = 0.0;
        else
            [nHat, inwardWeight] = dense_cell_inward_bias(c, cntVec, Nsm, params);
        end

        [x, partsInCell, cntVec, movedNow, recvTargets] = redistribute_dense_cell_at_range( ...
            x, partsInCell, cntVec, c, lowThrActiveMap, highThrActiveMap, ...
            highThrLocal, nHat, inwardWeight, params, 2);
        if ~isempty(recvTargets), recvMask(recvTargets) = true; end
        [momCorrBulkBulkMask, momCorrSkipMask] = update_transfer_correction_masks( ...
            momCorrBulkBulkMask, momCorrSkipMask, c, recvTargets, bulkMaskDenseBefore, false);
        movedOut = movedOut + movedNow;
        movedIn  = movedIn  + movedNow;

        if flags.enableSurfaceTopology && wasBulkDense && (cntVec(c) > highThrLocal) && ...
                strcmpi(flags.bulkOverflowPolicy, 'nearest_interface')
            [x, partsInCell, cntVec, movedNow, recvTargets, saturated] = ...
                report_dense_bulk_surplus_to_nearest_interface( ...
                    x, partsInCell, cntVec, c, highThrLocal, params);
            if ~isempty(recvTargets), recvMask(recvTargets) = true; end
            [momCorrBulkBulkMask, momCorrSkipMask] = update_transfer_correction_masks( ...
                momCorrBulkBulkMask, momCorrSkipMask, c, recvTargets, bulkMaskDenseBefore, true);
            movedOut = movedOut + movedNow;
            movedIn  = movedIn  + movedNow;

            if saturated
                fprintf(['ALERTE saturation redistribution bulk->interface : pass=%d, cell=%d, cnt=%g, highThr=%.3f\n'], ...
                        pass, c, cntVec(c), highThrLocal);
            end
        end
    end

    % ================================================================
    % 1b) Etalement tangent des couches de mouillage le long des bords
    % explicitement declarés "wall_wetting"
    % ================================================================
    if flags.wallWettingEnabled
        wettingTarget = resolve_wall_wetting_target(lowThrBulk, params);
        [x, partsInCell, cntVec, movedWet, recvWet, srcWet] = apply_wall_wetting_layers_general( ...
            x, partsInCell, cntVec, wettingTarget, params, flags);
        if ~isempty(recvWet), recvMask(recvWet) = true; end
        [momCorrBulkBulkMask, momCorrSkipMask] = update_transfer_correction_masks( ...
            momCorrBulkBulkMask, momCorrSkipMask, srcWet, recvWet, false(size(cntVec)), true);
        movedOut = movedOut + movedWet;
        movedIn  = movedIn  + movedWet;
    end

    % Réévaluation des seuils avant les cellules peu peuplées.
    Nmat = reshape(cntVec, [Ny, Nx]);
    Nsm  = smooth_occupancy_field(Nmat);
    [~, ~, lowThrLocMap, highThrLocMap] = ...
        build_local_fluid_fraction_threshold_maps(Nsm, params);

    bulkMask = build_topological_bulk_mask(cntVec, params) & activeMaskVec;
    if ~flags.enableSurfaceTopology
        bulkMask = (cntVec > 0) & activeMaskVec;
    end
    lowThrActiveMap  = lowThrLocMap;
    highThrActiveMap = highThrLocMap;
    if ~flags.useActiveVolumeTargets
        lowThrActiveMap(bulkMask)  = lowThrBulk;
        highThrActiveMap(bulkMask) = highThr;
    end
    bulkMaskLowBefore = bulkMask;

    stats(1) = stats(1) + nHigh;
    stats(2) = stats(2) + nLow;
    stats(3) = stats(3) + movedOut;
    stats(4) = stats(4) + movedIn;

    Nmat = reshape(cntVec, [Ny, Nx]);
    Nsm  = smooth_occupancy_field(Nmat);
    [~, ~, lowThrLocMap, highThrLocMap] = ...
        build_local_fluid_fraction_threshold_maps(Nsm, params);

    bulkMask = build_topological_bulk_mask(cntVec, params) & activeMaskVec;
    if ~flags.enableSurfaceTopology
        bulkMask = (cntVec > 0) & activeMaskVec;
    end
    highThrActiveMap = highThrLocMap;
    if ~flags.useActiveVolumeTargets
        highThrActiveMap(bulkMask) = highThr;
    end

    highThrVec = reshape(highThrActiveMap, [], 1);
    denseRemain = find(activeMaskVec & (cntVec > highThrVec));

    if ~isempty(denseRemain) && (pass >= params.maxRedistribPasses)
        fprintf('ALERTE saturation redistribution : pass=%d, nDenseRemain=%d\n', ...
            pass, numel(denseRemain));
    end

    keepRedistributing = ~isempty(denseRemain) && (movedOut > 0); %#ok<NASGU>
end

% ================================================================
% 2) Cellules peu peuplées : bulk / interface traités via des drapeaux
% explicites et les opérateurs de bord.
% ================================================================
lowThrVec = reshape(lowThrActiveMap, [], 1);
lowList = find(activeMaskVec & cntVec > 0 & cntVec < lowThrVec);
if ~isempty(lowList)
    lowList = lowList(randperm(numel(lowList)));
end

nHigh = 0; nLow = 0; movedOut = 0; movedIn = 0;
for ii = 1:numel(lowList)
    c = lowList(ii);
    nc = cntVec(c);
    if nc <= 0, continue; end
    nLow = nLow + 1;

    idsSrc = partsInCell{c};
    if isempty(idsSrc), continue; end

    neigh = neighborsCore{c};
    if isempty(neigh), continue; end

    isSurface = flags.enableSurfaceTopology && ~is_bulk_like_cell(c, cntVec, params);
    [iyc, ixc] = ind2sub([Ny, Nx], c);
    if flags.useActiveVolumeTargets
        lowThrLocal = lowThrActiveMap(iyc, ixc);
    elseif isSurface
        lowThrLocal = lowThrLocMap(iyc, ixc);
    else
        lowThrLocal = lowThrBulk;
    end
    if nc >= lowThrLocal
        continue;
    end

    if isSurface
        if flags.preserveWettingLayer && is_cell_adjacent_to_wetting_boundary(c, params, flags)
            continue;
        end

        wallAdj = is_wall_adjacent_cell(c, params);
        if wallAdj
            dir = wall_normal_inward_unit(c, params);
        else
            dir = inward_direction_from_local_gradient(c, Nsm, params);
        end

        if wallAdj
            admiss = admissible_wall_surface_targets_range2(c, cntVec, highThrActiveMap, params);
        else
            admiss = admissible_inward_targets_range2(c, cntVec, highThrActiveMap, params);
        end
        if isempty(admiss), continue; end

        if wallAdj
            ordered = rank_wall_surface_targets_range2(c, admiss, dir, cntVec, Nsm, params);
        elseif has_only_diagonal_support(c, cntVec, params)
            ordered = rank_admissible_neighbors_inward_range2(c, admiss, 2.0*dir, Nsm, params);
        else
            ordered = rank_admissible_neighbors_inward_range2(c, admiss, dir, Nsm, params);
        end

        if wallAdj
            targetMapSurf = reshape(highThrActiveMap, [], 1);
            needMove = numel(idsSrc);
            totalCap = 0;
            for kk = 1:numel(ordered)
                nb = ordered(kk);
                totalCap = totalCap + max(0, floor(targetMapSurf(nb) - cntVec(nb)));
            end
            if totalCap < needMove
                continue;
            end

            [x, partsInCell, cntVec, movedNow, recvTargets] = move_all_particles_to_targets_round_robin( ...
                x, partsInCell, cntVec, c, ordered, targetMapSurf, params);
            if ~isempty(recvTargets), recvMask(recvTargets) = true; end
            [momCorrBulkBulkMask, momCorrSkipMask] = update_transfer_correction_masks( ...
                momCorrBulkBulkMask, momCorrSkipMask, c, recvTargets, bulkMaskLowBefore, true);
            if movedNow ~= needMove
                error('Surface wall redistribution incomplete despite sufficient cumulative capacity.');
            end
            movedOut = movedOut + movedNow;
            movedIn  = movedIn  + movedNow;
        else
            target = ordered(1);

            nt = numel(idsSrc);
            idsMove = idsSrc;

            x(idsMove,:) = sample_in_cell(target, nt, params);

            partsInCell{c} = zeros(0,1);
            cntVec(c) = 0;

            partsInCell{target} = [partsInCell{target}; idsMove(:)];
            cntVec(target) = cntVec(target) + nt;
            recvMask(target) = true;
            [momCorrBulkBulkMask, momCorrSkipMask] = update_transfer_correction_masks( ...
                momCorrBulkBulkMask, momCorrSkipMask, c, target, bulkMaskLowBefore, true);

            movedOut = movedOut + nt;
            movedIn  = movedIn  + nt;
        end

    else
        lowThrBulkInt = ceil(lowThrBulk);
        need = max(0, lowThrBulkInt - nc);
        if need <= 0, continue; end

        donorMask = false(numel(neigh),1);
        for jj = 1:numel(neigh)
            dtest = neigh(jj);
            donorMask(jj) = (cntVec(dtest) > lowThrBulkInt) && ...
                            ((~flags.enableSurfaceTopology) || is_bulk_like_cell(dtest, cntVec, params));
        end
        donors = neigh(donorMask);
        if isempty(donors), continue; end

        dir = inward_direction_from_local_gradient(c, Nsm, params);
        orderedDonors = rank_admissible_neighbors_inward(c, donors, dir, Nsm, params);

        for kk = 1:numel(orderedDonors)
            d = orderedDonors(kk);
            idsDonor = partsInCell{d};
            if isempty(idsDonor), continue; end
            if flags.enableSurfaceTopology && ~is_bulk_like_cell(d, cntVec, params), continue; end

            avail = cntVec(d) - lowThrBulkInt;
            if avail <= 0, continue; end

            nt = min([need, avail, numel(idsDonor)]);
            if nt <= 0, continue; end

            permTake = randperm(numel(idsDonor), nt);
            idsMove = idsDonor(permTake);

            x(idsMove,:) = sample_in_cell(c, nt, params);

            keepMask = true(numel(idsDonor),1);
            keepMask(permTake) = false;
            partsInCell{d} = idsDonor(keepMask);
            cntVec(d) = numel(partsInCell{d});

            partsInCell{c} = [partsInCell{c}; idsMove(:)];
            cntVec(c) = cntVec(c) + nt;
            recvMask(c) = true;
            [momCorrBulkBulkMask, momCorrSkipMask] = update_transfer_correction_masks( ...
                momCorrBulkBulkMask, momCorrSkipMask, d, c, bulkMaskLowBefore, false);

            movedOut = movedOut + nt;
            movedIn  = movedIn  + nt;
            need = need - nt;

            if need <= 0, break; end
        end
    end
end

stats(1) = stats(1) + nHigh;
stats(2) = stats(2) + nLow;
stats(3) = stats(3) + movedOut;
stats(4) = stats(4) + movedIn;

Nmat = reshape(cntVec, [Ny, Nx]);
Nsm  = smooth_occupancy_field(Nmat);
[~, ~, lowThrLocMap, highThrLocMap] = ...
    build_local_fluid_fraction_threshold_maps(Nsm, params);

bulkMask = build_topological_bulk_mask(cntVec, params) & activeMaskVec;
if ~flags.enableSurfaceTopology
    bulkMask = (cntVec > 0) & activeMaskVec;
end
highThrActiveMap = highThrLocMap;
if ~flags.useActiveVolumeTargets
    highThrActiveMap(bulkMask) = highThr;
end

highThrVec = reshape(highThrActiveMap, [], 1);
denseRemain = find(activeMaskVec & (cntVec > highThrVec)); %#ok<NASGU>

dU2 = zeros(Nx*Ny,1);
dUmax = 0;
pxErr2 = zeros(Nx*Ny,1);
pyErr2 = zeros(Nx*Ny,1);
uxNow = zeros(Nx*Ny,1);
uyNow = zeros(Nx*Ny,1);

for c = 1:(Nx*Ny)
    ids = partsInCell{c};
    nc = numel(ids);
    if nc == 0, continue; end

    Pafter = sum(v(ids,:), 1);
    uNow = Pafter / nc;
    uxNow(c) = uNow(1);
    uyNow(c) = uNow(2);

    dP = [PbeforeX(c), PbeforeY(c)] - Pafter;
    du = dP / nc;

    pxErr2(c) = dP(1)^2;
    pyErr2(c) = dP(2)^2;
    dU2(c) = du(1)^2 + du(2)^2;
    dUmax = max(dUmax, sqrt(dU2(c)));

    doCellMomentumCorrection = enableMomentumCorrection;
    if doCellMomentumCorrection && momCorrSkipMask(c)
        doCellMomentumCorrection = false;
    end
    if doCellMomentumCorrection && recvMask(c) && ~momCorrBulkBulkMask(c)
        doCellMomentumCorrection = false;
    end

    if doCellMomentumCorrection && (du(1) ~= 0 || du(2) ~= 0)
        v(ids,1) = v(ids,1) + du(1);
        v(ids,2) = v(ids,2) + du(2);
    end
end

nzdu = dU2 > 0;
if any(nzdu)
    stats(5) = sqrt(mean(dU2(nzdu)));
    stats(6) = dUmax;
end

Ncell = reshape(cntVec, [Ny,Nx]);
PglobAfter = sum(v,1);
Eafter = 0.5 * sum(v(:,1).^2 + v(:,2).^2);
N = Ncell(:);

Nmat = reshape(cntVec, [Ny, Nx]);
Nsm  = smooth_occupancy_field(Nmat);
[~, ~, lowThrLocMap, highThrLocMap] = ...
    build_local_fluid_fraction_threshold_maps(Nsm, params);

Nact = N(activeMaskVec);
if isempty(Nact)
    diag(1) = 0;
else
    diag(1) = std(double(Nact))/max(gamma, eps);
end
bulkMask = build_topological_bulk_mask(cntVec, params) & activeMaskVec;
if ~flags.enableSurfaceTopology
    bulkMask = (cntVec > 0) & activeMaskVec;
end
lowThrActiveMap  = lowThrLocMap;
highThrActiveMap = highThrLocMap;
if ~flags.useActiveVolumeTargets
    lowThrActiveMap(bulkMask)  = lowThrBulk;
    highThrActiveMap(bulkMask) = highThr;
end
outBandMask = activeMaskVec & (double(N(:)) > reshape(highThrActiveMap, [], 1) | ...
              (double(N(:)) > 0 & double(N(:)) < reshape(lowThrActiveMap, [], 1)));
diag(2) = mean(outBandMask(activeMaskVec));
diag(3) = norm(PglobAfter - PglobBefore);
diag(4) = sqrt(mean(pxErr2(pxErr2>0)));
diag(5) = sqrt(mean(pyErr2(pyErr2>0)));
if ~isfinite(diag(4)), diag(4) = 0; end
if ~isfinite(diag(5)), diag(5) = 0; end
diag(6) = Eafter - Ebefore;
diag(8) = sqrt(mean(uxNow(N>0).^2));
diag(9) = sqrt(mean(uyNow(N>0).^2));
if any(activeMaskVec)
    diag(10) = mean(N(activeMaskVec) == 0);
else
    diag(10) = 0;
end

if doDivDiag
    [Ux, Uy, ~, ~] = velocity_field_and_vorticity(x, v, params);
    [dUx_dx, ~] = gradient(Ux, Lx/Nx, Ly/Ny);
    [~, dUy_dy] = gradient(Uy, Lx/Nx, Ly/Ny);
    div = dUx_dx + dUy_dy;
    diag(7) = sqrt(mean(div(:).^2));
else
    diag(7) = NaN;
end
end

function nb = preferred_neighbor_from_local_gradient(c, Nsm, mode, params)
Nx = params.Nx;
Ny = params.Ny;
aX = params.Lx / params.Nx;
aY = params.Ly / params.Ny;
neighborsCore = params.neighborsCore;
[iy, ix] = ind2sub([Ny,Nx], c);
neigh = neighborsCore{c};
if isempty(neigh)
    nb = -1;
    return;
end
neigh = neigh(randperm(numel(neigh)));

dNdx = local_diff_x(Nsm, iy, ix, params);
dNdy = local_diff_y(Nsm, iy, ix, params);

if strcmpi(mode, 'outward')
    dir = [-dNdx, -dNdy];
else
    dir = [dNdx, dNdy];
end

ndir = hypot(dir(1), dir(2));
if ndir < 1e-14
    occ = zeros(numel(neigh),1);
    for k = 1:numel(neigh)
        [jy, jx] = ind2sub([Ny,Nx], neigh(k));
        occ(k) = Nsm(jy, jx);
    end
    if strcmpi(mode, 'outward')
        [~, kbest] = min(occ);
    else
        [~, kbest] = max(occ);
    end
    nb = neigh(kbest);
    return;
end

x0 = (ix - 0.5)*aX;
y0 = (iy - 0.5)*aY;

bestScore = -inf;
nb = neigh(1);
for kk = 1:numel(neigh)
    cn = neigh(kk);
    [jyn, jxn] = ind2sub([Ny,Nx], cn);
    [dxCell, dyCell] = minimal_cell_step_general(ix, iy, jxn, jyn, params);
    step = [dxCell*aX, dyCell*aY];
    nstep = hypot(step(1), step(2));
    if nstep < 1e-14
        continue;
    end
    score = (step(1)*dir(1) + step(2)*dir(2)) / (ndir*nstep);
    score = score + 1e-12*randn();
    if score > bestScore
        bestScore = score;
        nb = cn;
    end
end
end

function dir = inward_direction_from_local_gradient(c, Nsm, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny,Nx], c);
dNdx = local_diff_x(Nsm, iy, ix, params);
dNdy = local_diff_y(Nsm, iy, ix, params);
dir = [dNdx, dNdy];
end

function ordered = rank_admissible_neighbors_inward(c, admiss, dir, Nsm, params)
Nx = params.Nx;
Ny = params.Ny;
if isempty(admiss)
    ordered = zeros(0,1);
    return;
end

[iy, ix] = ind2sub([Ny,Nx], c);
perm = randperm(numel(admiss));
admiss = admiss(perm);

ndir = hypot(dir(1), dir(2));
scoreAlign = zeros(numel(admiss),1);
scoreDense = zeros(numel(admiss),1);

for k = 1:numel(admiss)
    cn = admiss(k);
    [jy, jx] = ind2sub([Ny,Nx], cn);
    [dx, dy] = minimal_cell_step_general(ix, iy, jx, jy, params);
    step = [dx, dy];
    nstep = hypot(step(1), step(2));
    if ndir < 1e-14 || nstep < 1e-14
        scoreAlign(k) = 0;
    else
        scoreAlign(k) = (step(1)*dir(1) + step(2)*dir(2)) / (ndir*nstep);
    end
    scoreDense(k) = Nsm(jy, jx);
end

score = 100*scoreAlign + scoreDense;
score = score + 1e-9*randn(size(score));
[~, ord] = sort(score, 'descend');
ordered = admiss(ord);
end

function ordered = rank_admissible_neighbors_inward_range2(c, admiss, dir, Nsm, params)
Nx = params.Nx;
Ny = params.Ny;
if isempty(admiss)
    ordered = zeros(0,1);
    return;
end

[iy, ix] = ind2sub([Ny,Nx], c);
perm = randperm(numel(admiss));
admiss = admiss(perm);

ndir = hypot(dir(1), dir(2));
scoreAlign = zeros(numel(admiss),1);
scoreDense = zeros(numel(admiss),1);
scoreRange = zeros(numel(admiss),1);

for k = 1:numel(admiss)
    cn = admiss(k);
    [jy, jx] = ind2sub([Ny,Nx], cn);
    [dx, dy] = minimal_cell_step_general(ix, iy, jx, jy, params);
    step = [dx, dy];
    cheb = max(abs(step(1)), abs(step(2)));
    nstep = hypot(step(1), step(2));

    if ndir < 1e-14 || nstep < 1e-14
        scoreAlign(k) = 0;
    else
        scoreAlign(k) = (step(1)*dir(1) + step(2)*dir(2)) / (ndir*nstep);
    end
    scoreDense(k) = Nsm(jy, jx);

    if cheb <= 1
        scoreRange(k) = 0.0;
    else
        scoreRange(k) = -0.2;
    end
end

score = 100*scoreAlign + scoreDense + scoreRange;
score = score + 1e-9*randn(size(score));
[~, ord] = sort(score, 'descend');
ordered = admiss(ord);
end

function admiss = admissible_inward_targets_range2(c, cntVec, highThrLocMap, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny,Nx], c);

lst = zeros(24,1);
m = 0;

for dy = -2:2
    for dx = -2:2
        if dx == 0 && dy == 0
            continue;
        end
        [ok, jy, jx] = shift_index_general(iy, ix, dy, dx, params);
        if ~ok
            continue;
        end
        cn = sub2ind([Ny,Nx], jy, jx);
        highThrNb = highThrLocMap(jy, jx);
        if cntVec(cn) < highThrNb
            m = m + 1;
            lst(m) = cn;
        end
    end
end

admiss = unique(lst(1:m).', 'stable');
if ~isempty(admiss)
    admiss = admiss(randperm(numel(admiss)));
end
end

function tf = is_wall_adjacent_cell(c, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

tf = ((ix == 1)  && is_physical_boundary_side(params, 'left'))  || ...
     ((ix == Nx) && is_physical_boundary_side(params, 'right')) || ...
     ((iy == 1)  && is_physical_boundary_side(params, 'bottom')) || ...
     ((iy == Ny) && is_physical_boundary_side(params, 'top'));
end

function wettingTarget = resolve_bottom_wetting_target(lowThrBulk, params)
if isfield(params,'bottomWettingTargetOverride') && isfinite(params.bottomWettingTargetOverride)
    wettingTarget = ceil(params.bottomWettingTargetOverride);
else
    wettingTarget = ceil(lowThrBulk);
end
wettingTarget = max(1, wettingTarget);
end

function [x, partsInCell, cntVec, movedCount, recvTargets, srcTouched] = apply_bottom_wall_wetting_layer( ...
    x, partsInCell, cntVec, wettingTarget, params)

Nx = params.Nx;
Ny = params.Ny;
movedCount = 0;
recvTargets = zeros(0,1);
srcTouched = zeros(0,1);

if Ny < 1 || wettingTarget <= 0
    return;
end

bottomCells = zeros(1, Nx);
for jx = 1:Nx
    bottomCells(jx) = sub2ind([Ny, Nx], 1, jx);
end

sourceList = bottomCells(cntVec(bottomCells) > wettingTarget);
if isempty(sourceList)
    return;
end
sourceList = sourceList(randperm(numel(sourceList)));

targetMap = zeros(Ny, Nx);
targetMap(1, :) = wettingTarget;

for kk = 1:numel(sourceList)
    c = sourceList(kk);
    if cntVec(c) <= wettingTarget
        continue;
    end

    needOut = cntVec(c) - wettingTarget;
    admiss = admissible_bottom_wetting_targets(c, cntVec, wettingTarget, params);
    if isempty(admiss)
        continue;
    end

    ordered = rank_bottom_wetting_targets(c, admiss, params);
    [x, partsInCell, cntVec, movedNow, recvNow] = spread_particles_from_dense_source_round_robin( ...
        x, partsInCell, cntVec, c, ordered, targetMap, needOut, params);

    movedCount = movedCount + movedNow;
    if movedNow > 0
        srcTouched = unique([srcTouched(:); c]);
    end
    if ~isempty(recvNow)
        recvTargets = unique([recvTargets(:); recvNow(:)]).';
    end
end
end

function admiss = admissible_bottom_wetting_targets(c, cntVec, wettingTarget, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

if iy ~= 1
    admiss = zeros(0,1);
    return;
end

if isfield(params,'bottomWettingMaxRange') && isfinite(params.bottomWettingMaxRange)
    maxRange = max(1, floor(params.bottomWettingMaxRange));
else
    maxRange = Nx;
end

lst = zeros(1, 2*maxRange);
m = 0;
for dj = 1:maxRange
    jL = ix - dj;
    if jL >= 1
        nb = sub2ind([Ny, Nx], 1, jL);
        if cntVec(nb) < wettingTarget
            m = m + 1;
            lst(m) = nb;
        end
    end
    jR = ix + dj;
    if jR <= Nx
        nb = sub2ind([Ny, Nx], 1, jR);
        if cntVec(nb) < wettingTarget
            m = m + 1;
            lst(m) = nb;
        end
    end
end
admiss = lst(1:m);
end

function ordered = rank_bottom_wetting_targets(c, admiss, params)
Nx = params.Nx;
Ny = params.Ny;
[~, ix] = ind2sub([Ny, Nx], c);

if isempty(admiss)
    ordered = zeros(0,1);
    return;
end

jxAll = zeros(size(admiss));
for k = 1:numel(admiss)
    [~, jx] = ind2sub([Ny, Nx], admiss(k));
    jxAll(k) = jx;
end

dAll = abs(jxAll - ix);
dUnique = unique(dAll(:)).';
ordered = zeros(1, 0);

for d = dUnique
    mask = (dAll == d);
    cand = admiss(mask);
    if numel(cand) > 1
        cand = cand(randperm(numel(cand)));
    end
    ordered = [ordered, cand]; %#ok<AGROW>
end
end

function admiss = admissible_wall_surface_targets_range2(c, cntVec, highThrLocMap, params)
Nx = params.Nx;
Ny = params.Ny;
% Pour résorber une cellule pauvre collée à une vraie paroi (non périodique)
% sans créer de glissement tangent, on n'autorise que des cibles situées
% strictement vers l'intérieur du domaine relativement aux côtés physiques
% effectivement présents.
[iy, ix] = ind2sub([Ny, Nx], c);

leftWall   = (ix == 1)  && is_physical_boundary_side(params, 'left');
rightWall  = (ix == Nx) && is_physical_boundary_side(params, 'right');
bottomWall = (iy == 1)  && is_physical_boundary_side(params, 'bottom');
topWall    = (iy == Ny) && is_physical_boundary_side(params, 'top');

lst = zeros(24,1);
m = 0;

for dy = -2:2
    for dx = -2:2
        if dx == 0 && dy == 0
            continue;
        end
        [ok, jy, jx] = shift_index_general(iy, ix, dy, dx, params);
        if ~ok
            continue;
        end

        inwardOk = true;
        if bottomWall, inwardOk = inwardOk && (jy >= iy + 1); end
        if topWall,    inwardOk = inwardOk && (jy <= iy - 1); end
        if leftWall,   inwardOk = inwardOk && (jx >= ix + 1); end
        if rightWall,  inwardOk = inwardOk && (jx <= ix - 1); end
        if ~inwardOk
            continue;
        end

        cn = sub2ind([Ny, Nx], jy, jx);
        highThrNb = highThrLocMap(jy, jx);
        freeCap = floor(highThrNb - cntVec(cn));
        if freeCap <= 0
            continue;
        end

        m = m + 1;
        lst(m) = cn;
    end
end

admiss = lst(1:m).';
end

function dir = wall_normal_inward_unit(c, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

dir = [0, 0];
if (iy == 1)  && is_physical_boundary_side(params, 'bottom'), dir = dir + [0, 1]; end
if (iy == Ny) && is_physical_boundary_side(params, 'top'),    dir = dir + [0, -1]; end
if (ix == 1)  && is_physical_boundary_side(params, 'left'),   dir = dir + [1, 0]; end
if (ix == Nx) && is_physical_boundary_side(params, 'right'),  dir = dir + [-1, 0]; end

ndir = hypot(dir(1), dir(2));
if ndir > 0
    dir = dir / ndir;
end
end

function ordered = rank_wall_surface_targets_range2(c, admiss, dir, cntVec, Nsm, params)
% Classement symétrique au mur :
% 1) progression vers le bulk
% 2) écart tangentiel minimal en valeur absolue
% 3) cellules les moins occupées / les plus capacitaires
if isempty(admiss)
    ordered = zeros(0,1);
    return;
end

Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

score = zeros(numel(admiss),1);
for k = 1:numel(admiss)
    cn = admiss(k);
    [jy, jx] = ind2sub([Ny, Nx], cn);
    [dx, dy] = minimal_cell_step_general(ix, iy, jx, jy, params);
    step = [dx, dy];

    inwardScore = 0;
    nstep = hypot(step(1), step(2));
    ndir = hypot(dir(1), dir(2));
    if ndir > 0 && nstep > 0
        inwardScore = (step(1)*dir(1) + step(2)*dir(2)) / (ndir*nstep);
    end

    tangentialPenalty = abs(step(1)*dir(2) - step(2)*dir(1));
    score(k) = 100*inwardScore - 5*tangentialPenalty - 0.1*cntVec(cn) + 0.01*Nsm(jy,jx) + 1e-9*randn();
end

[~, ord] = sort(score, 'descend');
ordered = admiss(ord);
end

function [x, partsInCell, cntVec, movedCount, recvTargets] = move_all_particles_to_targets_round_robin( ...
    x, partsInCell, cntVec, c, orderedTargets, targetMap, params)

movedCount = 0;
recvTargets = zeros(0,1);
if isempty(orderedTargets)
    return;
end

idsAvail = partsInCell{c};
needOut = numel(idsAvail);
if needOut <= 0
    return;
end

totalCap = 0;
for kk = 1:numel(orderedTargets)
    nb = orderedTargets(kk);
    totalCap = totalCap + max(0, floor(targetMap(nb) - cntVec(nb)));
end
if totalCap < needOut
    return;
end

while needOut > 0 && ~isempty(idsAvail)
    progressed = false;

    for kk = 1:numel(orderedTargets)
        nb = orderedTargets(kk);
        deficit = max(0, floor(targetMap(nb) - cntVec(nb)));
        if deficit <= 0
            continue;
        end

        remainingTargets = 0;
        for jj = kk:numel(orderedTargets)
            nbj = orderedTargets(jj);
            defj = max(0, floor(targetMap(nbj) - cntVec(nbj)));
            if defj > 0
                remainingTargets = remainingTargets + 1;
            end
        end

        quota = max(1, ceil(needOut / max(1, remainingTargets)));
        nt = min([quota, deficit, needOut, numel(idsAvail)]);
        if nt <= 0
            continue;
        end

        idsMove = idsAvail(1:nt);
        idsAvail(1:nt) = [];

        x(idsMove,:) = sample_in_cell(nb, nt, params);
        partsInCell{nb} = [partsInCell{nb}; idsMove(:)];
        cntVec(nb) = cntVec(nb) + nt;
        recvTargets = [recvTargets; repmat(nb, nt>0, 1)]; %#ok<AGROW>

        movedCount = movedCount + nt;
        needOut = needOut - nt;
        progressed = true;

        if needOut <= 0 || isempty(idsAvail)
            break;
        end
    end

    if ~progressed
        break;
    end
end

if movedCount > 0
    partsInCell{c} = idsAvail;
    cntVec(c) = numel(idsAvail);
    recvTargets = unique(recvTargets);
else
    recvTargets = zeros(0,1);
end
end

function [bulkBulkMask, skipMask] = update_transfer_correction_masks( ...
    bulkBulkMask, skipMask, srcCells, recvCells, bulkMaskBefore, forceSkip)

srcCells = unique(srcCells(:));
recvCells = unique(recvCells(:));
touched = unique([srcCells; recvCells]);
if isempty(touched)
    return;
end

if forceSkip
    skipMask(touched) = true;
    bulkBulkMask(touched) = false;
    return;
end

allBulk = true;
for k = 1:numel(touched)
    ck = touched(k);
    if ck < 1 || ck > numel(bulkMaskBefore) || ~bulkMaskBefore(ck)
        allBulk = false;
        break;
    end
end

if allBulk
    bulkBulkMask(touched) = true;
else
    skipMask(touched) = true;
    bulkBulkMask(touched) = false;
end
end

function v = apply_interface_velocity_reorientation_step(v, x, params)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
% Réorientation douce des vitesses des particules d'interface vers l'intérieur,
% avec cartes interfaciales cohérentes avec le critère réel d'application.

if params.interfaceReorientBeta <= 0 || params.interfaceReorientFrac <= 0
    return;
end

[cntVec, partsInCell] = rebuild_parts_in_cell(x, params);
[~, ~, ~, ~, betaMap, fracMap, maskInterface] = build_interface_reorientation_maps( ...
    x, params);

Nmat = reshape(cntVec, [Ny, Nx]);
Nsm  = smooth_occupancy_field(Nmat);

moveIds = [];
vOldSel = zeros(0,2);
vTarSel = zeros(0,2);
betaSel = zeros(0,1);

nCellsOccupied  = 0;
nCellsInterface = 0;
nCellsMoved     = 0;
nPartsMoved     = 0;

for c = 1:(Nx*Ny)
    if cntVec(c) <= 0
        continue;
    end
    nCellsOccupied = nCellsOccupied + 1;

    [iy, ix] = ind2sub([Ny, Nx], c);
    if ~maskInterface(iy, ix)
        continue;
    end
    nCellsInterface = nCellsInterface + 1;

    ids = partsInCell{c};
    if isempty(ids)
        continue;
    end

    betaLoc = betaMap(iy, ix);
    fracLoc = fracMap(iy, ix);
    if betaLoc <= 0 || fracLoc <= 0
        continue;
    end

    dir = inward_direction_from_local_gradient(c, Nsm, params);
    nd = hypot(dir(1), dir(2));
    if nd < 1e-14
        continue;
    end
    dir = dir / nd;

    np = numel(ids);
    nm = max(1, round(fracLoc * np));
    nm = min(nm, np);
    if nm <= 0
        continue;
    end

    pick = randperm(np, nm);
    idsMove = ids(pick);

    vCell = v(idsMove, :);
    speed = sqrt(sum(vCell.^2, 2));
    vTarget = [speed * dir(1), speed * dir(2)];

    moveIds = [moveIds; idsMove(:)]; %#ok<AGROW>
    vOldSel = [vOldSel; vCell]; %#ok<AGROW>
    vTarSel = [vTarSel; vTarget]; %#ok<AGROW>
    betaSel = [betaSel; betaLoc * ones(nm, 1)]; %#ok<AGROW>

    nCellsMoved = nCellsMoved + 1;
    nPartsMoved = nPartsMoved + nm;
end

if isempty(moveIds)
    if isfield(params, 'interfaceReorientVerbose') && params.interfaceReorientVerbose
        fprintf('interface reorient: occupied=%d interface=%d movedCells=%d movedParts=%d', ...
            nCellsOccupied, nCellsInterface, nCellsMoved, nPartsMoved);
    end
    return;
end

vNewSel = zeros(size(vOldSel));
vNewSel(:,1) = (1 - betaSel) .* vOldSel(:,1) + betaSel .* vTarSel(:,1);
vNewSel(:,2) = (1 - betaSel) .* vOldSel(:,2) + betaSel .* vTarSel(:,2);

dV = vNewSel - vOldSel;
meanDV = mean(dV, 1);
dV = dV - meanDV;
vNewSel = vOldSel + dV;

s0 = sqrt(sum(vOldSel.^2, 2));
s1 = sqrt(sum(vNewSel.^2, 2));
ok = s1 > 1e-14;
vNewSel(ok,1) = vNewSel(ok,1) .* (s0(ok) ./ s1(ok));
vNewSel(ok,2) = vNewSel(ok,2) .* (s0(ok) ./ s1(ok));

bad = ~ok & (s0 > 0);
if any(bad)
    vNewSel(bad,:) = vTarSel(bad,:);
end

v(moveIds, :) = vNewSel;

if isfield(params, 'interfaceReorientVerbose') && params.interfaceReorientVerbose
    fprintf('interface reorient: occupied=%d interface=%d movedCells=%d movedParts=%d', ...
        nCellsOccupied, nCellsInterface, nCellsMoved, nPartsMoved);
end
end


function plot_interface_reorientation_diag(x, it, params)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
aX = Lx / Nx;
aY = Ly / Ny;

[Nsm, nx, ny, angMap, betaMap, fracMap, maskInterface] = ...
    build_interface_reorientation_maps(x, params);

xc = ((0:Nx-1)+0.5) * aX;
yc = ((0:Ny-1)+0.5) * aY;
[Xc, Yc] = meshgrid(xc, yc);

figure(3); clf;

subplot(2,3,1);
imagesc(xc, yc, Nsm); set(gca,'YDir','normal'); hold on;
stepQ = max(1, round(min(Nx,Ny)/30));
qmask = false(size(maskInterface));
qmask(1:stepQ:end, 1:stepQ:end) = true;
qmask = qmask & maskInterface;
quiver(Xc(qmask), Yc(qmask), nx(qmask), ny(qmask), 0.4, 'k');
axis equal tight; colorbar;
xlabel('x'); ylabel('y');
title(sprintf('Normales interface | it=%d', it));

subplot(2,3,2);
imagesc(xc, yc, angMap); set(gca,'YDir','normal');
axis equal tight; colorbar;
xlabel('x'); ylabel('y');
title('angMis');

subplot(2,3,3);
imagesc(xc, yc, double(maskInterface)); set(gca,'YDir','normal');
axis equal tight; colorbar;
xlabel('x'); ylabel('y');
title('Masque interface');

subplot(2,3,4);
imagesc(xc, yc, betaMap); set(gca,'YDir','normal');
axis equal tight; colorbar;
xlabel('x'); ylabel('y');
title(sprintf('beta local | base=%.3g', params.interfaceReorientBeta));

subplot(2,3,5);
imagesc(xc, yc, fracMap); set(gca,'YDir','normal');
axis equal tight; colorbar;
xlabel('x'); ylabel('y');
title(sprintf('frac locale | base=%.3g', params.interfaceReorientFrac));

subplot(2,3,6);
imagesc(xc, yc, angMap .* maskInterface); set(gca,'YDir','normal');
axis equal tight; colorbar;
xlabel('x'); ylabel('y');
title('angMis sur interface');
end


function [Nsm, nx, ny, angMap, betaMap, fracMap, maskInterface] = ...
    build_interface_reorientation_maps(x, params)

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;

[cntVec, ~] = rebuild_parts_in_cell(x, params);
Nmat = reshape(cntVec, [Ny, Nx]);
Nsm  = smooth_occupancy_field(Nmat);

aX = Lx / Nx;
aY = Ly / Ny;
[dNdx, dNdy] = gradient(Nsm, aX, aY);
ng = sqrt(dNdx.^2 + dNdy.^2) + 1e-12;
nx = dNdx ./ ng;
ny = dNdy ./ ng;

maskInterface = false(Ny, Nx);
angMap = nan(Ny, Nx);
betaMap = nan(Ny, Nx);
fracMap = nan(Ny, Nx);

for c = 1:(Nx*Ny)
    if cntVec(c) <= 0
        continue;
    end
    if is_bulk_like_cell(c, cntVec, params)
        continue;
    end

    [iy, ix] = ind2sub([Ny, Nx], c);
    maskInterface(iy, ix) = true;

    angMis = local_normal_misalignment(c, nx, ny, cntVec, params);
    angMap(iy, ix) = angMis;

    if params.useCurvatureWeightedReorientation
        angEff = min(angMis, params.curvatureAngleMax) ./ max(params.curvatureAngleMax, 1e-12);
        betaLoc = params.interfaceReorientBeta .* (1 + params.curvatureReorientBetaGain .* angEff);
        fracLoc = params.interfaceReorientFrac .* (1 + params.curvatureReorientFracGain .* angEff);
        betaLoc = min(betaLoc, params.interfaceReorientBetaMax);
        fracLoc = min(fracLoc, params.interfaceReorientFracMax);
    else
        betaLoc = params.interfaceReorientBeta;
        fracLoc = params.interfaceReorientFrac;
    end

    betaMap(iy, ix) = betaLoc;
    fracMap(iy, ix) = fracLoc;
end
end


function [nHat, inwardWeight] = dense_cell_inward_bias(c, cntVec, Nsm, params)
if is_bulk_like_cell(c, cntVec, params)
    nHat = [0, 0];
    inwardWeight = 0.0;
    return;
end

if is_wall_adjacent_cell(c, params)
    nHat = wall_normal_inward_unit(c, params);
else
    nHat = inward_direction_from_local_gradient(c, Nsm, params);
end
nNorm = hypot(nHat(1), nHat(2));
if nNorm > 1e-14
    nHat = nHat / nNorm;
else
    nHat = [0, 0];
end
inwardWeight = 0.75;
end

function [x, partsInCell, cntVec, movedCount, recvTargets] = redistribute_dense_cell_at_range( ...
    x, partsInCell, cntVec, c, lowThrActiveMap, highThrActiveMap, ...
    highThrLocal, nHat, inwardWeight, params, rangeR)

movedCount = 0;
recvTargets = zeros(0,1);
needOut = max(0, ceil(cntVec(c) - highThrLocal));
if needOut <= 0
    return;
end

poorTargets = admissible_dense_targets_rangeR(c, cntVec, lowThrActiveMap, params, rangeR);
orderedPoor = rank_dense_targets_unified(c, poorTargets, cntVec, lowThrActiveMap, nHat, inwardWeight, params);
[x, partsInCell, cntVec, movedNow, recvNow] = spread_particles_from_dense_source_round_robin( ...
    x, partsInCell, cntVec, c, orderedPoor, lowThrActiveMap, needOut, params);
if ~isempty(recvNow)
    recvTargets = unique([recvTargets(:); recvNow(:)]).';
end
movedCount = movedCount + movedNow;

remainingOut = max(0, ceil(cntVec(c) - highThrLocal));
if remainingOut <= 0
    return;
end

reliefTargets = admissible_dense_targets_rangeR(c, cntVec, highThrActiveMap, params, rangeR);
orderedRelief = rank_dense_targets_unified(c, reliefTargets, cntVec, highThrActiveMap, nHat, inwardWeight, params);
[x, partsInCell, cntVec, movedNow, recvNow] = spread_particles_from_dense_source_round_robin( ...
    x, partsInCell, cntVec, c, orderedRelief, highThrActiveMap, remainingOut, params);
if ~isempty(recvNow)
    recvTargets = unique([recvTargets(:); recvNow(:)]).';
end
movedCount = movedCount + movedNow;
end

function [x, partsInCell, cntVec, movedCount, recvTargets, saturated] = report_dense_bulk_surplus_to_nearest_interface( ...
    x, partsInCell, cntVec, c, highThrLocal, params)

movedCount = 0;
recvTargets = zeros(0,1);
saturated = false;

needOut = max(0, ceil(cntVec(c) - highThrLocal));
if needOut <= 0
    return;
end

idsAvail = partsInCell{c};
if isempty(idsAvail)
    return;
end

cand = geometric_interface_candidate_cells(c, cntVec, params);
if isempty(cand)
    saturated = true;
    return;
end

Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);
[jy, jx] = ind2sub([Ny, Nx], cand);
dist2 = zeros(size(cand));
for kk = 1:numel(cand)
    [dxCell, dyCell] = minimal_cell_step_general(ix, iy, jx(kk), jy(kk), params);
    dist2(kk) = dxCell.^2 + dyCell.^2;
end;
[dist2s, ord] = sort(dist2(:), 'ascend');
cand = cand(ord);

levels = unique(dist2s, 'stable');
for kk = 1:numel(levels)
    if needOut <= 0 || isempty(idsAvail)
        break;
    end

    thisLevel = cand(dist2s == levels(kk));
    if isempty(thisLevel)
        continue;
    end

    % Mélange isotrope des ex aequo à même distance géométrique.
    thisLevel = thisLevel(randperm(numel(thisLevel)));

    % Remplissage round-robin sur la coquille de distance minimale,
    % puis seulement ensuite sur la coquille suivante si nécessaire.
    progressed = true;
    while needOut > 0 && ~isempty(idsAvail) && progressed
        progressed = false;
        for jj = 1:numel(thisLevel)
            if needOut <= 0 || isempty(idsAvail)
                break;
            end
            nb = thisLevel(jj);
            idMove = idsAvail(1);
            idsAvail(1) = [];

            x(idMove,:) = sample_in_cell(nb, 1, params);
            partsInCell{nb} = [partsInCell{nb}; idMove];
            cntVec(nb) = cntVec(nb) + 1;

            recvTargets(end+1,1) = nb; %#ok<AGROW>
            movedCount = movedCount + 1;
            needOut = needOut - 1;
            progressed = true;
        end
    end
end

if movedCount > 0
    partsInCell{c} = idsAvail;
    cntVec(c) = numel(idsAvail);
    recvTargets = unique(recvTargets).';
else
    recvTargets = zeros(0,1);
end

saturated = cntVec(c) > highThrLocal;
end

function cand = geometric_interface_candidate_cells(cSource, cntVec, params)
Nx = params.Nx;
Ny = params.Ny;
activeMask = get_active_cell_mask(params, [Ny, Nx]);

candMask = false(Ny, Nx);
for iy = 1:Ny
    for ix = 1:Nx
        if ~activeMask(iy, ix)
            continue;
        end
        c = sub2ind([Ny, Nx], iy, ix);
        if c == cSource
            continue;
        end

        if is_geometric_interface_cell(c, cntVec, params)
            candMask(iy, ix) = true;
        end
    end
end
cand = find(candMask(:));
end

function tf = is_geometric_interface_cell(c, cntVec, params)
Nx = params.Nx;
Ny = params.Ny;
activeMask = get_active_cell_mask(params, [Ny, Nx]);
[iy, ix] = ind2sub([Ny, Nx], c);
if ~activeMask(iy, ix)
    tf = false;
    return;
end
if cntVec(c) <= 0
    tf = false;
    return;
end
tf = ~is_bulk_like_cell(c, cntVec, params);
end

function admiss = admissible_dense_targets_rangeR(c, cntVec, targetMap, params, rangeR)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

maxCount = (2*rangeR + 1)^2 - 1;
lst = zeros(maxCount,1);
m = 0;
for dy = -rangeR:rangeR
    for dx = -rangeR:rangeR
        if dx == 0 && dy == 0
            continue;
        end
        [ok, jy, jx] = shift_index_general(iy, ix, dy, dx, params);
        if ~ok
            continue;
        end
        nb = sub2ind([Ny, Nx], jy, jx);
        if cntVec(nb) < targetMap(jy, jx)
            m = m + 1;
            lst(m) = nb;
        end
    end
end

admiss = unique(lst(1:m).', 'stable');
end

function ordered = rank_dense_targets_unified(c, admiss, cntVec, targetMap, nHat, inwardWeight, params)
Nx = params.Nx;
Ny = params.Ny;
if isempty(admiss)
    ordered = zeros(0,1);
    return;
end

[iy, ix] = ind2sub([Ny, Nx], c);
nBins = 8;
bins = cell(nBins, 1);
startBin = randi(nBins);

for k = 1:numel(admiss)
    nb = admiss(k);
    [jy, jx] = ind2sub([Ny, Nx], nb);
    [dx, dy] = minimal_cell_step_general(ix, iy, jx, jy, params);
    ang = atan2(dy, dx);
    ib = mod(floor((ang + pi) / (2*pi/nBins)), nBins) + 1;

    deficit = max(0, targetMap(nb) - cntVec(nb));
    dist2 = dx*dx + dy*dy;
    cheb = max(abs(dx), abs(dy));

    nstep = hypot(dx, dy);
    inward = 0;
    if inwardWeight > 0 && nstep > 0
        inward = (dx*nHat(1) + dy*nHat(2)) / nstep;
    end

    bins{ib} = [bins{ib}; ...
        [nb, -deficit, -inwardWeight*inward, cntVec(nb), cheb, dist2, rand()]]; %#ok<AGROW>
end

for ib = 1:nBins
    if isempty(bins{ib})
        continue;
    end
    bins{ib} = sortrows(bins{ib}, [2 3 4 5 6 7]);
end

ordered = zeros(0,1);
keepGoing = true;
while keepGoing
    keepGoing = false;
    for s = 0:(nBins-1)
        ib = mod(startBin - 1 + s, nBins) + 1;
        if isempty(bins{ib})
            continue;
        end
        ordered(end+1,1) = bins{ib}(1,1); %#ok<AGROW>
        bins{ib}(1,:) = [];
        keepGoing = true;
    end
end
end

function [x, partsInCell, cntVec, movedCount, recvTargets] = spread_particles_from_dense_source_round_robin( ...
    x, partsInCell, cntVec, c, orderedTargets, targetMap, needOut, params)

movedCount = 0;
recvTargets = zeros(0,1);
if needOut <= 0 || isempty(orderedTargets)
    return;
end

idsAvail = partsInCell{c};
if isempty(idsAvail)
    return;
end

while needOut > 0 && ~isempty(idsAvail)
    progressed = false;

    for kk = 1:numel(orderedTargets)
        nb = orderedTargets(kk);
        targetLevel = targetMap(nb);
        deficit = max(0, floor(targetLevel - cntVec(nb)));
        if deficit <= 0
            continue;
        end

        remainingTargets = 0;
        for jj = kk:numel(orderedTargets)
            nbj = orderedTargets(jj);
            defj = max(0, floor(targetMap(nbj) - cntVec(nbj)));
            if defj > 0
                remainingTargets = remainingTargets + 1;
            end
        end

        quota = max(1, ceil(needOut / max(1, remainingTargets)));

        %nt = floor(mean([quota, deficit, needOut, numel(idsAvail)])); %%min
        %nt = min([quota, deficit, needOut, numel(idsAvail)]);
        nt = max(1, min([quota, deficit, needOut, numel(idsAvail)]));
        if nt <= 0
            continue;
        end

        idsMove = idsAvail(1:nt);
        idsAvail(1:nt) = [];

        x(idsMove,:) = sample_in_cell(nb, nt, params);
        partsInCell{nb} = [partsInCell{nb}; idsMove(:)];
        cntVec(nb) = cntVec(nb) + nt;
        recvTargets = [recvTargets; repmat(nb, nt>0, 1)]; %#ok<AGROW>

        movedCount = movedCount + nt;
        needOut = needOut - nt;
        progressed = true;

        if needOut <= 0 || isempty(idsAvail)
            break;
        end
    end

    if ~progressed
        break;
    end
end

if movedCount > 0
    partsInCell{c} = idsAvail;
    cntVec(c) = numel(idsAvail);
    recvTargets = unique(recvTargets);
else
    recvTargets = zeros(0,1);
end
end

function tf = has_nonempty_neighbor_excluding(c, cntVec, cExclude, params)
Nx = params.Nx;
Ny = params.Ny;

[iy, ix] = ind2sub([Ny, Nx], c);
[iy0, ix0] = ind2sub([Ny, Nx], cExclude);

% Si la cellule candidate est contiguë à la source (portée 1),
% on l'autorise même si la source est son seul voisin non vide.
if max(abs(ix - ix0), abs(iy - iy0)) == 1
    tf = true;
    return;
end

tf = false;

for dj = -1:1
    for di = -1:1
        if di == 0 && dj == 0
            continue;
        end

        jx = ix + di;
        jy = iy + dj;

        if jx < 1 || jx > Nx || jy < 1 || jy > Ny
            continue;
        end

        nb = sub2ind([Ny, Nx], jy, jx);
        if nb == cExclude
            continue;
        end

        if cntVec(nb) > 0
            tf = true;
            return;
        end
    end
end
end
function Ns = smooth_occupancy_field(N)
K = [1 2 1; 2 4 2; 1 2 1] / 16;
[Ny, Nx] = size(N);
Ns = zeros(size(N));
for iy = 1:Ny
    for ix = 1:Nx
        acc = 0;
        wsum = 0;
        for dy = -1:1
            for dx = -1:1
                jy = iy + dy;
                jx = ix + dx;
                if jy < 1 || jy > Ny || jx < 1 || jx > Nx
                    continue;
                end
                w = K(dy+2, dx+2);
                acc = acc + w * N(jy, jx);
                wsum = wsum + w;
            end
        end
        Ns(iy, ix) = acc / max(wsum, eps);
    end
end
end

function bulkMask = build_topological_bulk_mask(cntVec, params)
nCells = numel(cntVec);
bulkMask = false(nCells,1);
nz = find(cntVec > 0);
for kk = 1:numel(nz)
    c = nz(kk);
    bulkMask(c) = is_bulk_like_cell(c, cntVec, params);
end
end

function tf = is_bulk_like_cell(c, cntVec, params)
if is_wall_adjacent_cell(c, params)
    tf = false;
    return;
end

neighborsCore = params.neighborsCore;
neigh = neighborsCore{c};
if isempty(neigh)
    tf = false;
    return;
end

tf = all(cntVec(neigh) > 0);
end

function tf = has_only_diagonal_support(c, cntVec, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny,Nx], c);

orth = zeros(0,1);
for sh = [-1 1]
    [ok, jy, jx] = shift_index_general(iy, ix, sh, 0, params);
    if ok, orth(end+1,1) = sub2ind([Ny,Nx], jy, jx); end %#ok<AGROW>
    [ok, jy, jx] = shift_index_general(iy, ix, 0, sh, params);
    if ok, orth(end+1,1) = sub2ind([Ny,Nx], jy, jx); end %#ok<AGROW>
end
orth = unique(orth, 'stable');

diag = zeros(0,1);
for dy = [-1 1]
    for dx = [-1 1]
        [ok, jy, jx] = shift_index_general(iy, ix, dy, dx, params);
        if ok, diag(end+1,1) = sub2ind([Ny,Nx], jy, jx); end %#ok<AGROW>
    end
end
diag = unique(diag, 'stable');

tf = (sum(cntVec(orth) > 0) == 0) && (sum(cntVec(diag) > 0) >= 1);
end

function d = local_diff_x(Nmat, iy, ix, params)
Nx = params.Nx;
if is_periodic_pair_general(params, 'x')
    ixm = mod(ix - 2, Nx) + 1;
    ixp = mod(ix,     Nx) + 1;
    d = 0.5*(Nmat(iy, ixp) - Nmat(iy, ixm));
elseif ix == 1
    d = Nmat(iy, ix+1) - Nmat(iy, ix);
elseif ix == Nx
    d = Nmat(iy, ix) - Nmat(iy, ix-1);
else
    d = 0.5*(Nmat(iy, ix+1) - Nmat(iy, ix-1));
end
end

function d = local_diff_y(Nmat, iy, ix, params)
Ny = params.Ny;
if is_periodic_pair_general(params, 'y')
    iym = mod(iy - 2, Ny) + 1;
    iyp = mod(iy,     Ny) + 1;
    d = 0.5*(Nmat(iyp, ix) - Nmat(iym, ix));
elseif iy == 1
    d = Nmat(iy+1, ix) - Nmat(iy, ix);
elseif iy == Ny
    d = Nmat(iy, ix) - Nmat(iy-1, ix);
else
    d = 0.5*(Nmat(iy+1, ix) - Nmat(iy-1, ix));
end
end

function params = update_runtime_moving_piston(params, it)
Ly = params.Ly;
Ny = params.Ny;
Nx = params.Nx;

if isfield(params,'useMovingPiston') && params.useMovingPiston
    % Position du piston au pas courant
    if isfield(params,'pistonYCurrent') && ~isempty(params.pistonYCurrent)
        yPrev = params.pistonYCurrent;
    else
        yPrev = params.pistonY0;
    end

    yCur = params.pistonY0 + (it-1) * params.dt * params.pistonVy;
    if isfield(params,'pistonStopOnMin') && params.pistonStopOnMin
        if params.pistonVy < 0
            yCur = max(params.pistonYMin, yCur);
        else
            yCur = min(params.pistonYMin, yCur);
        end
    end

    vCur = params.pistonVy;
    if ((params.pistonVy < 0) && (yCur <= params.pistonYMin + eps)) || ...
       ((params.pistonVy > 0) && (yCur >= params.pistonYMin - eps))
        if isfield(params,'pistonStopOnMin') && params.pistonStopOnMin
            vCur = 0;
        end
    end

    params.pistonYPrev = yPrev;
    params.pistonYCurrent = yCur;
    params.pistonVyCurrent = vCur;

    fracMap = get_active_cell_fraction_map(params, [Ny, Nx]);
    fCut = 0.20;
    if isfield(params,'pistonFracCut') && ~isempty(params.pistonFracCut)
        fCut = params.pistonFracCut;
    end
    params.activeCellFracCurrent = fracMap;
    params.activeCellMaskCurrent = fracMap > fCut;

    nActiveCellsEff = max(sum(fracMap(:)), eps);
    gammaGeom = params.n / nActiveCellsEff;
    params.gammaGeomCurrent = gammaGeom;

    useConstTarget = true;
    if isfield(params,'pistonUseConstantGammaTarget') && ~isempty(params.pistonUseConstantGammaTarget)
        useConstTarget = logical(params.pistonUseConstantGammaTarget);
    end
    if useConstTarget && isfield(params,'gamma') && ~isempty(params.gamma) && isfinite(params.gamma)
        params.gammaTargetCurrent = params.gamma;
    else
        params.gammaTargetCurrent = gammaGeom;
    end
    % Compatibilite descendante : gammaCurrent garde ici le sens geometrique
    % "occupation requise si toute la masse reste sous le piston".
    params.gammaCurrent = gammaGeom;
else
    params.pistonYPrev = Ly;
    params.pistonYCurrent = Ly;
    params.pistonVyCurrent = 0;
    params.activeCellFracCurrent = ones(Ny, Nx);
    params.activeCellMaskCurrent = true(Ny, Nx);
    if isfield(params,'gamma') && ~isempty(params.gamma) && isfinite(params.gamma)
        params.gammaGeomCurrent = params.gamma;
        params.gammaTargetCurrent = params.gamma;
        params.gammaCurrent = params.gamma;
    else
        gammaDef = params.n / max(params.Nx*params.Ny,1);
        params.gammaGeomCurrent = gammaDef;
        params.gammaTargetCurrent = gammaDef;
        params.gammaCurrent = gammaDef;
    end
end
end

function [x, r0] = apply_affine_piston_compression(x, r0, params)
% Compression géométrique quasi-statique : le piston réduit la hauteur
% active et on recompose immédiatement les positions y du fluide dans
% cette hauteur, ce qui évite que seule la couche supérieure encaisse la
% compression.
if ~(isfield(params,'useMovingPiston') && params.useMovingPiston)
    return;
end
doAffine = true;
if isfield(params,'pistonAffineReposition') && ~isempty(params.pistonAffineReposition)
    doAffine = params.pistonAffineReposition;
end
if ~doAffine
    return;
end

yPrev = params.pistonYPrev;
yCur  = params.pistonYCurrent;
if ~(isfinite(yPrev) && isfinite(yCur)) || yPrev <= 0 || abs(yCur - yPrev) < eps
    return;
end

scale = yCur / yPrev;

% On compacte tout le fluide contenu sous l'ancienne position du piston.
% Toute particule numériquement au-dessus est rabattue juste sous yPrev.
yOld = min(max(x(:,2), 0), max(yPrev - eps, 0));
x(:,2) = min(max(yOld .* scale, 0), max(yCur - eps, 0));

if nargin >= 2 && ~isempty(r0)
    y0Old = min(max(r0(:,2), 0), max(yPrev - eps, 0));
    r0(:,2) = min(max(y0Old .* scale, 0), max(yCur - eps, 0));
end
end

function gammaEff = get_effective_gamma(params)
% Occupation cible effectivement utilisee par la redistribution.
if isfield(params,'gammaTargetCurrent') && ~isempty(params.gammaTargetCurrent) && isfinite(params.gammaTargetCurrent)
    gammaEff = params.gammaTargetCurrent;
elseif isfield(params,'gammaCurrent') && ~isempty(params.gammaCurrent) && isfinite(params.gammaCurrent)
    gammaEff = params.gammaCurrent;
elseif isfield(params,'gamma') && ~isempty(params.gamma) && isfinite(params.gamma)
    gammaEff = params.gamma;
else
    gammaEff = params.n / max(params.Nx*params.Ny,1);
end
end

function fracMap = get_active_cell_fraction_map(params, sz)
Ny = sz(1); Nx = sz(2);

% IMPORTANT :
% si le piston bouge, il faut RECALCULER la fraction active a chaque step.
% Le cache activeCellFracCurrent ne doit pas etre reutilise dans ce cas,
% sinon gammaCurrent reste fige.
if ~(isfield(params,'useMovingPiston') && params.useMovingPiston)
    if isfield(params,'activeCellFracCurrent') && ~isempty(params.activeCellFracCurrent)
        fracMap = params.activeCellFracCurrent;
        if ~isequal(size(fracMap), [Ny Nx])
            fracMap = ones(Ny, Nx);
        end
        return;
    end
end

fracMap = ones(Ny, Nx);
if isfield(params,'useMovingPiston') && params.useMovingPiston
    Ly = params.Ly;
    yTop = Ly;
    if isfield(params,'pistonYCurrent') && ~isempty(params.pistonYCurrent)
        yTop = params.pistonYCurrent;
    elseif isfield(params,'pistonY0') && ~isempty(params.pistonY0)
        yTop = params.pistonY0;
    end
    margin = 0;
    if isfield(params,'pistonActiveMargin') && ~isempty(params.pistonActiveMargin)
        margin = params.pistonActiveMargin;
    end
    yTopEff = min(max(yTop - margin, 0), Ly);
    aY = Ly / Ny;

    fracRows = zeros(Ny,1);
    for iy = 1:Ny
        y0 = (iy-1) * aY;
        y1 = iy * aY;
        fracRows(iy) = max(0, min(yTopEff, y1) - y0) / aY;
    end
    fracMap = repmat(fracRows, 1, Nx);
end
end

function activeMask = get_active_cell_mask(params, sz)
Ny = sz(1); Nx = sz(2);

% IMPORTANT :
% Pour un piston mobile, le masque actif doit etre RECALCULE a chaque step
% a partir de la position courante du piston. Si on reutilise
% activeCellMaskCurrent en cache, gammaCurrent reste fige.
if ~(isfield(params,'useMovingPiston') && params.useMovingPiston)
    if isfield(params,'activeCellMaskCurrent') && ~isempty(params.activeCellMaskCurrent)
        activeMask = logical(params.activeCellMaskCurrent);
        if ~isequal(size(activeMask), [Ny Nx])
            activeMask = true(Ny, Nx);
        end
        return;
    end
end

activeMask = true(Ny, Nx);
if isfield(params,'useMovingPiston') && params.useMovingPiston
    Ly = params.Ly;
    yTop = Ly;
    if isfield(params,'pistonYCurrent') && ~isempty(params.pistonYCurrent)
        yTop = params.pistonYCurrent;
    elseif isfield(params,'pistonY0') && ~isempty(params.pistonY0)
        yTop = params.pistonY0;
    end
    aY = Ly / Ny;
    yc = ((0:Ny-1) + 0.5)' * aY;
    margin = 0;
    if isfield(params,'pistonActiveMargin') && ~isempty(params.pistonActiveMargin)
        margin = params.pistonActiveMargin;
    end
    rowActive = yc <= (yTop - margin);
    activeMask = repmat(rowActive, 1, Nx);
end
end

function [cntVec, partsInCell] = rebuild_parts_in_cell(x, params)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
aX = Lx/Nx; aY = Ly/Ny;
ix = floor(x(:,1)/aX) + 1;
iy = floor(x(:,2)/aY) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
ic = sub2ind([Ny,Nx], iy, ix);
cntVec = accumarray(ic, 1, [Nx*Ny,1], @sum, 0);
partsInCell = accumarray(ic, (1:size(x,1))', [Nx*Ny,1], @(z){z}, {[]});
end

function [Px, Py] = cell_total_momentum(partsInCell, v, params)
Nc = params.Nc;
Px = zeros(Nc,1); Py = zeros(Nc,1);
for c = 1:Nc
    ids = partsInCell{c};
    if ~isempty(ids)
        pp = sum(v(ids,:),1);
        Px(c) = pp(1); Py(c) = pp(2);
    end
end
end


function [fMap, gammaLocMap, lowThrLocMap, highThrLocMap] = build_local_fluid_fraction_threshold_maps(Nsm, params)
gamma = get_effective_gamma(params);
coef = params.coef;
useLocalFluidFractionThresholds = params.useLocalFluidFractionThresholds;
fluidFracGain = params.fluidFracGain;
lowThrFloor = params.lowThrFloor;
highThrFloor = params.highThrFloor;
occAbsFloor = 1;
if isfield(params,'occAbsFloor') && ~isempty(params.occAbsFloor)
    occAbsFloor = params.occAbsFloor;
end

activeFrac = get_active_cell_fraction_map(params, size(Nsm));
fCut = 0.20;
if isfield(params,'pistonFracCut') && ~isempty(params.pistonFracCut)
    fCut = params.pistonFracCut;
end

% Cible locale pilotée explicitement par la fraction de volume active
% (par exemple sous un piston mobile, mais sans codage implicite de cas).
if isfield(params,'redistributionUseActiveVolumeTargets') && params.redistributionUseActiveVolumeTargets
    fMap = activeFrac;
    fMap(fMap <= fCut) = 0;

    gammaLocMap = gamma * fMap;
    dNloc = max(occAbsFloor, round(coef * max(gammaLocMap, 1)));
    lowThrLocMap = max(0, gammaLocMap - dNloc);
    highThrLocMap = gammaLocMap + dNloc;
    lowThrLocMap(fMap <= 0) = 0;
    highThrLocMap(fMap <= 0) = 0;
    return;
end

activeMask = activeFrac > fCut;

% Estimation minimale d'une fraction locale de fluide a partir de Nsm/gamma,
% puis conversion en cible locale gammaLoc et seuils locaux d'occupation.
fMap = zeros(size(Nsm));
if ~useLocalFluidFractionThresholds
    fMap(activeMask) = 1;
else
    fRaw = zeros(size(Nsm));
    fRaw(activeMask) = Nsm(activeMask) / max(gamma, 1e-12);
    fMap(activeMask) = min(max(fluidFracGain * fRaw(activeMask), 0), 1);
end

gammaLocMap = zeros(size(Nsm));
gammaLocMap(activeMask) = gamma * fMap(activeMask);

lowFrac  = 1 / (1 + coef);
highFrac = 1 + coef;

lowThrLocMap = zeros(size(Nsm));
highThrLocMap = zeros(size(Nsm));
lowThrLocMap(activeMask)  = max(lowFrac  * gammaLocMap(activeMask), lowThrFloor);
highThrLocMap(activeMask) = max(highFrac * gammaLocMap(activeMask), highThrFloor);
end

function [highThr, lowThrBulk, lowThrInterface] = occupancy_thresholds(params)
gamma = get_effective_gamma(params);
coef = params.coef;
highMode = params.highMode;
lowMode = params.lowMode;
lowThrBulkOverride = params.lowThrBulkOverride;
lowThrInterfaceOverride = params.lowThrInterfaceOverride;
occAbsFloor = 1;
if isfield(params,'occAbsFloor') && ~isempty(params.occAbsFloor)
    occAbsFloor = params.occAbsFloor;
end
% Seuil haut commun ; seuil bas différencié bulk / interface.
% En mode de cibles actives locales, on emploie une bande mixte avec
% plancher absolu autour de la cible gamma locale.
if isfield(params,'redistributionUseActiveVolumeTargets') && params.redistributionUseActiveVolumeTargets
    dN = max(occAbsFloor, round(coef * gamma));
    highThr = gamma + dN;
    lowThrBulk = max(0, gamma - dN);
else
    if strcmpi(highMode, '1+coef')
        highThr = gamma * (1 + coef);
    else
        highThr = gamma + gamma*coef;
    end

    if strcmpi(lowMode, '1+coef')
        lowThrBulk = gamma / (1 + coef);
    else
        lowThrBulk = gamma - gamma*coef;
    end
end

% Par défaut, l'interface est un peu plus permissive que le bulk.
lowThrInterface = max(0, lowThrBulk - 1);

if isfinite(lowThrBulkOverride)
    lowThrBulk = lowThrBulkOverride;
end
if isfinite(lowThrInterfaceOverride)
    lowThrInterface = lowThrInterfaceOverride;
end

lowThrInterface = min(lowThrInterface, lowThrBulk);
end

function xx = sample_in_cell(c, n, params)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny,Nx], c);
aX = Lx/Nx; aY = Ly/Ny;
x0 = (ix-1)*aX; y0 = (iy-1)*aY;
y1 = y0 + aY;
if isfield(params,'useMovingPiston') && params.useMovingPiston
    yTop = min(max(params.pistonYCurrent, 0), Ly);
    y1 = min(y1, max(yTop - eps, y0 + eps));
end
xx = [x0 + aX*rand(n,1), y0 + (y1-y0)*rand(n,1)];
end


function Ncell = cell_occupancy(x, params)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
aX = Lx/Nx; aY = Ly/Ny;
ix = floor(x(:,1)/aX) + 1;
iy = floor(x(:,2)/aY) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
ic = sub2ind([Ny,Nx], iy, ix);
Nvec = accumarray(ic, 1, [Nx*Ny,1], @sum, 0);
Ncell = reshape(Nvec, [Ny,Nx]);
end

function diag = default_redist_diag(Ncell, params)
gamma = get_effective_gamma(params);
N = Ncell(:);
activeMask = get_active_cell_mask(params, size(Ncell));
activeMaskVec = activeMask(:);
diag = zeros(1,10);

if any(activeMaskVec)
    Nact = double(N(activeMaskVec));
    diag(1) = std(Nact)/max(gamma, eps);
    diag(2) = mean(abs(Nact-gamma) > max(1, round(0.2*gamma)));
    diag(10) = mean(N(activeMaskVec) == 0);
else
    diag(1) = 0;
    diag(2) = 0;
    diag(10) = 0;
end
end
function angMis = local_normal_misalignment(c, nx, ny, cntVec, params)
Nx = params.Nx;
Ny = params.Ny;
% Ecart angulaire moyen entre la normale locale et celles des voisines
% non vides. Sert de proxy de courbure / angularité locale.

[iy, ix] = ind2sub([Ny, Nx], c);

n0 = [nx(iy,ix), ny(iy,ix)];
nn0 = hypot(n0(1), n0(2));
if nn0 < 1e-12
    angMis = 0;
    return;
end
n0 = n0 / nn0;

angVals = zeros(8,1);
m = 0;

for dy = -1:1
    for dx = -1:1
        if dx == 0 && dy == 0
            continue;
        end

        jy = iy + dy;
        jx = ix + dx;
        if jy < 1 || jy > Ny || jx < 1 || jx > Nx
            continue;
        end

        cn = sub2ind([Ny, Nx], jy, jx);
        if cntVec(cn) <= 0
            continue;
        end

        nv = [nx(jy,jx), ny(jy,jx)];
        nnv = hypot(nv(1), nv(2));
        if nnv < 1e-12
            continue;
        end
        nv = nv / nnv;

        cs = max(-1, min(1, n0(1)*nv(1) + n0(2)*nv(2)));
        m = m + 1;
        angVals(m) = acos(cs);
    end
end

if m == 0
    angMis = 0;
else
    angMis = mean(angVals(1:m));
end
end




function flags = resolve_redistribution_flags(params)
flags = struct();
flags.enableSurfaceTopology = isfield(params,'redistributionEnableSurfaceTopology') && ...
                              logical(params.redistributionEnableSurfaceTopology);
flags.useActiveVolumeTargets = isfield(params,'redistributionUseActiveVolumeTargets') && ...
                               logical(params.redistributionUseActiveVolumeTargets);

if isfield(params,'redistributionBulkOverflowPolicy') && ~isempty(params.redistributionBulkOverflowPolicy)
    flags.bulkOverflowPolicy = lower(char(params.redistributionBulkOverflowPolicy));
else
    flags.bulkOverflowPolicy = 'none';
end

flags.wallWettingEnabled = isfield(params,'redistributionWallWettingEnabled') && ...
                           logical(params.redistributionWallWettingEnabled);
if isfield(params,'redistributionPreserveWettingLayer') && ~isempty(params.redistributionPreserveWettingLayer)
    flags.preserveWettingLayer = logical(params.redistributionPreserveWettingLayer);
else
    flags.preserveWettingLayer = flags.wallWettingEnabled;
end
flags.wettingSides = resolve_redistribution_wetting_sides(params);
if flags.wettingSides.count == 0
    flags.wallWettingEnabled = false;
    flags.preserveWettingLayer = false;
end
end

function sides = resolve_redistribution_wetting_sides(params)
sides = struct();
sides.left   = strcmpi(get_boundary_mode_general(params,'left'),   'wall_wetting');
sides.right  = strcmpi(get_boundary_mode_general(params,'right'),  'wall_wetting');
sides.bottom = strcmpi(get_boundary_mode_general(params,'bottom'), 'wall_wetting');
sides.top    = strcmpi(get_boundary_mode_general(params,'top'),    'wall_wetting');

if isfield(params,'useBottomWallWetting') && ~isempty(params.useBottomWallWetting) && params.useBottomWallWetting
    sides.bottom = true;
end
sides.count = double(sides.left) + double(sides.right) + double(sides.bottom) + double(sides.top);
end

function tf = is_physical_boundary_side(params, side)
tf = ~strcmpi(get_boundary_mode_general(params, side), 'periodic');
end

function tf = is_cell_adjacent_to_wetting_boundary(c, params, flags)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);
tf = (ix == 1  && flags.wettingSides.left)  || ...
     (ix == Nx && flags.wettingSides.right) || ...
     (iy == 1  && flags.wettingSides.bottom) || ...
     (iy == Ny && flags.wettingSides.top);
end

function wettingTarget = resolve_wall_wetting_target(lowThrBulk, params)
if isfield(params,'redistributionWettingTargetOverride') && isfinite(params.redistributionWettingTargetOverride)
    wettingTarget = ceil(params.redistributionWettingTargetOverride);
elseif isfield(params,'bottomWettingTargetOverride') && isfinite(params.bottomWettingTargetOverride)
    wettingTarget = ceil(params.bottomWettingTargetOverride);
else
    wettingTarget = ceil(lowThrBulk);
end
wettingTarget = max(1, wettingTarget);
end

function [x, partsInCell, cntVec, movedCount, recvTargets, srcTouched] = apply_wall_wetting_layers_general( ...
    x, partsInCell, cntVec, wettingTarget, params, flags)

movedCount = 0;
recvTargets = zeros(0,1);
srcTouched = zeros(0,1);

sideOrder = {'bottom','top','left','right'};
for ks = 1:numel(sideOrder)
    side = sideOrder{ks};
    if ~flags.wettingSides.(side)
        continue;
    end
    [x, partsInCell, cntVec, movedNow, recvNow, srcNow] = apply_wall_wetting_side_layer( ...
        x, partsInCell, cntVec, wettingTarget, side, params);
    movedCount = movedCount + movedNow;
    if ~isempty(recvNow)
        recvTargets = unique([recvTargets(:); recvNow(:)]);
    end
    if ~isempty(srcNow)
        srcTouched = unique([srcTouched(:); srcNow(:)]);
    end
end
recvTargets = recvTargets(:).';
srcTouched = srcTouched(:).';
end

function [x, partsInCell, cntVec, movedCount, recvTargets, srcTouched] = apply_wall_wetting_side_layer( ...
    x, partsInCell, cntVec, wettingTarget, side, params)

Nx = params.Nx;
Ny = params.Ny;
movedCount = 0;
recvTargets = zeros(0,1);
srcTouched = zeros(0,1);

sideCells = boundary_side_cells(side, params);
if isempty(sideCells) || wettingTarget <= 0
    return;
end

sourceList = sideCells(cntVec(sideCells) > wettingTarget);
if isempty(sourceList)
    return;
end
sourceList = sourceList(randperm(numel(sourceList)));

targetMap = zeros(Ny, Nx);
targetMap(sideCells) = wettingTarget;

for kk = 1:numel(sourceList)
    c = sourceList(kk);
    if cntVec(c) <= wettingTarget
        continue;
    end

    needOut = cntVec(c) - wettingTarget;
    admiss = admissible_wall_wetting_targets_side(c, cntVec, wettingTarget, side, params);
    if isempty(admiss)
        continue;
    end

    ordered = rank_wall_wetting_targets_side(c, admiss, side, params);
    [x, partsInCell, cntVec, movedNow, recvNow] = spread_particles_from_dense_source_round_robin( ...
        x, partsInCell, cntVec, c, ordered, targetMap, needOut, params);

    movedCount = movedCount + movedNow;
    if movedNow > 0
        srcTouched = unique([srcTouched(:); c]);
    end
    if ~isempty(recvNow)
        recvTargets = unique([recvTargets(:); recvNow(:)]);
    end
end
recvTargets = recvTargets(:).';
srcTouched = srcTouched(:).';
end

function cells = boundary_side_cells(side, params)
Nx = params.Nx;
Ny = params.Ny;
switch lower(side)
    case 'bottom'
        cells = arrayfun(@(jx) sub2ind([Ny, Nx], 1,  jx), 1:Nx);
    case 'top'
        cells = arrayfun(@(jx) sub2ind([Ny, Nx], Ny, jx), 1:Nx);
    case 'left'
        cells = arrayfun(@(jy) sub2ind([Ny, Nx], jy, 1), 1:Ny);
    case 'right'
        cells = arrayfun(@(jy) sub2ind([Ny, Nx], jy, Nx), 1:Ny);
    otherwise
        cells = zeros(1,0);
end
cells = cells(:).';
end

function admiss = admissible_wall_wetting_targets_side(c, cntVec, wettingTarget, side, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

if isfield(params,'redistributionWettingMaxRange') && isfinite(params.redistributionWettingMaxRange)
    maxRange = max(1, floor(params.redistributionWettingMaxRange));
elseif isfield(params,'bottomWettingMaxRange') && isfinite(params.bottomWettingMaxRange)
    maxRange = max(1, floor(params.bottomWettingMaxRange));
else
    maxRange = max(Nx, Ny);
end

lst = zeros(1, 2*maxRange);
m = 0;
for dj = 1:maxRange
    switch lower(side)
        case {'bottom','top'}
            jxCand = [ix - dj, ix + dj];
            jyCand = [iy, iy];
        case {'left','right'}
            jxCand = [ix, ix];
            jyCand = [iy - dj, iy + dj];
        otherwise
            admiss = zeros(0,1);
            return;
    end

    for k = 1:2
        jy = jyCand(k);
        jx = jxCand(k);
        if jy < 1 || jy > Ny || jx < 1 || jx > Nx
            continue;
        end
        nb = sub2ind([Ny, Nx], jy, jx);
        if cntVec(nb) < wettingTarget
            m = m + 1;
            lst(m) = nb;
        end
    end
end
admiss = unique(lst(1:m), 'stable');
end

function ordered = rank_wall_wetting_targets_side(c, admiss, side, params)
Nx = params.Nx;
Ny = params.Ny;
[iy, ix] = ind2sub([Ny, Nx], c);

if isempty(admiss)
    ordered = zeros(0,1);
    return;
end

dist = zeros(numel(admiss),1);
for k = 1:numel(admiss)
    [jy, jx] = ind2sub([Ny, Nx], admiss(k));
    switch lower(side)
        case {'bottom','top'}
            dist(k) = abs(jx - ix);
        case {'left','right'}
            dist(k) = abs(jy - iy);
        otherwise
            dist(k) = inf;
    end
end

dUnique = unique(dist(:)).';
ordered = zeros(1,0);
for d = dUnique
    mask = (dist == d);
    cand = admiss(mask);
    if numel(cand) > 1
        cand = cand(randperm(numel(cand)));
    end
    ordered = [ordered, cand(:).']; %#ok<AGROW>
end
end

function [dxStep, dyStep] = minimal_cell_step_general(ix, iy, jx, jy, params)
Nx = params.Nx;
Ny = params.Ny;

dxStep = jx - ix;
if is_periodic_pair_general(params, 'x')
    if dxStep >  Nx/2, dxStep = dxStep - Nx; end
    if dxStep < -Nx/2, dxStep = dxStep + Nx; end
end

dyStep = jy - iy;
if is_periodic_pair_general(params, 'y')
    if dyStep >  Ny/2, dyStep = dyStep - Ny; end
    if dyStep < -Ny/2, dyStep = dyStep + Ny; end
end
end

function [ok, jy, jx] = shift_index_general(iy, ix, dy, dx, params)
Ny = params.Ny;
Nx = params.Nx;

jy = iy + dy;
if is_periodic_pair_general(params, 'y')
    jy = mod(jy - 1, Ny) + 1;
elseif jy < 1 || jy > Ny
    ok = false;
    jx = NaN;
    return;
end

jx = ix + dx;
if is_periodic_pair_general(params, 'x')
    jx = mod(jx - 1, Nx) + 1;
elseif jx < 1 || jx > Nx
    ok = false;
    return;
end

ok = true;
end




function [Ux, Uy, Cnt, Om] = velocity_field_and_vorticity(x, v, params)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
aX = Lx/Nx; aY = Ly/Ny;
ix = floor(x(:,1)/aX) + 1;
iy = floor(x(:,2)/aY) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
ic = sub2ind([Ny,Nx], iy, ix);

Cvec  = accumarray(ic, 1,      [Nx*Ny,1], @sum, 0);
Uxvec = accumarray(ic, v(:,1), [Nx*Ny,1], @sum, 0);
Uyvec = accumarray(ic, v(:,2), [Nx*Ny,1], @sum, 0);

Uxv = zeros(Nx*Ny,1); Uyv = zeros(Nx*Ny,1);
nz = Cvec > 0;
Uxv(nz) = Uxvec(nz)./Cvec(nz);
Uyv(nz) = Uyvec(nz)./Cvec(nz);

Ux = reshape(Uxv, [Ny,Nx]);
Uy = reshape(Uyv, [Ny,Nx]);
Cnt = reshape(Cvec, [Ny,Nx]);

[dUy_dx, ~] = gradient(Uy, aX, aY);
[~, dUx_dy] = gradient(Ux, aX, aY);
Om = dUy_dx - dUx_dy;
end

