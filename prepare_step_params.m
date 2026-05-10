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

if ~isfield(params, 'useLiquidClosure') || isempty(params.useLiquidClosure), params.useLiquidClosure = false; end
if ~isfield(params, 'useOptimalBetaRepair') || isempty(params.useOptimalBetaRepair), params.useOptimalBetaRepair = true; end
if ~isfield(params, 'betaRepair') || isempty(params.betaRepair), params.betaRepair = 0.0; end
if ~isfield(params, 'betaEOS') || isempty(params.betaEOS), params.betaEOS = 1.0; end
if ~isfield(params, 'Kvirial') || isempty(params.Kvirial), params.Kvirial = 0.0; end
if ~isfield(params, 'liquidClosureDumpStages') || isempty(params.liquidClosureDumpStages), params.liquidClosureDumpStages = false; end

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

if ~isfield(params,'accumulateMeanFields') || isempty(params.accumulateMeanFields)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(params,'useLiquidClosure') && ~isempty(params.useLiquidClosure)
        params.accumulateMeanFields = params.useLiquidClosure;
    else
        params.accumulateMeanFields = false;
    end
end

if ~isfield(params,'meanFieldsStartStep') || isempty(params.meanFieldsStartStep)
    params.meanFieldsStartStep = 1;
end
if ~isfield(params,'meanFieldsStride') || isempty(params.meanFieldsStride)
    params.meanFieldsStride = 1;
end

params.stepParamsPrepared = true;
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
            mode = lower(char(params.xBoundary));
        else
            mode = 'specular';
        end
    case 'bottom'
        if isfield(params,'yWallMode') && ~isempty(params.yWallMode)
            mode = lower(char(params.yWallMode));
        else
            mode = 'specular';
        end
    case 'top'
        if isfield(params,'useMovingPiston') && params.useMovingPiston
            mode = 'piston';
        elseif isfield(params,'yWallMode') && ~isempty(params.yWallMode)
            mode = lower(char(params.yWallMode));
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

function neighbors = precompute_neighbors_general(params)
Nx = params.Nx;
Ny = params.Ny;
periodicX = is_periodic_pair_general(params,'x');
periodicY = is_periodic_pair_general(params,'y');
neighbors = cell(Nx*Ny,1);
for c = 1:(Nx*Ny)
    [iy, ix] = ind2sub([Ny,Nx], c);
    list = zeros(8,1);
    k = 0;
    for dy = -1:1
        for dx = -1:1
            if dx == 0 && dy == 0, continue; end
            jy = iy + dy;
            if periodicY
                jy = mod(jy-1, Ny) + 1;
            else
                if jy < 1 || jy > Ny, continue; end
            end
            jx = ix + dx;
            if periodicX
                jx = mod(jx-1, Nx) + 1;
            else
                if jx < 1 || jx > Nx, continue; end
            end
            k = k + 1;
            list(k) = sub2ind([Ny,Nx], jy, jx);
        end
    end
    nb = unique(list(1:k).', 'stable');
    if ~isempty(nb)
        nb = nb(randperm(numel(nb)));
    end
    neighbors{c} = nb;
end
end
