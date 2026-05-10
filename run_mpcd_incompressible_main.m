function runOut = run_mpcd_incompressible_main(cfg)
%RUN_MPCD_INCOMPRESSIBLE_MAIN_BOUNDARY_OPS
% Runner principal compatible avec le noyau refactorise boundary_ops.
% Il conserve l'entree par caseType pour definir des presets coherents,
% mais construit maintenant un jeu de params explicite pilote par :
%   - boundary_left / boundary_right / boundary_bottom / boundary_top
%   - useMovingPiston + parametres piston
%   - drapeaux explicites de redistribution et de diagnostics
%
% Usage:
%   runOut = run_mpcd_incompressible_main_boundary_ops(cfg)

if nargin < 1 || isempty(cfg), cfg = struct(); end
if ~isstruct(cfg) || ~isscalar(cfg), error('cfg doit etre un struct scalaire.'); end

% Par defaut, on conserve un mode de compatibilite stricte avec l'ancien
% runner : on construit les anciens champs historiques, et on n'injecte les
% nouveaux champs explicites (boundary_*, piston, redistribution*, etc.)
% que si l'utilisateur desactive ce mode. Cela evite qu'un vieux cfg qui
% traine des champs explicites modifie la dynamique d'un cas legacy.
strictLegacyCaseSemantics = getv(cfg,'strictLegacyCaseSemantics',true);

caseType = lower(string(getv(cfg,'caseType','diffusion')));
rngSeed  = getv(cfg,'rngSeed',1);
rng(rngSeed);

% Gestion noyau / sorties
coreFcnName      = getv(cfg,'coreFcnName','');
if strlength(string(coreFcnName)) == 0
    coreFcnName = 'mpcd_incompressible_simulator_pistondiag';
end
saveRunResults   = getv(cfg,'saveRunResults',true);
resultsRoot      = getv(cfg,'resultsRoot',fullfile(pwd,'mpcd_results_densitywave'));
runLabel         = getv(cfg,'runLabel','');
saveInitialState = getv(cfg,'saveInitialState',true);
saveFinalState   = getv(cfg,'saveFinalState',false);
saveSimOutStruct = getv(cfg,'saveSimOutStruct',false);
appendSummaryCSV = getv(cfg,'appendSummaryCSV',true);
runNotes         = getv(cfg,'runNotes','');

% Domaine / maillage
Lx = getv(cfg,'Lx',10.0);
Ly = getv(cfg,'Ly',10.0);
Nx = getv(cfg,'Nx',50);
Ny = getv(cfg,'Ny',50);
Nc = Nx*Ny;

% MPCD commun
gamma = getv(cfg,'gamma',10);
nTotalParticles = getv(cfg,'nTotalParticles',round(gamma*Nc)); %#ok<NASGU>
dt = getv(cfg,'dt',5e-3);
nSteps = getv(cfg,'nSteps',750);
alphaDeg = getv(cfg,'alphaDeg',170);
alpha = deg2rad(alphaDeg);
kBT = getv(cfg,'kBT',1.0);
g = getv(cfg,'g',getv(cfg,'gy',getv(cfg,'gravityY',0.0)));
bodyForceX = getv(cfg,'bodyForceX',0.0);
useThermostat = getv(cfg,'useThermostat',true);
keepMeanFlow = getv(cfg,'keepMeanFlow',false);

% BL
xBoundary = getv(cfg,'xBoundary','specular');
yWallMode = getv(cfg,'yWallMode','specular');
Utop = getv(cfg,'Utop',0.0);
Ubottom = getv(cfg,'Ubottom',0.0);
wallSigma = getv(cfg,'wallSigma',sqrt(max(kBT,0)));

% CL explicites / piston
boundary_left_cfg = getv(cfg,'boundary_left','');
boundary_right_cfg = getv(cfg,'boundary_right','');
boundary_bottom_cfg = getv(cfg,'boundary_bottom','');
boundary_top_cfg = getv(cfg,'boundary_top','');
useMovingPiston = getv(cfg,'useMovingPiston',false);
pistonY0 = getv(cfg,'pistonY0',Ly);
pistonVy = getv(cfg,'pistonVy',0.0);
pistonYMin = getv(cfg,'pistonYMin',0.5*Ly);
pistonStopOnMin = getv(cfg,'pistonStopOnMin',true);
pistonFracCut = getv(cfg,'pistonFracCut',0.20);
pistonAffineReposition = getv(cfg,'pistonAffineReposition',true);
pistonActiveMargin = getv(cfg,'pistonActiveMargin',0.0);
occAbsFloor = getv(cfg,'occAbsFloor',1);

% Redistribution
useIncompressibleRedistribution = getv(cfg,'useIncompressibleRedistribution',true);
coef = getv(cfg,'coef',0.10);
highMode = getv(cfg,'highMode','coef');
lowMode = getv(cfg,'lowMode','coef');
redistribAfterCollision = getv(cfg,'redistribAfterCollision',true);
maxRedistribPasses = getv(cfg,'maxRedistribPasses',2);
enableMomentumCorrectionPostRedistribution = getv(cfg,'enableMomentumCorrectionPostRedistribution',true);
useLocalFluidFractionThresholds = getv(cfg,'useLocalFluidFractionThresholds',true);
fluidFracGain = getv(cfg,'fluidFracGain',0.8);
lowThrFloor = getv(cfg,'lowThrFloor',20);
highThrFloor = getv(cfg,'highThrFloor',25);
lowThrBulkOverride = getv(cfg,'lowThrBulkOverride',NaN);
lowThrInterfaceOverride = getv(cfg,'lowThrInterfaceOverride',10);
useBottomWallWetting = getv(cfg,'useBottomWallWetting',false);
bottomWettingTargetOverride = getv(cfg,'bottomWettingTargetOverride',NaN);
bottomWettingMaxRange = getv(cfg,'bottomWettingMaxRange',inf);
noSurfaceCase = getv(cfg,'noSurfaceCase',[]);
useLocalTargetsForPiston = getv(cfg,'useLocalTargetsForPiston',[]);
redistributionEnableSurfaceTopology = getv(cfg,'redistributionEnableSurfaceTopology',[]);
redistributionUseActiveVolumeTargets = getv(cfg,'redistributionUseActiveVolumeTargets',[]);
redistributionBulkOverflowPolicy = getv(cfg,'redistributionBulkOverflowPolicy','');
redistributionWallWettingEnabled = getv(cfg,'redistributionWallWettingEnabled',[]);
redistributionPreserveWettingLayer = getv(cfg,'redistributionPreserveWettingLayer',[]);
redistributionWettingTargetOverride = getv(cfg,'redistributionWettingTargetOverride',NaN);
redistributionWettingMaxRange = getv(cfg,'redistributionWettingMaxRange',NaN);

% Reorientation interfaciale
useInterfaceVelocityReorientation = getv(cfg,'useInterfaceVelocityReorientation',false);
interfaceReorientBeta = getv(cfg,'interfaceReorientBeta',0.36);
interfaceReorientFrac = getv(cfg,'interfaceReorientFrac',0.35);
useCurvatureWeightedReorientation = getv(cfg,'useCurvatureWeightedReorientation',true);
curvatureReorientBetaGain = getv(cfg,'curvatureReorientBetaGain',5);
curvatureReorientFracGain = getv(cfg,'curvatureReorientFracGain',5);
interfaceReorientBetaMax = getv(cfg,'interfaceReorientBetaMax',1);
interfaceReorientFracMax = getv(cfg,'interfaceReorientFracMax',1);
curvatureAngleMax = getv(cfg,'curvatureAngleMax',pi/6);

% Fermeture liquide
useLiquidClosure = getv(cfg,'useLiquidClosure',true);
stepFcnName = getv(cfg,'stepFcnName','mpcd_incompressible_step_liquidclosure');
useOptimalBetaRepair = getv(cfg,'useOptimalBetaRepair',true);
betaRepair = getv(cfg,'betaRepair',0.0);
betaEOS = getv(cfg,'betaEOS',1.0);
Kvirial = getv(cfg,'Kvirial',0.0);
liquidClosureDumpStages = getv(cfg,'liquidClosureDumpStages',false);

% Diagnostics / runtime
makePlots = getv(cfg,'makePlots',true);
realtimeStride = getv(cfg,'realtimeStride',50);
nplot = getv(cfg,'nplot',realtimeStride);
dumpStride = getv(cfg,'dumpStride',8);
quiverStrideX = getv(cfg,'quiverStrideX',8);
quiverStrideY = getv(cfg,'quiverStrideY',8);
quiverScale = getv(cfg,'quiverScale',2.0);
vortClipSigma = getv(cfg,'vortClipSigma',3.0);
diagStride = getv(cfg,'diagStride',20);
viscoStride = getv(cfg,'viscoStride',20);
movAvgWindow = getv(cfg,'movAvgWindow',200);
fitWindowFrac = getv(cfg,'fitWindowFrac',0.5);
steadyFracStart = getv(cfg,'steadyFracStart',0.6);
fitMarginFrac = getv(cfg,'fitMarginFrac',0.10);
fitSmoothWindowY = getv(cfg,'fitSmoothWindowY',7);
nRadialBins = getv(cfg,'nRadialBins',50);

% Diagnostics explicites
computeShapeDiagnostics = getv(cfg,'computeShapeDiagnostics',[]);
computeTracerDiagnostics = getv(cfg,'computeTracerDiagnostics',[]);
computeDensityWaveDiagnostics = getv(cfg,'computeDensityWaveDiagnostics',[]);
densityWaveDiagMode = getv(cfg,'densityWaveDiagMode',[]);
computePressureSummary = getv(cfg,'computePressureSummary',[]);
computePoiseuilleDiagnostics = getv(cfg,'computePoiseuilleDiagnostics',[]);
computeViscosityDiagnostics = getv(cfg,'computeViscosityDiagnostics',[]);
printSummary = getv(cfg,'printSummary',true);

% Auto-stop
enableAutoStop = getv(cfg,'enableAutoStop',true);
convCheckEvery = getv(cfg,'convCheckEvery',max([10, diagStride, viscoStride]));
convPatienceChecks = getv(cfg,'convPatienceChecks',2);
ellipseMinSteps = getv(cfg,'ellipseMinSteps',max(200, round(0.12*nSteps)));
ellipseConvWindow = getv(cfg,'ellipseConvWindow',max(40, round(0.05*nSteps)));
ellipseAtol = getv(cfg,'ellipseAtol',5e-2);
ellipseRelRTol = getv(cfg,'ellipseRelRTol',30.0);
ellipseDeltaATol = getv(cfg,'ellipseDeltaATol',8e-3);
ellipseRangeATol = getv(cfg,'ellipseRangeATol',15e-3);
diffusionMinSteps = getv(cfg,'diffusionMinSteps',max(200, round(0.10*nSteps)));
diffusionConvWindow = getv(cfg,'diffusionConvWindow',max(80, round(0.08*nSteps)));
diffusionRelDTol = getv(cfg,'diffusionRelDTol',3e-2);
diffusionR2linMin = getv(cfg,'diffusionR2linMin',0.995);
diffusionR2maxFrac = getv(cfg,'diffusionR2maxFrac',0.08);
poiseuilleMinSteps = getv(cfg,'poiseuilleMinSteps',max(500, round(0.20*nSteps)));
poiseuilleConvSamples = getv(cfg,'poiseuilleConvSamples',5);
poiseuilleRelNuTol = getv(cfg,'poiseuilleRelNuTol',10e-2);
poiseuilleRelQTol = getv(cfg,'poiseuilleRelQTol',5e-2);
poiseuilleR2fitMin = getv(cfg,'poiseuilleR2fitMin',0.90);
poiseuilleSlipTol = getv(cfg,'poiseuilleSlipTol',2e-2);
densitywaveMinSteps = getv(cfg,'densitywaveMinSteps',max(300, round(0.15*nSteps)));
densitywaveConvWindow = getv(cfg,'densitywaveConvWindow',max(60, round(0.08*nSteps)));
densitywaveAtol = getv(cfg,'densitywaveAtol',1e-2);
densitywaveDeltaATol = getv(cfg,'densitywaveDeltaATol',2e-3);
densitywaveRangeATol = getv(cfg,'densitywaveRangeATol',2e-3);
densitywaveDecayFrac = getv(cfg,'densitywaveDecayFrac',0.05);
compressibilityMinSteps = getv(cfg,'compressibilityMinSteps',max(400, round(0.20*nSteps)));
compressibilityConvWindow = getv(cfg,'compressibilityConvWindow',max(80, round(0.10*nSteps)));
compressibilityRelPTol = getv(cfg,'compressibilityRelPTol',1e-2);
compressibilityRelPRangeTol = getv(cfg,'compressibilityRelPRangeTol',2e-2);
compressibilityOccTol = getv(cfg,'compressibilityOccTol',2e-3);
compressibilityPmin = getv(cfg,'compressibilityPmin',1e-6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accumulateMeanFields = getv(cfg,'accumulateMeanFields',[]);
meanFieldsStartStep  = getv(cfg,'meanFieldsStartStep',[]);
meanFieldsStride     = getv(cfg,'meanFieldsStride',[]);
% Cas-specifique + initialisation
switch char(caseType)
    case 'ellipse'
        n = getv(cfg,'n',10000);
        xDrop = getv(cfg,'xDrop',0.5*Lx);
        yDrop = getv(cfg,'yDrop',0.5*Ly);
        aDrop = getv(cfg,'aDrop',1.8);
        bDrop = getv(cfg,'bDrop',0.9);
        thetaEllipseDeg = getv(cfg,'thetaEllipseDeg',0.0);
        dropInitNoise = getv(cfg,'dropInitNoise',0.0);
        dropVelNoise = getv(cfg,'dropVelNoise',0.2);
        nuEffForSigma = getv(cfg,'nuEffForSigma',8.30851);
        rhoEff = n/(Lx*Ly);
        xBoundary = 'specular'; yWallMode = 'specular'; bodyForceX = 0.0; useInterfaceVelocityReorientation = true;
        [x, v, type, r0, gammaInit, nCellsEllipse] = initialize_particles_in_ellipse_main(n, xDrop, yDrop, aDrop, bDrop, thetaEllipseDeg, kBT, dropInitNoise, dropVelNoise, Lx, Ly, Nx, Ny, rngSeed);
        gamma = gammaInit;
        initInfo = struct('nCellsEllipse',nCellsEllipse,'gammaInit',gammaInit,'nType0',sum(type==0),'nType1',sum(type==1));

    case 'diffusion'
        n = getv(cfg,'n',round(gamma*Nc));
        xcTracer = getv(cfg,'xcTracer',0.5*Lx);
        ycTracer = getv(cfg,'ycTracer',0.5*Ly);
        rTracer = getv(cfg,'rTracer',0.25*min(Lx,Ly));
        xBoundary = 'specular'; yWallMode = 'specular'; bodyForceX = 0.0;
        useInterfaceVelocityReorientation = false; useBottomWallWetting = false; useLocalFluidFractionThresholds = false;
        [x, v, type, r0] = initialize_particles_diffusion_main(n, Lx, Ly, kBT, xcTracer, ycTracer, rTracer);
        initInfo = struct('Ntype1',sum(type==1),'rTracer',rTracer);

    case 'densitywave'
        n = getv(cfg,'n',round(gamma*Nc));
        densityWaveAmp = getv(cfg,'densityWaveAmp',0.10);
        densityWaveModeX = getv(cfg,'densityWaveModeX',1);
        xBoundary = 'periodic'; yWallMode = 'specular'; bodyForceX = 0.0;
        useInterfaceVelocityReorientation = false; useBottomWallWetting = false;
        [x, v, type, r0] = initialize_particles_density_wave_main(n, Lx, Ly, Nx, Ny, kBT, densityWaveAmp, densityWaveModeX, rngSeed);
        initInfo = struct('densityWaveAmp',densityWaveAmp,'densityWaveModeX',densityWaveModeX,'nParticles',size(x,1));

    case 'compressibility'
        n = getv(cfg,'n',round(gamma*Nc));
        xBoundary = 'periodic';
        yWallMode = getv(cfg,'yWallMode','specular');
        bodyForceX = 0.0;
        Utop = 0.0; Ubottom = 0.0;
g = getv(cfg,'g',-1.0);
useInterfaceVelocityReorientation = false;
        useBottomWallWetting = false;
        [x, v, type, r0, initCellCounts] = initialize_particles_uniform_by_cell_main(n, Lx, Ly, Nx, Ny, kBT, rngSeed);
        initInfo = struct('nParticles',size(x,1), 'rhoMeanInit',n/(Lx*Ly), 'gamma',gamma, ...
            'minNcell0', min(initCellCounts(:)), 'maxNcell0', max(initCellCounts(:)), ...
            'meanNcell0', mean(initCellCounts(:)), 'stdNcell0', std(double(initCellCounts(:))));

    case 'piston'
        n = getv(cfg,'n',round(gamma*Nc));
        xBoundary = 'periodic';
        yWallMode = getv(cfg,'yWallMode','specular');
        bodyForceX = 0.0;
        Utop = 0.0; Ubottom = 0.0;
        g = getv(cfg,'g',-1.0);
        useInterfaceVelocityReorientation = false;
        useBottomWallWetting = false;
        useMovingPiston = true;
        [x, v, type, r0, initCellCounts] = initialize_particles_uniform_by_cell_main(n, Lx, Ly, Nx, Ny, kBT, rngSeed);
        initInfo = struct('nParticles',size(x,1), 'rhoMeanInit',n/(Lx*Ly), 'gamma',gamma, ...
            'minNcell0', min(initCellCounts(:)), 'maxNcell0', max(initCellCounts(:)), ...
            'meanNcell0', mean(initCellCounts(:)), 'stdNcell0', std(double(initCellCounts(:))), ...
            'pistonY0', pistonY0, 'pistonVy', pistonVy, 'pistonYMin', pistonYMin, ...
            'pistonStopOnMin', pistonStopOnMin, 'pistonFracCut', pistonFracCut, ...
            'pistonAffineReposition', pistonAffineReposition);

    case 'poiseuille'
        n = getv(cfg,'n',round(gamma*Nc));
        fxBody = getv(cfg,'fxBody',getv(cfg,'bodyForceX',0.1));
        poiseuilleInit = getv(cfg,'poiseuilleInit','uniform');
        bodyForceX = fxBody;
        xBoundary = 'periodic'; yWallMode = getv(cfg,'yWallMode','thermalize');
        useInterfaceVelocityReorientation = false; useBottomWallWetting = false; useLocalFluidFractionThresholds = false;
        [x, v, type, r0] = initialize_particles_poiseuille_main(n, Lx, Ly, kBT, fxBody, poiseuilleInit);
        initInfo = struct('fxBody',fxBody,'poiseuilleInit',poiseuilleInit);
    case 'dam'
        % Rupture de barrage :
        % masse initialement confinee dans un rectangle a gauche,
        % puis evolution libre sous gravite dans le domaine.

        xDam = getv(cfg,'xDam',0.35*Lx);   % longueur initiale de la colonne
        hDam = getv(cfg,'hDam',0.60*Ly);   % hauteur initiale de la colonne
        n = getv(cfg,'n',[]);              % si vide, deduit de gamma et du nb de cellules remplies

        bodyForceX = 0.0;
        g = getv(cfg,'g',-10.0);
        kBt = getv(cfg,'kBT',0.001);

        % Presets legacy raisonnables
        xBoundary = 'specular';
        yWallMode = 'thermalize';%'specular';

        % Cas avec surface libre : on garde la topologie de surface active
        useInterfaceVelocityReorientation = true;
        useBottomWallWetting = false;
        noSurfaceCase = false;
interfaceReorientBeta = 0.36;
interfaceReorientFrac = 0.35;
        [x, v, type, r0, gammaInit, nCellsDam] = ...
            initialize_particles_dam_main( ...
                n, gamma, xDam, hDam, kBT, Lx, Ly, Nx, Ny, rngSeed);

        % Si n etait impose manuellement, gammaInit est l'occupation effective
        % dans la zone initialement remplie
        gamma = gammaInit;

        initInfo = struct( ...
            'xDam', xDam, ...
            'hDam', hDam, ...
            'nCellsDam', nCellsDam, ...
            'gammaInit', gammaInit, ...
            'nParticles', size(x,1), ...
            'rhoMeanFilledInit', size(x,1)/(xDam*hDam), ...
            'nType0', sum(type==0), ...
            'nType1', sum(type==1));
    otherwise
        error('caseType inconnu: %s', caseType);
end

% Resolution finale des conditions aux limites explicites.
% Les presets de caseType fixent xBoundary / yWallMode par defaut ; les
% champs boundary_* de cfg peuvent ensuite les surcharger cote par cote.
boundary_left  = string(getv(cfg,'boundary_left',  xBoundary));
boundary_right = string(getv(cfg,'boundary_right', xBoundary));
boundary_bottom = string(getv(cfg,'boundary_bottom', yWallMode));
if useMovingPiston || strcmpi(string(boundary_top_cfg), 'piston')
    boundary_top_default = 'piston';
    useMovingPiston = true;
else
    boundary_top_default = yWallMode;
end
boundary_top = string(getv(cfg,'boundary_top', boundary_top_default));

% Valeurs par defaut des nouveaux drapeaux explicites, avec surcharge cfg.
if isempty(noSurfaceCase)
    noSurfaceCase = ~strcmpi(caseType,'ellipse');
end
if isempty(useLocalTargetsForPiston)
    useLocalTargetsForPiston = useMovingPiston;
end
if isempty(redistributionEnableSurfaceTopology)
    redistributionEnableSurfaceTopology = ~noSurfaceCase;
end
if isempty(redistributionUseActiveVolumeTargets)
    redistributionUseActiveVolumeTargets = useLocalTargetsForPiston;
end
if strlength(string(redistributionBulkOverflowPolicy)) == 0
    if redistributionEnableSurfaceTopology
        redistributionBulkOverflowPolicy = 'nearest_interface';
    else
        redistributionBulkOverflowPolicy = 'none';
    end
end
if isempty(redistributionWallWettingEnabled)
    redistributionWallWettingEnabled = ...
        strcmpi(boundary_left,'wall_wetting') || strcmpi(boundary_right,'wall_wetting') || ...
        strcmpi(boundary_bottom,'wall_wetting') || strcmpi(boundary_top,'wall_wetting') || ...
        logical(useBottomWallWetting);
end
if isempty(redistributionPreserveWettingLayer)
    redistributionPreserveWettingLayer = redistributionWallWettingEnabled;
end
if ~isfinite(redistributionWettingTargetOverride)
    redistributionWettingTargetOverride = bottomWettingTargetOverride;
end
if ~isfinite(redistributionWettingMaxRange)
    redistributionWettingMaxRange = bottomWettingMaxRange;
end

% Valeurs par defaut des diagnostics explicites.
if isempty(computeShapeDiagnostics), computeShapeDiagnostics = strcmpi(caseType,'ellipse'); end
if isempty(computeTracerDiagnostics), computeTracerDiagnostics = strcmpi(caseType,'diffusion'); end
if isempty(computeDensityWaveDiagnostics), computeDensityWaveDiagnostics = strcmpi(caseType,'densitywave'); end
if isempty(densityWaveDiagMode)
    if exist('densityWaveModeX','var')
        densityWaveDiagMode = densityWaveModeX;
    else
        densityWaveDiagMode = 1;
    end
end
if isempty(computePressureSummary)
    computePressureSummary = strcmpi(caseType,'compressibility') || ...
                             strcmpi(caseType,'piston') || ...
                             strcmpi(caseType,'dam') || ...
                             useMovingPiston || strcmpi(boundary_top,'piston');
end
if isempty(computePoiseuilleDiagnostics), computePoiseuilleDiagnostics = strcmpi(caseType,'poiseuille'); end
if isempty(computeViscosityDiagnostics), computeViscosityDiagnostics = computePoiseuilleDiagnostics; end

% Compatibilite descendante : conserver aussi les anciens champs globaux.
if strcmpi(boundary_left, boundary_right)
    xBoundaryResolved = char(boundary_left);
else
    xBoundaryResolved = char(string(xBoundary));
end
if strcmpi(boundary_bottom, boundary_top) && ~strcmpi(boundary_top,'piston')
    yWallModeResolved = char(boundary_bottom);
else
    yWallModeResolved = char(string(yWallMode));
end

exportExplicitBoundaryFields = ~strictLegacyCaseSemantics || strcmpi(caseType,'piston');

x0 = x; v0 = v; type0 = type; r00 = r0;

params = struct();
params.caseType = char(caseType);
params.Lx = Lx; params.Ly = Ly; params.Nx = Nx; params.Ny = Ny; params.a0 = Lx/Nx; params.Nc = Nc;
params.gamma = gamma; params.n = size(x,1); params.dt = dt; params.nSteps = nSteps;
params.alphaDeg = alphaDeg; params.alpha = alpha; params.kBT = kBT; params.g = g; params.bodyForceX = bodyForceX;
params.useThermostat = useThermostat; params.keepMeanFlow = keepMeanFlow;
params.xBoundary = xBoundaryResolved; params.yWallMode = yWallModeResolved; params.Utop = Utop; params.Ubottom = Ubottom; params.wallSigma = wallSigma;
if exportExplicitBoundaryFields
    params.boundary_left = char(boundary_left); params.boundary_right = char(boundary_right);
    params.boundary_bottom = char(boundary_bottom); params.boundary_top = char(boundary_top);
    params.useMovingPiston = useMovingPiston; params.pistonY0 = pistonY0; params.pistonVy = pistonVy; params.pistonYMin = pistonYMin;
    params.pistonStopOnMin = pistonStopOnMin; params.pistonFracCut = pistonFracCut; params.pistonAffineReposition = pistonAffineReposition;
    params.pistonActiveMargin = pistonActiveMargin; params.occAbsFloor = occAbsFloor;
end
params.useIncompressibleRedistribution = useIncompressibleRedistribution; params.redistribAfterCollision = redistribAfterCollision;
params.coef = coef; params.highMode = highMode; params.lowMode = lowMode; params.maxRedistribPasses = maxRedistribPasses;
params.enableMomentumCorrectionPostRedistribution = enableMomentumCorrectionPostRedistribution;
params.useLocalFluidFractionThresholds = useLocalFluidFractionThresholds; params.fluidFracGain = fluidFracGain;
params.lowThrFloor = lowThrFloor; params.highThrFloor = highThrFloor; params.lowThrBulkOverride = lowThrBulkOverride; params.lowThrInterfaceOverride = lowThrInterfaceOverride;
params.useBottomWallWetting = useBottomWallWetting; params.bottomWettingTargetOverride = bottomWettingTargetOverride; params.bottomWettingMaxRange = bottomWettingMaxRange;
if exportExplicitBoundaryFields
    params.noSurfaceCase = noSurfaceCase; params.useLocalTargetsForPiston = useLocalTargetsForPiston;
    params.redistributionEnableSurfaceTopology = redistributionEnableSurfaceTopology;
    params.redistributionUseActiveVolumeTargets = redistributionUseActiveVolumeTargets;
    params.redistributionBulkOverflowPolicy = char(string(redistributionBulkOverflowPolicy));
    params.redistributionWallWettingEnabled = redistributionWallWettingEnabled;
    params.redistributionPreserveWettingLayer = redistributionPreserveWettingLayer;
    params.redistributionWettingTargetOverride = redistributionWettingTargetOverride;
    params.redistributionWettingMaxRange = redistributionWettingMaxRange;
end
params.useInterfaceVelocityReorientation = useInterfaceVelocityReorientation;
params.interfaceReorientBeta = interfaceReorientBeta; params.interfaceReorientFrac = interfaceReorientFrac;
params.useCurvatureWeightedReorientation = useCurvatureWeightedReorientation;
params.curvatureReorientBetaGain = curvatureReorientBetaGain; params.curvatureReorientFracGain = curvatureReorientFracGain;
params.interfaceReorientBetaMax = interfaceReorientBetaMax; params.interfaceReorientFracMax = interfaceReorientFracMax; params.curvatureAngleMax = curvatureAngleMax;
params.stepFcnName = char(string(stepFcnName));
params.useLiquidClosure = useLiquidClosure;
params.useOptimalBetaRepair = useOptimalBetaRepair;
params.betaRepair = betaRepair;
params.betaEOS = betaEOS;
params.Kvirial = Kvirial;
%params.liquidClosureDumpStages = logical(liquidClosureDumpStages);
params.liquidClosureDumpStages = liquidClosureDumpStages;
params.makePlots = makePlots; params.realtimeStride = realtimeStride; params.nplot = nplot; params.dumpStride = dumpStride;
params.quiverStrideX = quiverStrideX; params.quiverStrideY = quiverStrideY; params.quiverScale = quiverScale; params.vortClipSigma = vortClipSigma;
params.diagStride = diagStride; params.viscoStride = viscoStride; params.movAvgWindow = movAvgWindow; params.fitWindowFrac = fitWindowFrac;
params.steadyFracStart = steadyFracStart; params.fitMarginFrac = fitMarginFrac; params.fitSmoothWindowY = fitSmoothWindowY; params.nRadialBins = nRadialBins;
% Les diagnostics explicites sont transmis meme en mode legacy strict :
% ils ne changent pas la physique du coeur, seulement les sorties et
% analyses. Sans cela, le core refactorise retombe a false par defaut.
params.computeShapeDiagnostics = computeShapeDiagnostics; params.computeTracerDiagnostics = computeTracerDiagnostics;
params.computeDensityWaveDiagnostics = computeDensityWaveDiagnostics; params.densityWaveDiagMode = densityWaveDiagMode;
params.computePressureSummary = computePressureSummary; params.computePoiseuilleDiagnostics = computePoiseuilleDiagnostics;
params.computeViscosityDiagnostics = computeViscosityDiagnostics; params.printSummary = printSummary;
params.enableAutoStop = enableAutoStop; params.convCheckEvery = convCheckEvery; params.convPatienceChecks = convPatienceChecks;
params.ellipseMinSteps = ellipseMinSteps; params.ellipseConvWindow = ellipseConvWindow; params.ellipseAtol = ellipseAtol;
params.ellipseRelRTol = ellipseRelRTol; params.ellipseDeltaATol = ellipseDeltaATol; params.ellipseRangeATol = ellipseRangeATol;
params.diffusionMinSteps = diffusionMinSteps; params.diffusionConvWindow = diffusionConvWindow; params.diffusionRelDTol = diffusionRelDTol;
params.diffusionR2linMin = diffusionR2linMin; params.diffusionR2maxFrac = diffusionR2maxFrac;
params.poiseuilleMinSteps = poiseuilleMinSteps; params.poiseuilleConvSamples = poiseuilleConvSamples;
params.poiseuilleRelNuTol = poiseuilleRelNuTol; params.poiseuilleRelQTol = poiseuilleRelQTol; params.poiseuilleR2fitMin = poiseuilleR2fitMin; params.poiseuilleSlipTol = poiseuilleSlipTol;
params.densitywaveMinSteps = densitywaveMinSteps; params.densitywaveConvWindow = densitywaveConvWindow;
params.densitywaveAtol = densitywaveAtol; params.densitywaveDeltaATol = densitywaveDeltaATol; params.densitywaveRangeATol = densitywaveRangeATol; params.densitywaveDecayFrac = densitywaveDecayFrac;
params.compressibilityMinSteps = compressibilityMinSteps; params.compressibilityConvWindow = compressibilityConvWindow;
params.compressibilityRelPTol = compressibilityRelPTol; params.compressibilityRelPRangeTol = compressibilityRelPRangeTol;
params.compressibilityOccTol = compressibilityOccTol; params.compressibilityPmin = compressibilityPmin;
params.coreFcnName = char(string(coreFcnName)); params.runNotes = runNotes;
params.strictLegacyCaseSemantics = strictLegacyCaseSemantics;

if ~isempty(accumulateMeanFields), params.accumulateMeanFields = accumulateMeanFields; end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(meanFieldsStartStep),  params.meanFieldsStartStep  = meanFieldsStartStep;  end
if ~isempty(meanFieldsStride),     params.meanFieldsStride     = meanFieldsStride;     end


if strcmp(caseType,'ellipse'), params.nuEffForSigma = nuEffForSigma; params.rhoEff = rhoEff; end
if strcmp(caseType,'diffusion'), params.xcTracer = xcTracer; params.ycTracer = ycTracer; params.rTracer = rTracer; end
if strcmp(caseType,'densitywave'), params.densityWaveAmp = densityWaveAmp; params.densityWaveModeX = densityWaveModeX; end
if strcmp(caseType,'poiseuille'), params.fxBody = bodyForceX; end
params.neighborsCore = [];
if strcmp(caseType,'dam')
    params.xDam = xDam;
    params.hDam = hDam;
end
runClock = tic;
simOut = feval(coreFcnName, params, x, v, type, r0);
elapsedWallClock_s = toc(runClock);

runOut = struct();
paramsFinal = params;
if isfield(simOut,'params') && ~isempty(simOut.params)
    paramsFinal = simOut.params;
end

runOut.params = paramsFinal;
runOut.simOut = simOut;
runOut.elapsedWallClock_s = elapsedWallClock_s;
runOut.lastCompletedStep = simOut.lastCompletedStep;
runOut.autoStopTriggered = simOut.autoStop.triggered;
runOut.autoStopInfo = simOut.autoStop.info;
runOut.initInfo = initInfo;
runOut.initialState = struct('x',x0,'v',v0,'type',type0,'r0',r00);
runOut.finalState = struct('x',simOut.x,'v',simOut.v,'type',simOut.type,'r0',simOut.r0);
runOut.summary = build_run_summary_main(paramsFinal, simOut, '', '', elapsedWallClock_s);

if isfield(paramsFinal,'lastLiquidClosureDump')
    runOut.lastLiquidClosureDump = paramsFinal.lastLiquidClosureDump;
end

if saveRunResults
    saveInfo = save_mpcd_run_outputs_main(resultsRoot, runLabel, paramsFinal, simOut, x0, v0, type0, r00, simOut.x, simOut.v, simOut.type, simOut.r0, saveInitialState, saveFinalState, saveSimOutStruct, appendSummaryCSV, elapsedWallClock_s);
    runOut.outputDir = saveInfo.resultsDir;
    runOut.summary = saveInfo.summary;
else
    runOut.outputDir = '';
end
end

function val = getv(s, name, default)
if isfield(s,name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end

function [x, v, type, r0, gammaEff, nCellsIn] = initialize_particles_in_ellipse_main(n, xc, yc, a0, b0, thetaDeg, kBT, posNoise, velNoise, Lx, Ly, Nx, Ny, seed)
theta = deg2rad(thetaDeg); ct = cos(theta); st = sin(theta);
dx = Lx / Nx; dy = Ly / Ny;
cellList = zeros(Nx*Ny, 1); m = 0;
for iy = 1:Ny
    yC = (iy - 0.5) * dy;
    for ix = 1:Nx
        xC = (ix - 0.5) * dx;
        xr =  ct*(xC - xc) + st*(yC - yc);
        yr = -st*(xC - xc) + ct*(yC - yc);
        if (xr/a0)^2 + (yr/b0)^2 <= 1
            m = m + 1; cellList(m) = sub2ind([Ny, Nx], iy, ix); end
    end
end
cellList = cellList(1:m); nCellsIn = numel(cellList);
if nCellsIn == 0, error('Ellipse initiale vide.'); end
gammaEff = n / nCellsIn; baseCount = floor(gammaEff); extra = n - baseCount * nCellsIn;
counts = baseCount * ones(nCellsIn, 1); if extra > 0, counts(1:extra) = counts(1:extra) + 1; end
x = zeros(n, 2); v = velNoise * sqrt(max(kBT, 0)) * randn(n, 2); k0 = 0;
for j = 1:nCellsIn
    c = cellList(j); [iy, ix] = ind2sub([Ny, Nx], c); np = counts(j); if np <= 0, continue; end
    x1 = (ix - 1) * dx; y1 = (iy - 1) * dy;
    xi = [x1 + dx * rand(np, 1), y1 + dy * rand(np, 1)];
    if posNoise > 0
        xi = xi + posNoise * [dx * randn(np, 1), dy * randn(np, 1)];
        xi(:, 1) = min(max(xi(:, 1), x1), x1 + dx); xi(:, 2) = min(max(xi(:, 2), y1), y1 + dy);
    end
    x(k0 + 1:k0 + np, :) = xi; k0 = k0 + np;
end
v(:,1) = v(:,1) - mean(v(:,1)); v(:,2) = v(:,2) - mean(v(:,2)); r0 = x; type = split_particle_types_main(n, seed + 12345);
end

function [x, v, type, r0] = initialize_particles_diffusion_main(n, Lx, Ly, kBT, xcTracer, ycTracer, rTracer)
x = [Lx*rand(n,1), Ly*rand(n,1)];
v = sqrt(kBT) * randn(n,2);
v(:,1) = v(:,1) - mean(v(:,1)); v(:,2) = v(:,2) - mean(v(:,2));
rr0 = hypot(x(:,1)-xcTracer, x(:,2)-ycTracer);
type = zeros(n,1,'uint8'); type(rr0 <= rTracer) = uint8(1);
r0 = x;
end

function [x, v, type, r0] = initialize_particles_density_wave_main(n, Lx, Ly, Nx, Ny, kBT, densityWaveAmp, densityWaveModeX, seed)
densityWaveAmp = max(0, min(0.95, densityWaveAmp));
aX = Lx / Nx; aY = Ly / Ny;
xc = ((0:Nx-1)+0.5) * aX;
wcol = 1 + densityWaveAmp * cos(2*pi*densityWaveModeX*xc/Lx);
wcol = max(wcol, 1e-12);
wmat = repmat(wcol, Ny, 1);
w = wmat(:) / sum(wmat(:));
counts = floor(n * w);
remCount = n - sum(counts);
if remCount > 0
    rs = RandStream('mt19937ar','Seed', double(seed) + 34567);
    addIdx = randperm(rs, numel(counts), remCount);
    counts(addIdx) = counts(addIdx) + 1;
end
x = zeros(n,2); k0 = 0;
for c = 1:numel(counts)
    np = counts(c); if np <= 0, continue; end
    [iy, ix] = ind2sub([Ny, Nx], c);
    x1 = (ix - 1) * aX; y1 = (iy - 1) * aY;
    x(k0+1:k0+np,1) = x1 + aX * rand(np,1);
    x(k0+1:k0+np,2) = y1 + aY * rand(np,1);
    k0 = k0 + np;
end
v = sqrt(kBT) * randn(n,2);
v(:,1) = v(:,1) - mean(v(:,1));
v(:,2) = v(:,2) - mean(v(:,2));
type = zeros(n,1,'uint8');
r0 = x;
end

function [x, v, type, r0, initCellCounts] = initialize_particles_uniform_by_cell_main(n, Lx, Ly, Nx, Ny, kBT, seed)
% Initialisation homogène par cellule pour le cas compressibility.
% Répartit d'abord le nombre total de particules presque uniformément sur
% toutes les cellules du maillage, puis tire les positions uniformément
% à l'intérieur de chaque cellule. Cela évite les surcharges locales
% artificielles qui ralentissent fortement la redistribution au démarrage.

if nargin >= 7 && ~isempty(seed)
    rs = RandStream('mt19937ar', 'Seed', double(seed) + 98765);
else
    rs = RandStream.getGlobalStream();
end

Nc = Nx*Ny;
aX = Lx / Nx;
aY = Ly / Ny;

baseCount = floor(n / Nc);
remCount = n - baseCount*Nc;
counts = baseCount * ones(Nc,1);

if remCount > 0
    addIdx = randperm(rs, Nc, remCount);
    counts(addIdx) = counts(addIdx) + 1;
end

x = zeros(n,2);
k0 = 0;
for c = 1:Nc
    np = counts(c);
    if np <= 0, continue; end
    [iy, ix] = ind2sub([Ny, Nx], c);
    x1 = (ix - 1) * aX;
    y1 = (iy - 1) * aY;
    x(k0+1:k0+np,1) = x1 + aX * rand(rs, np, 1);
    x(k0+1:k0+np,2) = y1 + aY * rand(rs, np, 1);
    k0 = k0 + np;
end

v = sqrt(kBT) * randn(rs, n, 2);
v(:,1) = v(:,1) - mean(v(:,1));
v(:,2) = v(:,2) - mean(v(:,2));
type = zeros(n,1,'uint8');
r0 = x;
initCellCounts = reshape(counts, [Ny, Nx]);
end

function [x, v, type, r0] = initialize_particles_poiseuille_main(n, Lx, Ly, kBT, fxBody, poiseuilleInit)
x = [Lx*rand(n,1), Ly*rand(n,1)];
v = sqrt(kBT) * randn(n,2);
if strcmpi(poiseuilleInit, 'parabolicguess')
    y = x(:,2); nuGuess = 1.0; v(:,1) = v(:,1) + 0.5*fxBody/nuGuess * y .* (Ly - y);
end
type = zeros(n,1,'uint8');
r0 = x;
end
function [x, v, type, r0, gammaEff, nCellsIn] = ...
    initialize_particles_dam_main(nWanted, gammaTarget, xDam, hDam, kBT, Lx, Ly, Nx, Ny, seed)

% Initialise une "colonne d'eau" rectangulaire dans
%   0 <= x <= xDam
%   0 <= y <= hDam
% avec remplissage quasi uniforme par cellule.

if nargin < 10 || isempty(seed)
    seed = 1;
end
rs = RandStream('mt19937ar','Seed', double(seed) + 24680);

dx = Lx / Nx;
dy = Ly / Ny;

% Liste des cellules dont le centre est dans le reservoir initial
cellList = zeros(Nx*Ny,1);
m = 0;
for iy = 1:Ny
    yC = (iy - 0.5) * dy;
    for ix = 1:Nx
        xC = (ix - 0.5) * dx;
        if (xC <= xDam) && (yC <= hDam)
            m = m + 1;
            cellList(m) = sub2ind([Ny, Nx], iy, ix);
        end
    end
end
cellList = cellList(1:m);
nCellsIn = numel(cellList);

if nCellsIn == 0
    error('La zone initiale dam est vide : augmenter xDam et/ou hDam.');
end

% Si le nombre de particules n'est pas impose, on vise gammaTarget
if isempty(nWanted) || ~isscalar(nWanted) || nWanted <= 0
    nWanted = round(gammaTarget * nCellsIn);
end

gammaEff = nWanted / nCellsIn;

baseCount = floor(gammaEff);
extra = nWanted - baseCount * nCellsIn;

counts = baseCount * ones(nCellsIn,1);
if extra > 0
    addIdx = randperm(rs, nCellsIn, extra);
    counts(addIdx) = counts(addIdx) + 1;
end

x = zeros(nWanted, 2);
k0 = 0;

for j = 1:nCellsIn
    c = cellList(j);
    np = counts(j);
    if np <= 0, continue; end

    [iy, ix] = ind2sub([Ny, Nx], c);
    x1 = (ix - 1) * dx;
    y1 = (iy - 1) * dy;

    x(k0+1:k0+np, 1) = x1 + dx * rand(rs, np, 1);
    x(k0+1:k0+np, 2) = y1 + dy * rand(rs, np, 1);

    k0 = k0 + np;
end

% Vitesses thermiques initiales
v = sqrt(kBT) * randn(rs, nWanted, 2);
v(:,1) = v(:,1) - mean(v(:,1));
v(:,2) = v(:,2) - mean(v(:,2));

% Deux sous-populations juste pour visualisation/eventuel suivi
type = split_particle_types_main(nWanted, seed + 54321);

r0 = x;
end
function type = split_particle_types_main(n, seed)
stream = RandStream('mt19937ar', 'Seed', seed);
perm = randperm(stream, n);
nhalf = floor(n/2);
type = zeros(n,1,'uint8');
type(perm(nhalf+1:end)) = uint8(1);
end

function saveInfo = save_mpcd_run_outputs_main(resultsRoot, runLabel, params, simOut, x0, v0, type0, r00, xF, vF, typeF, r0F, saveInitialState, saveFinalState, saveSimOutStruct, appendSummaryCSV, elapsedWallClock_s)
if ~exist(resultsRoot, 'dir'), mkdir(resultsRoot); end
stamp = datestr(now, 'yyyymmdd_HHMMSSFFF');
if isempty(runLabel), runLabel = char(params.caseType); else, runLabel = char(runLabel); end
safeLabel = regexprep(runLabel, '[^a-zA-Z0-9_\-]', '_');
runDirName = sprintf('%s_%s', safeLabel, stamp);
resultsDir = fullfile(resultsRoot, runDirName);
mkdir(resultsDir);

summary = build_run_summary_main(params, simOut, resultsDir, runDirName, elapsedWallClock_s);
save(fullfile(resultsDir, 'summary.mat'), 'summary');

paramsSave = params; paramsSave.resultsDir = resultsDir; paramsSave.resultsTimestamp = stamp; paramsSave.elapsedWallClock_s = elapsedWallClock_s; %#ok<NASGU>
save(fullfile(resultsDir, 'params.mat'), 'paramsSave');

series = build_series_payload_main(params, simOut);
save(fullfile(resultsDir, 'series.mat'), '-struct', 'series', '-v7.3');

if saveInitialState
    initialState = struct('x0',x0,'v0',v0,'type0',type0,'r00',r00); %#ok<NASGU>
    save(fullfile(resultsDir,'initial_state.mat'),'initialState','-v7.3');
end
if saveFinalState
    finalState = struct('x',xF,'v',vF,'type',typeF,'r0',r0F); %#ok<NASGU>
    save(fullfile(resultsDir,'final_state.mat'),'finalState','-v7.3');
end
if saveSimOutStruct
    save(fullfile(resultsDir,'simOut.mat'),'simOut','-v7.3');
end

write_run_readme_main(fullfile(resultsDir,'README.txt'), summary, params);
if appendSummaryCSV
    append_summary_csv_main(fullfile(resultsRoot,'campaign_index.csv'), summary);
end
saveInfo = struct('resultsDir',resultsDir,'summary',summary);
end

function summary = build_run_summary_main(params, simOut, resultsDir, runDirName, elapsedWallClock_s)
summary = struct();
summary.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
summary.runDirName = runDirName;
summary.resultsDir = resultsDir;
summary.caseType = char(params.caseType);
summary.coreFcnName = params.coreFcnName;
summary.elapsedWallClock_s = elapsedWallClock_s;
summary.lastCompletedStep = safe_getfield_main(simOut,'lastCompletedStep',NaN);

summary.autoStopTriggered = false; summary.autoStopReason=''; summary.autoStopMetric=''; summary.autoStopValue=NaN;
if isfield(simOut,'autoStop') && isstruct(simOut.autoStop)
    if isfield(simOut.autoStop,'triggered'), summary.autoStopTriggered = logical(simOut.autoStop.triggered); end
    if isfield(simOut.autoStop,'info') && isstruct(simOut.autoStop.info)
        info = simOut.autoStop.info;
        if isfield(info,'reason'), summary.autoStopReason = char(string(info.reason)); end
        if isfield(info,'metric'), summary.autoStopMetric = char(string(info.metric)); end
        if isfield(info,'value') && isnumeric(info.value) && isscalar(info.value), summary.autoStopValue = info.value; end
    end
end

summary.nRequested = params.n;
summary.gamma = params.gamma;
summary.Lx = params.Lx; summary.Ly = params.Ly; summary.Nx = params.Nx; summary.Ny = params.Ny;
summary.dt = params.dt; summary.nSteps = params.nSteps; summary.alphaDeg = params.alphaDeg; summary.kBT = params.kBT; summary.g = params.g; summary.bodyForceX = params.bodyForceX;
summary.useIncompressibleRedistribution = logical(params.useIncompressibleRedistribution);
summary.useInterfaceVelocityReorientation = logical(params.useInterfaceVelocityReorientation);
summary.meanOccStd = nanmean(simOut.occStd);
summary.finalOccStd = safe_last_numeric_main(simOut.occStd);
summary.meanOutBand = nanmean(simOut.outBand);
summary.finalOutBand = safe_last_numeric_main(simOut.outBand);
summary.finalEnergy = safe_last_numeric_main(simOut.E);
summary.nParticlesFinal = size(simOut.x,1);

% Champs cas-specifiques (preinitialises)
summary.A_final = NaN; summary.A_eq = NaN; summary.tauA = NaN; summary.sigmaProxyA = NaN; summary.Rx_final = NaN; summary.Ry_final = NaN;
summary.Dglobal = NaN; summary.DinstMean = NaN; summary.fracDiskMean = NaN; summary.R2final = NaN;
summary.nuCurve_mean = NaN; summary.nuTau_mean = NaN; summary.qFlow_mean = NaN; summary.R2fit_mean = NaN; summary.slipBot_final = NaN; summary.slipTop_final = NaN; summary.uxMean_final = NaN;
summary.densityAmpRel_final = NaN; summary.densityAmpRel_eq = NaN; summary.tauDensity = NaN; summary.densityColStd_final = NaN;
summary.rhoMean = NaN; summary.pBot_mean = NaN; summary.pTop_mean = NaN; summary.pMean_mean = NaN; summary.pMean_std = NaN;

switch lower(char(params.caseType))
    case 'ellipse'
        if isfield(simOut,'Ashape'), summary.A_final = safe_last_numeric_main(simOut.Ashape); end
        if isfield(simOut,'RxCtrl'), summary.Rx_final = safe_last_numeric_main(simOut.RxCtrl); end
        if isfield(simOut,'RyCtrl'), summary.Ry_final = safe_last_numeric_main(simOut.RyCtrl); end
        if isfield(simOut,'sigma') && isstruct(simOut.sigma)
            s = simOut.sigma;
            if isfield(s,'Aeq'), summary.A_eq = s.Aeq; end
            if isfield(s,'tauA'), summary.tauA = s.tauA; end
            if isfield(s,'sigmaProxyA'), summary.sigmaProxyA = s.sigmaProxyA; end
        end

    case 'diffusion'
        if isfield(simOut,'diffusion') && isstruct(simOut.diffusion)
            d = simOut.diffusion;
            if isfield(d,'Dglobal'), summary.Dglobal = d.Dglobal; end
            if isfield(d,'DinstMean'), summary.DinstMean = d.DinstMean; end
            if isfield(d,'fracDiskMean'), summary.fracDiskMean = d.fracDiskMean; end
            if isfield(d,'diffDiag') && ~isempty(d.diffDiag), summary.R2final = d.diffDiag(end,6); end
        end

    case 'densitywave'
        if isfield(simOut,'densitywave') && isstruct(simOut.densitywave)
            d = simOut.densitywave;
            if isfield(d,'densityWaveDiag') && ~isempty(d.densityWaveDiag)
                summary.densityAmpRel_final = d.densityWaveDiag(end,2);
                summary.densityColStd_final = d.densityWaveDiag(end,4);
            end
            if isfield(d,'ArelEq'), summary.densityAmpRel_eq = d.ArelEq; end
            if isfield(d,'tauDensity'), summary.tauDensity = d.tauDensity; end
        end

    case {'compressibility','piston'}
        if isfield(simOut,'compressibility') && isstruct(simOut.compressibility)
            c = simOut.compressibility;
            if isfield(c,'rhoMean'),   summary.rhoMean = c.rhoMean; end
            if isfield(c,'pBot_mean'), summary.pBot_mean = c.pBot_mean; end
            if isfield(c,'pTop_mean'), summary.pTop_mean = c.pTop_mean; end
            if isfield(c,'pMean_mean'), summary.pMean_mean = c.pMean_mean; end
            if isfield(c,'pMean_std'),  summary.pMean_std  = c.pMean_std; end
        end

    case 'poiseuille'
        if isfield(simOut,'poiseuille') && isstruct(simOut.poiseuille)
            p = simOut.poiseuille;
            if isfield(p,'nuCurve_mean'), summary.nuCurve_mean = p.nuCurve_mean; end
            if isfield(p,'nuTau_mean'), summary.nuTau_mean = p.nuTau_mean; end
            if isfield(p,'qFlow_mean'), summary.qFlow_mean = p.qFlow_mean; end
            if isfield(p,'viscoDiag') && ~isempty(p.viscoDiag) && size(p.viscoDiag,2)>=9, summary.R2fit_mean = nanmean(p.viscoDiag(:,9)); end
            if isfield(p,'poiseuilleDiag') && ~isempty(p.poiseuilleDiag)
                pd = p.poiseuilleDiag;
                if size(pd,2)>=1, summary.slipBot_final = pd(end,1); end
                if size(pd,2)>=2, summary.slipTop_final = pd(end,2); end
                if size(pd,2)>=3, summary.uxMean_final = pd(end,3); end
            end
        end
end
end

function series = build_series_payload_main(params, simOut)
series = struct();
series.caseType = params.caseType;
series.lastCompletedStep = safe_getfield_main(simOut,'lastCompletedStep',NaN);
series.E = simOut.E; series.occStd = simOut.occStd; series.outBand = simOut.outBand; series.redistDiag = simOut.redistDiag; series.wallDiag = simOut.wallDiag;
series.msdAll = simOut.msdAll; series.msdType0 = simOut.msdType0; series.msdType1 = simOut.msdType1; series.nType0 = simOut.nType0; series.nType1 = simOut.nType1;

switch lower(char(params.caseType))
    case 'ellipse'
        series.xCM = safe_getfield_main(simOut,'xCM',[]);
        series.yCM = safe_getfield_main(simOut,'yCM',[]);
        series.rRMS = safe_getfield_main(simOut,'rRMS',[]);
        series.rMax = safe_getfield_main(simOut,'rMax',[]);
        series.RxCtrl = safe_getfield_main(simOut,'RxCtrl',[]);
        series.RyCtrl = safe_getfield_main(simOut,'RyCtrl',[]);
        series.Ashape = safe_getfield_main(simOut,'Ashape',[]);
        series.sigma = safe_getfield_main(simOut,'sigma',struct());

    case 'diffusion'
        if isfield(simOut,'diffusion') && isstruct(simOut.diffusion)
            d = simOut.diffusion;
            series.diffDiag = safe_getfield_main(d,'diffDiag',[]);
            series.radialDiag = safe_getfield_main(d,'radialDiag',[]);
            series.rBins = safe_getfield_main(d,'rBins',[]);
            series.cRadial = safe_getfield_main(d,'cRadial',[]);
        end

    case 'densitywave'
        if isfield(simOut,'densitywave') && isstruct(simOut.densitywave)
            d = simOut.densitywave;
            series.densityWaveDiag = safe_getfield_main(d,'densityWaveDiag',[]);
            series.Arel_f = safe_getfield_main(d,'Arel_f',[]);
            series.ArelEq = safe_getfield_main(d,'ArelEq',NaN);
            series.tauDensity = safe_getfield_main(d,'tauDensity',NaN);
        end

    case {'compressibility','piston'}
        if isfield(simOut,'compressibility') && isstruct(simOut.compressibility)
            c = simOut.compressibility;
            series.pressureDiag = safe_getfield_main(simOut,'pressureDiag',[]);
            series.compressibility = c;
        end

    case 'poiseuille'
        if isfield(simOut,'poiseuille') && isstruct(simOut.poiseuille)
            p = simOut.poiseuille;
            series.poiseuilleDiag = safe_getfield_main(p,'poiseuilleDiag',[]);
            series.viscoDiag = safe_getfield_main(p,'viscoDiag',[]);
            series.uxProf = safe_getfield_main(p,'uxProf',[]);
            series.yCenters = safe_getfield_main(p,'yCenters',[]);
            series.uAna = safe_getfield_main(p,'uAna',[]);
        end
end
end

function write_run_readme_main(filename, summary, params)
fid = fopen(filename,'w'); if fid<0, return; end
clean = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Run MPCD\n');
fprintf(fid,'Case: %s\n', summary.caseType);
fprintf(fid,'Core: %s\n', summary.coreFcnName);
fprintf(fid,'Timestamp: %s\n', summary.timestamp);
fprintf(fid,'Elapsed wall clock [s]: %.6g\n', summary.elapsedWallClock_s);
fprintf(fid,'Last completed step: %g\n', summary.lastCompletedStep);
fprintf(fid,'Auto-stop: %d\n', summary.autoStopTriggered);
fprintf(fid,'Auto-stop reason: %s\n', summary.autoStopReason);
fprintf(fid,'Auto-stop metric: %s\n', summary.autoStopMetric);
fprintf(fid,'Auto-stop value: %.6g\n', summary.autoStopValue);
fprintf(fid,'\nParams:\n');
try
    fns = fieldnames(params);
    for i=1:numel(fns)
        val = params.(fns{i});
        if isnumeric(val) && isscalar(val)
            fprintf(fid,'  %s = %.12g\n', fns{i}, val);
        elseif islogical(val) && isscalar(val)
            fprintf(fid,'  %s = %d\n', fns{i}, val);
        elseif ischar(val) || (isstring(val) && isscalar(val))
            fprintf(fid,'  %s = %s\n', fns{i}, char(string(val)));
        end
    end
catch
end
end

function append_summary_csv_main(csvPath, summary)
row = summary_struct_to_row_main(summary);
headers = fieldnames(row);

if ~exist(csvPath,'file')
    fid = fopen(csvPath,'w'); if fid<0, return; end
    clean = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', strjoin(headers, ','));
    fprintf(fid, '%s\n', row_to_csv_line_main(row, headers));
else
    fid = fopen(csvPath,'a'); if fid<0, return; end
    clean = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', row_to_csv_line_main(row, headers));
end
end

function row = summary_struct_to_row_main(summary)
row = summary;
fns = fieldnames(row);
for i=1:numel(fns)
    val = row.(fns{i});
    if islogical(val) && isscalar(val), row.(fns{i}) = double(val); end
end
end

function line = row_to_csv_line_main(row, headers)
vals = cell(1,numel(headers));
for i=1:numel(headers)
    v = row.(headers{i});
    if isnumeric(v) && isscalar(v)
        vals{i} = sprintf('%.12g', v);
    else
        vals{i} = csv_escape_main(char(string(v)));
    end
end
line = strjoin(vals, ',');
end

function s = csv_escape_main(s)
s = strrep(s, '"', '""');
if contains(s,',') || contains(s,'"')
    s = ['"' s '"'];
end
end

function v = safe_getfield_main(s, fld, default)
if isstruct(s) && isfield(s,fld), v = s.(fld); else, v = default; end
end

function v = safe_last_numeric_main(x)
if isempty(x) || ~isnumeric(x), v = NaN; else, v = x(end); end
end
