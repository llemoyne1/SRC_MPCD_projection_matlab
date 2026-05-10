function out = mpcd_incompressible_simulator_pistondiag(params, x, v, type, r0)
%MPCD_INCOMPRESSIBLE_SIMULATOR
% Boucle temporelle, diagnostics et affichages MPCD.
% Ce simulateur s'appuie sur un coeur "un pas de temps" appele a chaque
% iteration. Le coeur physique est isole dans mpcd_incompressible_step.
%
% Usage:
%   out = mpcd_incompressible_simulator(params, x, v, type, r0)
%
% Entree/sortie compatibles avec l'ancien noyau monolithique, de sorte que
% le runner principal n'a besoin que de pointer coreFcnName sur ce fichier.

if nargin < 3
    error('mpcd_incompressible_simulator requiert au moins params, x, v.');
end
if nargin < 4 || isempty(type)
    type = zeros(size(x,1),1,'uint8');
end
if nargin < 5 || isempty(r0)
    r0 = x;
end
if size(x,1) ~= size(v,1) || size(x,1) ~= numel(type) || size(x,1) ~= size(r0,1)
    error('x, v, type et r0 doivent avoir le meme nombre de particules.');
end

params = prepare_simulator_params(params, x, type);
params = prepare_step_params(params, x);
stepFcnName = getv_local(params, 'stepFcnName', 'mpcd_incompressible_step_liquidclosure');

Lx = params.Lx; Ly = params.Ly; %#ok<NASGU>
Nx = params.Nx; Ny = params.Ny;
dt = params.dt; nSteps = params.nSteps;
movAvgWindow = params.movAvgWindow;
fitWindowFrac = params.fitWindowFrac;
fitMarginFrac = params.fitMarginFrac;
fitSmoothWindowY = params.fitSmoothWindowY;
steadyFracStart = params.steadyFracStart;
nRadialBins = params.nRadialBins;
xcTracer = params.xcTracer; ycTracer = params.ycTracer; rTracer = params.rTracer;
densityWaveDiagMode = params.densityWaveDiagMode;
bodyForceX = params.bodyForceX;
makePlots = params.makePlots;
nplot = params.nplot;
printSummary = params.printSummary;
computeShapeDiagnostics = params.computeShapeDiagnostics;
computeTracerDiagnostics = params.computeTracerDiagnostics;
computeDensityWaveDiagnostics = params.computeDensityWaveDiagnostics;
computePressureSummary = params.computePressureSummary;
computePistonDiagnostics = params.computePistonDiagnostics;
computePoiseuilleDiagnostics = params.computePoiseuilleDiagnostics;
computeViscosityDiagnostics = params.computeViscosityDiagnostics;
viscoStride = params.viscoStride;
accumulateMeanFields = params.accumulateMeanFields;
meanFieldsStartStep = params.meanFieldsStartStep;
meanFieldsStride = params.meanFieldsStride;
accumulateMeanFields = params.accumulateMeanFields;
meanFieldsStride = params.meanFieldsStride;
meanFieldsStartStep = params.meanFieldsStartStep;

itStart = 1;
lastCompletedStep = 0;
nSeg = nSteps - itStart + 1;
if nSeg <= 0
    error('nSteps doit etre >= itStart');
end
segIter = (itStart:nSteps).';
timeSeg = (segIter - 1)*dt;

% --- diagnostics communs
E = nan(nSeg,1);
occStd = nan(nSeg,1);
outBand = nan(nSeg,1);
redistDiag = nan(nSeg,10);
redistribStats = nan(nSeg,6);
wallDiag = nan(nSeg,7); % [nBot nTop dEwall dPxBot dPxTop dPyBot dPyTop]
msdAll = nan(nSeg,1);
msdType0 = nan(nSeg,1);
msdType1 = nan(nSeg,1);
nType0 = nan(nSeg,1);
nType1 = nan(nSeg,1);

% --- ellipse
rRMS = nan(nSeg,1);
rMax = nan(nSeg,1);
xCM = nan(nSeg,1);
yCM = nan(nSeg,1);
RxCtrl = nan(nSeg,1);
RyCtrl = nan(nSeg,1);
Ashape = nan(nSeg,1);

% --- diffusion
traceDiag = nan(nSeg,13);
radialDiag = nan(nSeg,4);

% --- poiseuille
poiseuilleDiag = nan(nSeg,11);
viscoDiag = nan(nSeg,12);

% --- density wave / compressibilite
densityWaveDiag = nan(nSeg,8);
pressureDiag = nan(nSeg,8); % [dPyBot dPyTop pBotInst pTopInst pMeanInst pMeanCum occStd outBand]
pistonDiag = nan(nSeg,17); % [yTop Hactive compression activeVol activeFracMean nActiveEff gammaGeom gammaTarget rhoMeanActive capacityRatio pBot pTop pMean PkinActive PvirActive PtotActive overfillActive]

autoStop = struct('triggered', false, 'info', struct('reason','', 'metric', '', 'value', NaN, 'step', NaN, 'candidateCount', 0));
lastLiquidClosure = struct();

meanFieldAcc = init_mean_field_accumulator(params);

meanFieldSum = struct();
meanFieldNSamples = 0;
if accumulateMeanFields
    meanFieldTemplate = compute_liquidclosure_mean_fields(x, v, params);
    z = zeros(size(meanFieldTemplate.Ux));
    meanFieldSum.Ux = z;
    meanFieldSum.Uy = z;
    meanFieldSum.rho = z;
    meanFieldSum.T = z;
    meanFieldSum.P = z;
    meanFieldSum.Pkin = z;
    meanFieldSum.Pvir = z;
    meanFieldSum.rhoTarget = z;
    meanFieldGeom = rmfield(meanFieldTemplate, {'Ux','Uy','rho','T','P','Pkin','Pvir','rhoTarget'});
else
    meanFieldGeom = struct();
end

for it = itStart:nSteps
    iseg = it - itStart + 1;

    [x, v, type, r0, stepInfo, params] = feval(stepFcnName, params, x, v, type, r0, it);
    if isfield(stepInfo,'liquidClosure') && ~isempty(stepInfo.liquidClosure)
        lastLiquidClosure = stepInfo.liquidClosure;
    end
    if accumulateMeanFields && (it >= meanFieldsStartStep) && (mod(it - meanFieldsStartStep, max(1,meanFieldsStride)) == 0)
        meanFieldAcc = accumulate_mean_field_snapshot(meanFieldAcc, params, x, v);
    end

    wallInfo = stepInfo.wallInfo;
    Ncell = stepInfo.Ncell;
    diag = stepInfo.redistDiag;
    stats = stepInfo.redistribStats;
    gammaCell = stepInfo.gammaCell;

    wallDiag(iseg,:) = [wallInfo.nBot, wallInfo.nTop, wallInfo.dEwall, wallInfo.dPxBot, wallInfo.dPxTop, wallInfo.dPyBot, wallInfo.dPyTop];
    if iseg == 1
        pMeanCum = stepInfo.pMeanInst;
    else
        pMeanCum = mean([pressureDiag(1:iseg-1,5); stepInfo.pMeanInst], 'omitnan');
    end
    pressureDiag(iseg,:) = [stepInfo.pRow(1:5), pMeanCum, NaN, NaN];

    % Diagnostics communs
    E(iseg) = 0.5 * sum(v(:,1).^2 + v(:,2).^2);
    occStd(iseg) = diag(1);
    outBand(iseg) = diag(2);
    redistDiag(iseg,:) = diag;
    redistribStats(iseg,:) = stats;
    pressureDiag(iseg,7) = occStd(iseg);
    pressureDiag(iseg,8) = outBand(iseg);

    if computePistonDiagnostics
        pistonDiag(iseg,:) = piston_runtime_diagnostics(params, x, v, stepInfo);
    end

    dr2 = sum((x - r0).^2, 2);
    msdAll(iseg) = mean(dr2);
    mask0 = (type == uint8(0));
    mask1 = (type == uint8(1));
    nType0(iseg) = sum(mask0);
    nType1(iseg) = sum(mask1);
    if any(mask0), msdType0(iseg) = mean(dr2(mask0)); end
    if any(mask1), msdType1(iseg) = mean(dr2(mask1)); end

    lastCompletedStep = it;

    if accumulateMeanFields && (it >= meanFieldsStartStep) && (mod(it - meanFieldsStartStep, meanFieldsStride) == 0)
        Gmean = compute_liquidclosure_mean_fields(x, v, params);
        meanFieldSum.Ux = meanFieldSum.Ux + Gmean.Ux;
        meanFieldSum.Uy = meanFieldSum.Uy + Gmean.Uy;
        meanFieldSum.rho = meanFieldSum.rho + Gmean.rho;
        meanFieldSum.T = meanFieldSum.T + Gmean.T;
        meanFieldSum.P = meanFieldSum.P + Gmean.P;
        meanFieldSum.Pkin = meanFieldSum.Pkin + Gmean.Pkin;
        meanFieldSum.Pvir = meanFieldSum.Pvir + Gmean.Pvir;
        meanFieldSum.rhoTarget = meanFieldSum.rhoTarget + Gmean.rhoTarget;
        meanFieldNSamples = meanFieldNSamples + 1;
    end

    if computeShapeDiagnostics
        xcm = mean(x(:,1)); ycm = mean(x(:,2));
        xCM(iseg) = xcm; yCM(iseg) = ycm;
        dx = x(:,1) - xcm; dy = x(:,2) - ycm;
        rr = hypot(dx, dy);
        rRMS(iseg) = sqrt(mean(rr.^2));
        rMax(iseg) = max(rr);
        RxCtrl(iseg) = sqrt(2*mean(dx.^2));
        RyCtrl(iseg) = sqrt(2*mean(dy.^2));
        Ashape(iseg) = (RxCtrl(iseg) - RyCtrl(iseg)) / max(RxCtrl(iseg) + RyCtrl(iseg), eps);
    end

    if computeTracerDiagnostics
        [traceDiag(iseg,:), radialDiag(iseg,:)] = tracer_diffusion_diagnostics( ...
            x, v, type, xcTracer, ycTracer, rTracer, nRadialBins, params.Lx, params.Ly, dt, iseg);
    end

    if computeDensityWaveDiagnostics
        densityWaveDiag(iseg,:) = density_wave_diagnostics(Ncell, Nx, Ny, gammaCell, densityWaveDiagMode);
    end

    if computePistonDiagnostics
    E_f = movmean_omitnan(E, movAvgWindow);
    rhoStd_f = movmean_omitnan(redistDiag(:,1), movAvgWindow);
    outBand_f = movmean_omitnan(redistDiag(:,2), movAvgWindow);
    yTop_f = movmean_omitnan(pistonDiag(:,1), movAvgWindow);
    Hactive_f = movmean_omitnan(pistonDiag(:,2), movAvgWindow);
    compression_f = movmean_omitnan(pistonDiag(:,3), movAvgWindow);
    gammaGeom_f = movmean_omitnan(pistonDiag(:,7), movAvgWindow);
    gammaTarget_f = movmean_omitnan(pistonDiag(:,8), movAvgWindow);
    rhoMeanActive_f = movmean_omitnan(pistonDiag(:,9), movAvgWindow);
    capacityRatio_f = movmean_omitnan(pistonDiag(:,10), movAvgWindow);
    pBot_f = movmean_omitnan(pistonDiag(:,11), movAvgWindow);
    pTop_f = movmean_omitnan(pistonDiag(:,12), movAvgWindow);
    pMean_f = movmean_omitnan(pistonDiag(:,13), movAvgWindow);
    PkinActive_f = movmean_omitnan(pistonDiag(:,14), movAvgWindow);
    PvirActive_f = movmean_omitnan(pistonDiag(:,15), movAvgWindow);
    PtotActive_f = movmean_omitnan(pistonDiag(:,16), movAvgWindow);
    overfillActive_f = movmean_omitnan(pistonDiag(:,17), movAvgWindow);
    [Beff, kEff, dPdrhoEff, BfromdPdrho] = piston_effective_compressibility(pistonDiag);
    Beff_f = movmean_omitnan(Beff, movAvgWindow);
    kEff_f = movmean_omitnan(kEff, movAvgWindow);
    dPdrhoEff_f = movmean_omitnan(dPdrhoEff, movAvgWindow);
    BfromdPdrho_f = movmean_omitnan(BfromdPdrho, movAvgWindow);
    out.piston = struct('pistonDiag',pistonDiag,'yTop_f',yTop_f,'Hactive_f',Hactive_f,...
        'compression_f',compression_f,'gammaGeom_f',gammaGeom_f,'gammaTarget_f',gammaTarget_f,'rhoMeanActive_f',rhoMeanActive_f,'capacityRatio_f',capacityRatio_f,...
        'pBot_f',pBot_f,'pTop_f',pTop_f,'pMean_f',pMean_f,'PkinActive_f',PkinActive_f,'PvirActive_f',PvirActive_f,'PtotActive_f',PtotActive_f,'overfillActive_f',overfillActive_f,...
        'Beff',Beff,'kEff',kEff,'dPdrhoEff',dPdrhoEff,'BfromdPdrho',BfromdPdrho,...
        'Beff_f',Beff_f,'kEff_f',kEff_f,'dPdrhoEff_f',dPdrhoEff_f,'BfromdPdrho_f',BfromdPdrho_f,'E_f',E_f,'rhoStd_f',rhoStd_f,'outBand_f',outBand_f);
    if printSummary
        fprintf('\n=== Resume piston ===\n');
        fprintf('yTop final                 = %.6g\n', last_or_nan(pistonDiag(:,1)));
        fprintf('H active final             = %.6g\n', last_or_nan(pistonDiag(:,2)));
        fprintf('compression H/Ly final     = %.6g\n', last_or_nan(pistonDiag(:,3)));
        fprintf('gammaGeom final            = %.6g\n', last_or_nan(pistonDiag(:,7)));
        fprintf('gammaTarget final          = %.6g\n', last_or_nan(pistonDiag(:,8)));
        fprintf('rhoMean active final       = %.6g\n', last_or_nan(pistonDiag(:,9)));
        fprintf('capacity ratio final       = %.6g\n', last_or_nan(pistonDiag(:,10)));
        fprintf('pression bas / haut final  = %.6g / %.6g\n', last_or_nan(pistonDiag(:,11)), last_or_nan(pistonDiag(:,12)));
        fprintf('pression moyenne final     = %.6g\n', last_or_nan(pistonDiag(:,13)));
        fprintf('Pkin / Pvir / Ptot actives = %.6g / %.6g / %.6g\n', last_or_nan(pistonDiag(:,14)), last_or_nan(pistonDiag(:,15)), last_or_nan(pistonDiag(:,16)));
        fprintf('overfill actif final       = %.6g\n', last_or_nan(pistonDiag(:,17)));
        fprintf('B_eff final                = %.6g\n', last_or_nan(Beff));
        fprintf('k_eff final                = %.6g\n', last_or_nan(kEff));
        fprintf('dPtot/drho eff final       = %.6g\n', last_or_nan(dPdrhoEff));
        fprintf('rho*dPtot/drho final       = %.6g\n', last_or_nan(BfromdPdrho));
        fprintf('======================\n');
    end
end

if computePoiseuilleDiagnostics
        poiseuilleDiag(iseg,:) = poiseuille_diagnostics_fast(x, v, params.Ly, Ny, wallInfo, dt, params.Lx);
        if computeViscosityDiagnostics && (mod(it, viscoStride) == 0 || it == itStart || it == nSteps)
            [uxProf, ~, ~] = mean_profiles_y_fast(x, v, params.Lx, params.Ly, Nx, Ny);
            yCenters = ((0:Ny-1)+0.5)' * (params.Ly/Ny);
            uxProfFit = movmean_omitnan(uxProf, fitSmoothWindowY);
            j0 = max(2, round(fitMarginFrac*Ny));
            j1 = min(Ny-1, round((1-fitMarginFrac)*Ny));
            yy = yCenters(j0:j1);
            uu = uxProfFit(j0:j1);
            if numel(yy) >= 6
                p2 = polyfit(yy, uu, 2);
                uuFit = polyval(p2, yy);
                ssRes = sum((uu - uuFit).^2);
                ssTot = sum((uu - mean(uu)).^2);
                R2 = 1 - ssRes/max(ssTot, eps);
            else
                p2 = [NaN NaN NaN];
                R2 = NaN;
            end
            a2 = p2(1); a1 = p2(2); a0 = p2(3);
            if isfinite(a2) && abs(a2) > eps
                nuCurve = abs(-bodyForceX/(2*a2));
            else
                nuCurve = NaN;
            end
            shearBot = a1;
            shearTop = 2*a2*params.Ly + a1;
            tauBot = poiseuilleDiag(iseg,9);
            tauTop = poiseuilleDiag(iseg,10);
            rho2D = size(x,1)/(params.Lx*params.Ly);
            nuBotTau = abs(tauBot) / max(rho2D * abs(shearBot), eps);
            nuTopTau = abs(tauTop) / max(rho2D * abs(shearTop), eps);
            nuTauMean = mean([nuBotTau, nuTopTau], 'omitnan');
            qFlow = poiseuilleDiag(iseg,11);
            uFitAll = polyval(p2, yCenters);
            umaxFit = max(uFitAll);
            uCenterFit = polyval(p2, params.Ly/2);
            viscoDiag(iseg,:) = [a2, a1, a0, nuCurve, nuTauMean, tauBot, tauTop, rho2D, R2, qFlow, umaxFit, uCenterFit];
        end
    end

    if makePlots && (it == itStart || mod(it,nplot)==0 || it == nSteps)
        [Ux, Uy, Cnt, ~] = velocity_field_and_vorticity(x, v, params);
        plot_runtime_panels_generic(x, v, type, params, segIter(1:iseg), Ashape(1:iseg), ...
            traceDiag(1:iseg,:), densityWaveDiag(1:iseg,:), pressureDiag(1:iseg,:), poiseuilleDiag(1:iseg,:), viscoDiag(1:iseg,:), Ux, Uy, Cnt, it);
        drawnow;
    end
end

if lastCompletedStep == 0
    lastCompletedStep = min(nSteps, max(itStart-1, 0));
end
nDone = max(0, lastCompletedStep - itStart + 1);
segIter = segIter(1:nDone);
timeSeg = timeSeg(1:nDone);
E = E(1:nDone); occStd = occStd(1:nDone); outBand = outBand(1:nDone);
redistDiag = redistDiag(1:nDone,:); redistribStats = redistribStats(1:nDone,:); wallDiag = wallDiag(1:nDone,:);
msdAll = msdAll(1:nDone); msdType0 = msdType0(1:nDone); msdType1 = msdType1(1:nDone); nType0 = nType0(1:nDone); nType1 = nType1(1:nDone);
rRMS = rRMS(1:nDone); rMax = rMax(1:nDone); xCM = xCM(1:nDone); yCM = yCM(1:nDone); RxCtrl = RxCtrl(1:nDone); RyCtrl = RyCtrl(1:nDone); Ashape = Ashape(1:nDone);
traceDiag = traceDiag(1:nDone,:); radialDiag = radialDiag(1:nDone,:);
poiseuilleDiag = poiseuilleDiag(1:nDone,:); viscoDiag = viscoDiag(1:nDone,:);
densityWaveDiag = densityWaveDiag(1:nDone,:);
pressureDiag = pressureDiag(1:nDone,:);
pistonDiag = pistonDiag(1:nDone,:);
nSeg = nDone;

out = struct();
out.params = params; out.x = x; out.v = v; out.type = type; out.r0 = r0;
out.time = timeSeg; out.E = E; out.occStd = occStd; out.outBand = outBand;
out.redistDiag = redistDiag; out.redistribStats = redistribStats; out.wallDiag = wallDiag; out.pressureDiag = pressureDiag; out.pistonDiag = pistonDiag;
out.msdAll = msdAll; out.msdType0 = msdType0; out.msdType1 = msdType1; out.nType0 = nType0; out.nType1 = nType1;
out.lastCompletedStep = lastCompletedStep;
out.autoStop = autoStop;
out.lastLiquidClosure = lastLiquidClosure;
if isfield(params,'lastLiquidClosureDump')
    out.lastLiquidClosureDump = params.lastLiquidClosureDump;
end
if accumulateMeanFields
    out.meanFields = finalize_mean_fields(meanFieldSum, meanFieldNSamples, meanFieldGeom, params, meanFieldsStartStep, meanFieldsStride);
end
out.series = struct('E',E,'occStd',occStd,'outBand',outBand,'redistDiag',redistDiag,'redistribStats',redistribStats, ...
    'wallDiag',wallDiag,'pressureDiag',pressureDiag,'pistonDiag',pistonDiag,'msdAll',msdAll,'msdType0',msdType0,'msdType1',msdType1, ...
    'traceDiag',traceDiag,'radialDiag',radialDiag,'poiseuilleDiag',poiseuilleDiag,'viscoDiag',viscoDiag,'densityWaveDiag',densityWaveDiag, ...
    'Ashape',Ashape,'RxCtrl',RxCtrl,'RyCtrl',RyCtrl);

out.common = build_common_runtime_summary(E, occStd, outBand, pressureDiag, wallDiag, msdAll, msdType0, msdType1);

if computeShapeDiagnostics
    Rx_f = movmean_omitnan(RxCtrl, movAvgWindow);
    Ry_f = movmean_omitnan(RyCtrl, movAvgWindow);
    A_f  = movmean_omitnan(Ashape, movAvgWindow);
    nTail = max(20, round(0.10*nSeg));
    RxEq = mean(Rx_f(max(1,nSeg-nTail+1):nSeg), 'omitnan');
    RyEq = mean(Ry_f(max(1,nSeg-nTail+1):nSeg), 'omitnan');
    Aeq  = mean(A_f(max(1,nSeg-nTail+1):nSeg), 'omitnan');
    iFit0 = max(2, round(fitWindowFrac*nSeg));
    [tauRx, fitRxOk] = exp_relaxation_tau(timeSeg(iFit0:end), abs(Rx_f(iFit0:end) - RxEq));
    [tauRy, fitRyOk] = exp_relaxation_tau(timeSeg(iFit0:end), abs(Ry_f(iFit0:end) - RyEq));
    [tauA,  fitAOk ] = exp_relaxation_tau(timeSeg(iFit0:end), abs(A_f(iFit0:end)  - Aeq ));
    Rref = 0.5*(RxEq + RyEq);
    sigmaProxyRx = NaN; sigmaProxyRy = NaN; sigmaProxyA = NaN;
    if isfield(params,'nuEffForSigma') && isfield(params,'rhoEff') && isfinite(params.nuEffForSigma) && isfinite(params.rhoEff)
        if fitRxOk && tauRx > 0, sigmaProxyRx = params.rhoEff * params.nuEffForSigma * Rref / tauRx; end
        if fitRyOk && tauRy > 0, sigmaProxyRy = params.rhoEff * params.nuEffForSigma * Rref / tauRy; end
        if fitAOk  && tauA  > 0, sigmaProxyA  = params.rhoEff * params.nuEffForSigma * Rref / tauA; end
    end
    out.xCM = xCM; out.yCM = yCM; out.rRMS = rRMS; out.rMax = rMax;
    out.RxCtrl = RxCtrl; out.RyCtrl = RyCtrl; out.Ashape = Ashape;
    out.shape = struct('RxEq',RxEq,'RyEq',RyEq,'Aeq',Aeq,'tauRx',tauRx,'tauRy',tauRy,'tauA',tauA,...
        'sigmaProxyRx',sigmaProxyRx,'sigmaProxyRy',sigmaProxyRy,'sigmaProxyA',sigmaProxyA,...
        'fitOkRx',fitRxOk,'fitOkRy',fitRyOk,'fitOkA',fitAOk);
    out.sigma = out.shape;
    if printSummary
        fprintf('\n=== Resume shape diagnostics ===\n');
        fprintf('E finale                   = %.6g\n', last_or_nan(E));
        fprintf('A final / eq               = %.6g / %.6g\n', last_or_nan(Ashape), Aeq);
        fprintf('std(N)/gamma final         = %.6g\n', last_or_nan(occStd));
        fprintf('tauRx / tauRy / tauA       = %.6g / %.6g / %.6g\n', tauRx, tauRy, tauA);
        fprintf('sigmaProxyRx / Ry / A      = %.6g / %.6g / %.6g\n', sigmaProxyRx, sigmaProxyRy, sigmaProxyA);
        fprintf('===============================\n');
    end
end

if computeTracerDiagnostics
    E_f = movmean_omitnan(E, movAvgWindow);
    rhoStd_f = movmean_omitnan(redistDiag(:,1), movAvgWindow);
    outBand_f = movmean_omitnan(redistDiag(:,2), movAvgWindow);
    PglobErr_f = movmean_omitnan(redistDiag(:,3), movAvgWindow);
    R2tr_f = movmean_omitnan(traceDiag(:,6), movAvgWindow);
    Dinst_f = movmean_omitnan(traceDiag(:,7), movAvgWindow);
    fracDisk_f = movmean_omitnan(traceDiag(:,8), movAvgWindow);
    iFit0 = max(2, round(fitWindowFrac*nSeg));
    fitMask = false(nSeg,1); fitMask(iFit0:end) = true; fitMask = fitMask & isfinite(traceDiag(:,6));
    if nnz(fitMask) >= 5
        pR2 = polyfit(timeSeg(fitMask), traceDiag(fitMask,6), 1);
        slopeR2 = pR2(1); Dglobal = 0.25*slopeR2;
    else
        pR2 = [NaN NaN]; Dglobal = NaN;
    end
    Dinst_mean = mean(traceDiag(fitMask,7), 'omitnan');
    fracDisk_mean = mean(traceDiag(fitMask,8), 'omitnan');
    [rBins, cRadial] = tracer_radial_profile(x, type, xcTracer, ycTracer, nRadialBins);
    out.diffusion = struct('diffDiag',traceDiag,'radialDiag',radialDiag,'Dglobal',Dglobal,'DinstMean',Dinst_mean,...
        'fracDiskMean',fracDisk_mean,'rBins',rBins,'cRadial',cRadial,'E_f',E_f,'rhoStd_f',rhoStd_f,'outBand_f',outBand_f,...
        'PglobErr_f',PglobErr_f,'R2tr_f',R2tr_f,'Dinst_f',Dinst_f,'fracDisk_f',fracDisk_f,'fitMask',fitMask,'pR2',pR2);
    if printSummary
        fprintf('\n=== Resume diffusion traceur ===\n');
        fprintf('N type 1 final             = %.6g\n', last_or_nan(traceDiag(:,1)));
        fprintf('R2 traceur final           = %.6g\n', last_or_nan(traceDiag(:,6)));
        fprintf('D global                   = %.6g\n', Dglobal);
        fprintf('D instantane moyen         = %.6g\n', Dinst_mean);
        fprintf('frac type1 disque moyen    = %.6g\n', fracDisk_mean);
        fprintf('std(N)/gamma final         = %.6g\n', last_or_nan(redistDiag(:,1)));
        fprintf('=============================\n');
    end
end

if computeDensityWaveDiagnostics
    E_f = movmean_omitnan(E, movAvgWindow);
    rhoStd_f = movmean_omitnan(redistDiag(:,1), movAvgWindow);
    outBand_f = movmean_omitnan(redistDiag(:,2), movAvgWindow);
    Arel_f = movmean_omitnan(densityWaveDiag(:,2), movAvgWindow);
    Aabs_f = movmean_omitnan(densityWaveDiag(:,1), movAvgWindow);
    nTail = max(20, round(0.10*nSeg));
    ArelEq = mean(Arel_f(max(1,nSeg-nTail+1):nSeg), 'omitnan');
    iFit0 = max(2, round(fitWindowFrac*nSeg));
    [tauDensity, fitDensityOk] = exp_relaxation_tau(timeSeg(iFit0:end), abs(Arel_f(iFit0:end) - ArelEq));
    out.densitywave = struct('densityWaveDiag',densityWaveDiag,'E_f',E_f,'rhoStd_f',rhoStd_f,...
        'outBand_f',outBand_f,'Arel_f',Arel_f,'Aabs_f',Aabs_f,'ArelEq',ArelEq,...
        'tauDensity',tauDensity,'fitDensityOk',fitDensityOk);
    if printSummary
        fprintf('\n=== Resume onde de densite ===\n');
        fprintf('Amplitude relative finale  = %.6g\n', last_or_nan(densityWaveDiag(:,2)));
        fprintf('Amplitude relative eq      = %.6g\n', ArelEq);
        fprintf('tau densite                = %.6g\n', tauDensity);
        fprintf('std colonnes / gamma final = %.6g\n', last_or_nan(densityWaveDiag(:,4)));
        fprintf('std(N)/gamma final         = %.6g\n', last_or_nan(redistDiag(:,1)));
        fprintf('============================\n');
    end
end

if computePressureSummary
    E_f = movmean_omitnan(E, movAvgWindow);
    rhoStd_f = movmean_omitnan(redistDiag(:,1), movAvgWindow);
    outBand_f = movmean_omitnan(redistDiag(:,2), movAvgWindow);
    pMean_f = movmean_omitnan(pressureDiag(:,5), movAvgWindow);
    pBot_f = movmean_omitnan(pressureDiag(:,3), movAvgWindow);
    pTop_f = movmean_omitnan(pressureDiag(:,4), movAvgWindow);
    nTail = max(20, round(0.10*nSeg));
    idxTail = max(1,nSeg-nTail+1):nSeg;
    pBot_mean = mean(pressureDiag(idxTail,3), 'omitnan');
    pTop_mean = mean(pressureDiag(idxTail,4), 'omitnan');
    pMean_mean = mean(pressureDiag(idxTail,5), 'omitnan');
    pMean_std = std(pressureDiag(idxTail,5), 0, 'omitnan');
    if isfield(params,'useMovingPiston') && params.useMovingPiston
        rhoMean = size(x,1) / max(params.Lx * params.pistonYCurrent, eps);
    else
        rhoMean = size(x,1)/(params.Lx*params.Ly);
    end
    out.pressure = struct('pressureDiag',pressureDiag,'E_f',E_f,'rhoStd_f',rhoStd_f,...
        'outBand_f',outBand_f,'pMean_f',pMean_f,'pBot_f',pBot_f,'pTop_f',pTop_f,...
        'pBot_mean',pBot_mean,'pTop_mean',pTop_mean,'pMean_mean',pMean_mean,'pMean_std',pMean_std,'rhoMean',rhoMean);
    out.compressibility = out.pressure;
    if printSummary
        fprintf('\n=== Resume pression / compression ===\n');
        fprintf('rho moyenne                = %.6g\n', rhoMean);
        fprintf('pression bas / haut        = %.6g / %.6g\n', pBot_mean, pTop_mean);
        fprintf('pression moyenne           = %.6g\n', pMean_mean);
        fprintf('std pression (queue)       = %.6g\n', pMean_std);
        fprintf('std(N)/gamma final         = %.6g\n', last_or_nan(redistDiag(:,1)));
        fprintf('====================================\n');
    end
end

if computePoiseuilleDiagnostics
    E_f = movmean_omitnan(E, movAvgWindow);
    rhoStd_f = movmean_omitnan(redistDiag(:,1), movAvgWindow);
    outBand_f = movmean_omitnan(redistDiag(:,2), movAvgWindow);
    PglobErr_f = movmean_omitnan(redistDiag(:,3), movAvgWindow);
    nuCurve_f = movmean_omitnan(viscoDiag(:,4), movAvgWindow);
    nuTau_f = movmean_omitnan(viscoDiag(:,5), movAvgWindow);
    qFlow_f = movmean_omitnan(viscoDiag(:,10), movAvgWindow);
    R2_f = movmean_omitnan(viscoDiag(:,9), movAvgWindow);
    iSteady0 = max(1, round(steadyFracStart*nSeg));
    steadyMask = false(nSeg,1); steadyMask(iSteady0:end) = true;
    visMask = steadyMask & ~isnan(viscoDiag(:,4));
    nuCurve_mean = mean(viscoDiag(visMask,4), 'omitnan');
    nuTau_mean = mean(viscoDiag(visMask,5), 'omitnan');
    qFlow_mean = mean(viscoDiag(visMask,10), 'omitnan');
    R2_mean = mean(viscoDiag(visMask,9), 'omitnan');
    yCenters = ((0:Ny-1)+0.5)' * (params.Ly/Ny);
    [uxProf, ~, ~] = mean_profiles_y_fast(x, v, params.Lx, params.Ly, Nx, Ny);
    nuMean = nuCurve_mean;
    if isfinite(nuMean) && nuMean > 0
        uAna = 0.5*bodyForceX/nuMean * yCenters .* (params.Ly - yCenters);
    else
        uAna = nan(size(yCenters));
    end
    out.poiseuille = struct('poiseuilleDiag',poiseuilleDiag,'viscoDiag',viscoDiag,'nuCurve_f',nuCurve_f,'nuTau_f',nuTau_f,...
        'qFlow_f',qFlow_f,'R2_f',R2_f,'nuCurve_mean',nuCurve_mean,'nuTau_mean',nuTau_mean,'qFlow_mean',qFlow_mean,...
        'R2_mean',R2_mean,'uxProf',uxProf,'yCenters',yCenters,'uAna',uAna,'E_f',E_f,'rhoStd_f',rhoStd_f,...
        'outBand_f',outBand_f,'PglobErr_f',PglobErr_f);
    if printSummary
        fprintf('\n=== Resume Poiseuille ===\n');
        fprintf('ux moyen final             = %.6g\n', last_or_nan(poiseuilleDiag(:,3)));
        fprintf('slip bas / haut final      = %.6g / %.6g\n', last_or_nan(poiseuilleDiag(:,1)), last_or_nan(poiseuilleDiag(:,2)));
        fprintf('Q final                    = %.6g\n', last_or_nan(poiseuilleDiag(:,11)));
        fprintf('nuEff courbure moyen       = %.6g\n', nuCurve_mean);
        fprintf('nuEff paroi moyen          = %.6g\n', nuTau_mean);
        fprintf('Q moyen (regime etabli)    = %.6g\n', qFlow_mean);
        fprintf('R2 fit moyen               = %.6g\n', R2_mean);
        fprintf('==========================\n');
    end
end

if printSummary
    fprintf('\n=== Resume commun ===\n');
    fprintf('it final                   = %d\n', lastCompletedStep);
    fprintf('E finale                   = %.6g\n', out.common.Efinal);
    fprintf('std(N)/gamma final         = %.6g\n', out.common.occStdFinal);
    fprintf('outBand final              = %.6g\n', out.common.outBandFinal);
    fprintf('pMean final                = %.6g\n', out.common.pMeanFinal);
    fprintf('MSD totale finale          = %.6g\n', out.common.msdAllFinal);
    fprintf('=====================\n');
end

end

function params = prepare_simulator_params(params, x, type)
if nargin < 2, x = []; end
if nargin < 3 || isempty(type), type = []; end
if ~isfield(params,'n') || isempty(params.n)
    if ~isempty(x)
        params.n = size(x,1);
    else
        params.n = 0;
    end
end
if ~isfield(params,'makePlots') || isempty(params.makePlots), params.makePlots = true; end
if ~isfield(params,'nplot') || isempty(params.nplot)
    if isfield(params,'realtimeStride') && ~isempty(params.realtimeStride)
        params.nplot = params.realtimeStride;
    else
        params.nplot = 100;
    end
end
if ~isfield(params,'movAvgWindow') || isempty(params.movAvgWindow), params.movAvgWindow = 200; end
if ~isfield(params,'fitWindowFrac') || isempty(params.fitWindowFrac), params.fitWindowFrac = 0.5; end
if ~isfield(params,'steadyFracStart') || isempty(params.steadyFracStart), params.steadyFracStart = 0.6; end
if ~isfield(params,'fitMarginFrac') || isempty(params.fitMarginFrac), params.fitMarginFrac = 0.10; end
if ~isfield(params,'fitSmoothWindowY') || isempty(params.fitSmoothWindowY), params.fitSmoothWindowY = 7; end
if ~isfield(params,'nRadialBins') || isempty(params.nRadialBins), params.nRadialBins = 50; end
if ~isfield(params,'diagStride') || isempty(params.diagStride), params.diagStride = 20; end
if ~isfield(params,'viscoStride') || isempty(params.viscoStride), params.viscoStride = 20; end
if ~isfield(params,'xcTracer') || isempty(params.xcTracer), params.xcTracer = 0.5*params.Lx; end
if ~isfield(params,'ycTracer') || isempty(params.ycTracer), params.ycTracer = 0.5*params.Ly; end
if ~isfield(params,'rTracer') || isempty(params.rTracer), params.rTracer = 0.15*min(params.Lx, params.Ly); end
if ~isfield(params,'densityWaveAmp') || isempty(params.densityWaveAmp), params.densityWaveAmp = 0.10; end
if ~isfield(params,'densityWaveModeX') || isempty(params.densityWaveModeX), params.densityWaveModeX = 1; end
if ~isfield(params,'densityWaveDiagMode') || isempty(params.densityWaveDiagMode), params.densityWaveDiagMode = params.densityWaveModeX; end
if ~isfield(params,'computeShapeDiagnostics') || isempty(params.computeShapeDiagnostics), params.computeShapeDiagnostics = false; end
if ~isfield(params,'computeTracerDiagnostics') || isempty(params.computeTracerDiagnostics), params.computeTracerDiagnostics = false; end
if ~isfield(params,'computeDensityWaveDiagnostics') || isempty(params.computeDensityWaveDiagnostics), params.computeDensityWaveDiagnostics = false; end
if ~isfield(params,'computePistonDiagnostics') || isempty(params.computePistonDiagnostics)
    params.computePistonDiagnostics = ...
        (isfield(params,'useMovingPiston') && params.useMovingPiston) || ...
        (isfield(params,'boundary_top') && ~isempty(params.boundary_top) && strcmpi(char(string(params.boundary_top)),'piston')) || ...
        (isfield(params,'caseType') && ~isempty(params.caseType) && strcmpi(char(string(params.caseType)),'piston'));
end
if ~isfield(params,'computePressureSummary') || isempty(params.computePressureSummary)
    params.computePressureSummary = ...
        (isfield(params,'useMovingPiston') && params.useMovingPiston) || ...
        (isfield(params,'boundary_top') && ~isempty(params.boundary_top) && strcmpi(char(string(params.boundary_top)),'piston')) || ...
        (isfield(params,'caseType') && ~isempty(params.caseType) && strcmpi(char(string(params.caseType)),'compressibility'));
end
if ~isfield(params,'computePoiseuilleDiagnostics') || isempty(params.computePoiseuilleDiagnostics), params.computePoiseuilleDiagnostics = false; end
if ~isfield(params,'computeViscosityDiagnostics') || isempty(params.computeViscosityDiagnostics), params.computeViscosityDiagnostics = params.computePoiseuilleDiagnostics; end
if ~isfield(params,'printSummary') || isempty(params.printSummary), params.printSummary = true; end
if ~isfield(params,'quiverStrideX') || isempty(params.quiverStrideX), params.quiverStrideX = 8; end
if ~isfield(params,'quiverStrideY') || isempty(params.quiverStrideY), params.quiverStrideY = 8; end
if ~isfield(params,'quiverScale') || isempty(params.quiverScale), params.quiverScale = 2.0; end
if ~isfield(params,'dumpStride') || isempty(params.dumpStride), params.dumpStride = 8; end
if ~isfield(params,'stepFcnName') || isempty(params.stepFcnName), params.stepFcnName = 'mpcd_incompressible_step_liquidclosure'; end
if ~isfield(params,'nType0Init') || isempty(params.nType0Init)
    if ~isempty(type)
        params.nType0Init = sum(type == uint8(0));
        params.nType1Init = sum(type == uint8(1));
    else
        params.nType0Init = NaN;
        params.nType1Init = NaN;
    end
end
if ~isfield(params,'accumulateMeanFields') || isempty(params.accumulateMeanFields)
    if isfield(params,'useLiquidClosure') && ~isempty(params.useLiquidClosure)
        params.accumulateMeanFields = logical(params.useLiquidClosure);
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
end


function acc = init_mean_field_accumulator(params)
Nx = params.Nx; Ny = params.Ny;
acc = struct();
acc.Ux = zeros(Ny, Nx);
acc.Uy = zeros(Ny, Nx);
acc.rho = zeros(Ny, Nx);
acc.T = zeros(Ny, Nx);
acc.Pkin = zeros(Ny, Nx);
acc.Pvir = zeros(Ny, Nx);
acc.P = zeros(Ny, Nx);
acc.rhoTarget = zeros(Ny, Nx);
acc.nSamples = 0;
end

function acc = accumulate_mean_field_snapshot(acc, params, x, v)
snap = reconstruct_mean_field_snapshot(params, x, v);
acc.Ux = acc.Ux + snap.Ux;
acc.Uy = acc.Uy + snap.Uy;
acc.rho = acc.rho + snap.rho;
acc.T = acc.T + snap.T;
acc.Pkin = acc.Pkin + snap.Pkin;
acc.Pvir = acc.Pvir + snap.Pvir;
acc.P = acc.P + snap.P;
acc.rhoTarget = acc.rhoTarget + snap.rhoTarget;
acc.nSamples = acc.nSamples + 1;
end

function MF = finalize_mean_field_accumulator(acc, params)
MF = struct();
n = max(acc.nSamples, 1);
MF.Ux = acc.Ux / n;
MF.Uy = acc.Uy / n;
MF.rho = acc.rho / n;
MF.T = acc.T / n;
MF.Pkin = acc.Pkin / n;
MF.Pvir = acc.Pvir / n;
MF.P = acc.P / n;
MF.rhoTarget = acc.rhoTarget / n;
MF.nSamples = acc.nSamples;
MF.startStep = params.meanFieldsStartStep;
MF.stride = params.meanFieldsStride;
MF.dx = params.Lx / params.Nx;
MF.dy = params.Ly / params.Ny;
MF.xc = ((0:params.Nx-1)+0.5) * MF.dx;
MF.yc = ((0:params.Ny-1)+0.5) * MF.dy;
[MF.Xc, MF.Yc] = meshgrid(MF.xc, MF.yc);
end

function snap = reconstruct_mean_field_snapshot(params, x, v)
Nx = params.Nx; Ny = params.Ny;
dx = params.Lx / Nx; dy = params.Ly / Ny;
Vc = dx * dy;
ix = floor(x(:,1) / dx) + 1;
iy = floor(x(:,2) / dy) + 1;
ix = max(1, min(Nx, ix));
iy = max(1, min(Ny, iy));
ic = sub2ind([Ny, Nx], iy, ix);
Nc = Nx * Ny;

countVec = accumarray(ic, 1, [Nc 1], @sum, 0);
sumVx = accumarray(ic, v(:,1), [Nc 1], @sum, 0);
sumVy = accumarray(ic, v(:,2), [Nc 1], @sum, 0);

UxVec = zeros(Nc,1); UyVec = zeros(Nc,1);
nz = countVec > 0;
UxVec(nz) = sumVx(nz) ./ countVec(nz);
UyVec(nz) = sumVy(nz) ./ countVec(nz);

uxPart = UxVec(ic);
uyPart = UyVec(ic);
rel2 = (v(:,1) - uxPart).^2 + (v(:,2) - uyPart).^2;
sumRel2 = accumarray(ic, rel2, [Nc 1], @sum, 0);
dof = 2 * max(countVec - 1, 0);
Tvec = nan(Nc,1);
maskT = dof > 0;
Tvec(maskT) = sumRel2(maskT) ./ dof(maskT);

rhoVec = countVec ./ Vc;
PkinVec = nan(Nc,1);
PkinVec(maskT) = rhoVec(maskT) .* Tvec(maskT);

[~, rhoTarget] = target_maps_from_params_sim(params);
Kvirial = getv_local(params, 'Kvirial', NaN);
Pvir = Kvirial .* (reshape(rhoVec, [Ny, Nx]) - rhoTarget);
Pkin = reshape(PkinVec, [Ny, Nx]);

snap = struct();
snap.Ux = reshape(UxVec, [Ny, Nx]);
snap.Uy = reshape(UyVec, [Ny, Nx]);
snap.rho = reshape(rhoVec, [Ny, Nx]);
snap.T = reshape(Tvec, [Ny, Nx]);
snap.Pkin = Pkin;
snap.Pvir = Pvir;
snap.P = Pkin + Pvir;
snap.rhoTarget = rhoTarget;
end

function [targetOcc, rhoTarget] = target_maps_from_params_sim(params)
Nx = params.Nx; Ny = params.Ny;
aX = params.Lx / params.Nx;
aY = params.Ly / params.Ny;
Vc = aX * aY;
gammaEff = getv_local(params, 'gammaTargetCurrent', getv_local(params, 'gammaCurrent', getv_local(params,'gamma', [])));
activeFrac = active_cell_fraction_map_light_sim(params, [Ny, Nx]);
targetOcc = gammaEff .* activeFrac;
rhoTarget = targetOcc ./ Vc;
end

function fracMap = active_cell_fraction_map_light_sim(params, sz)
Ny = sz(1); Nx = sz(2);
if isfield(params,'activeCellFracCurrent') && ~isempty(params.activeCellFracCurrent)
    fracMap = params.activeCellFracCurrent;
    if isequal(size(fracMap), [Ny Nx]), return; end
end
fracMap = ones(Ny, Nx);
if isfield(params,'useMovingPiston') && params.useMovingPiston
    Ly = params.Ly;
    if isfield(params,'pistonYCurrent') && ~isempty(params.pistonYCurrent)
        yTop = params.pistonYCurrent;
    elseif isfield(params,'pistonY0') && ~isempty(params.pistonY0)
        yTop = params.pistonY0;
    else
        yTop = Ly;
    end
    margin = getv_local(params, 'pistonActiveMargin', 0.0);
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



function row = piston_runtime_diagnostics(params, x, v, stepInfo)
Nx = params.Nx; Ny = params.Ny;
fracMap = active_cell_fraction_map_light_sim(params, [Ny, Nx]);
activeFracMean = mean(fracMap(:), 'omitnan');
nActiveEff = sum(fracMap(:), 'omitnan');
activeVol = params.Lx * params.Ly * activeFracMean;
Hactive = activeVol / max(params.Lx, eps);
compression = Hactive / max(params.Ly, eps);
gammaGeom = getv_local(params,'gammaGeomCurrent', getv_local(params,'gammaCurrent', getv_local(params,'gamma', NaN)));
gammaTarget = getv_local(params,'gammaTargetCurrent', getv_local(params,'gamma', NaN));
rhoMeanActive = size(x,1) / max(activeVol, eps);
yTop = getv_local(params,'pistonYCurrent', getv_local(params,'pistonY0', params.Ly));
capacityRatio = gammaGeom / max(gammaTarget, eps);
pBot = NaN; pTop = NaN; pMean = NaN;
if nargin >= 3 && isstruct(stepInfo)
    pBot = getv_local(stepInfo,'pBotInst', NaN);
    pTop = getv_local(stepInfo,'pTopInst', NaN);
    pMean = getv_local(stepInfo,'pMeanInst', NaN);
end
[PkinActive, PvirActive, PtotActive, overfillActive] = piston_active_pressures(params, x, v, stepInfo, fracMap);
row = [yTop, Hactive, compression, activeVol, activeFracMean, nActiveEff, gammaGeom, gammaTarget, rhoMeanActive, capacityRatio, pBot, pTop, pMean, PkinActive, PvirActive, PtotActive, overfillActive];
end

function [PkinActive, PvirActive, PtotActive, overfillActive] = piston_active_pressures(params, x, v, stepInfo, fracMap)
PkinActive = NaN; PvirActive = NaN; PtotActive = NaN; overfillActive = NaN;
if nargin < 4 || isempty(fracMap)
    fracMap = active_cell_fraction_map_light_sim(params, [params.Ny, params.Nx]);
end
if nargin >= 3 && isstruct(stepInfo)
    if isfield(stepInfo,'Ncell') && ~isempty(stepInfo.Ncell)
        Ncell = double(stepInfo.Ncell);
        Tcell = [];
        if isfield(stepInfo,'Tcell') && ~isempty(stepInfo.Tcell)
            Tcell = double(stepInfo.Tcell);
        elseif isfield(stepInfo,'T') && ~isempty(stepInfo.T)
            Tcell = double(stepInfo.T);
        end
        if ~isempty(Tcell)
            dx = params.Lx / params.Nx;
            dy = params.Ly / params.Ny;
            Vc = dx*dy;
            rho = Ncell / Vc;
            Pkin = rho .* Tcell;
            [~, rhoTarget] = target_maps_from_params_sim(params);
            Kvirial = getv_local(params, 'Kvirial', NaN);
            Pvir = Kvirial .* (rho - rhoTarget);
            Ptot = Pkin + Pvir;
            w = double(fracMap);
            wsum = sum(w(:), 'omitnan');
            if wsum > 0
                PkinActive = sum(Pkin(:).*w(:), 'omitnan') / wsum;
                PvirActive = sum(Pvir(:).*w(:), 'omitnan') / wsum;
                PtotActive = sum(Ptot(:).*w(:), 'omitnan') / wsum;
                overfillActive = sum(max(0, rho(:) - rhoTarget(:)).*w(:), 'omitnan') / wsum;
                return;
            end
        end
    end
end
% Fallback robuste: reconstruction a partir de l'etat particulaire courant.
if nargin >= 3 && ~isempty(v)
    snap = reconstruct_mean_field_snapshot(params, x, v);
    w = double(fracMap);
    wsum = sum(w(:), 'omitnan');
    if wsum > 0
        PkinActive = sum(snap.Pkin(:).*w(:), 'omitnan') / wsum;
        PvirActive = sum(snap.Pvir(:).*w(:), 'omitnan') / wsum;
        PtotActive = sum(snap.P(:).*w(:), 'omitnan') / wsum;
        overfillActive = sum(max(0, snap.rho(:) - snap.rhoTarget(:)).*w(:), 'omitnan') / wsum;
    end
end
end

function [Beff, kEff, dPdrhoEff, BfromdPdrho] = piston_effective_compressibility(pistonDiag)
Beff = nan(size(pistonDiag,1),1);
kEff = nan(size(pistonDiag,1),1);
dPdrhoEff = nan(size(pistonDiag,1),1);
BfromdPdrho = nan(size(pistonDiag,1),1);
if isempty(pistonDiag) || size(pistonDiag,2) < 16
    return;
end
activeVol0 = pistonDiag(1,4);
Ptot0 = pistonDiag(1,16);
rho0 = pistonDiag(1,9);
if ~isfinite(activeVol0) || activeVol0 <= 0 || ~isfinite(Ptot0)
    return;
end
strain = 1 - pistonDiag(:,4) ./ activeVol0;   % = -Delta V / V0
dP = pistonDiag(:,16) - Ptot0;
mask = isfinite(strain) & isfinite(dP) & (strain > 0);
Beff(mask) = dP(mask) ./ strain(mask);
mask2 = mask & isfinite(Beff) & (abs(Beff) > eps);
kEff(mask2) = 1 ./ Beff(mask2);
if isfinite(rho0)
    drho = pistonDiag(:,9) - rho0;
    maskR = isfinite(drho) & isfinite(dP) & (abs(drho) > eps);
    dPdrhoEff(maskR) = dP(maskR) ./ drho(maskR);
    rhoUse = pistonDiag(:,9);
    maskB = maskR & isfinite(rhoUse);
    BfromdPdrho(maskB) = rhoUse(maskB) .* dPdrhoEff(maskB);
end
end

function val = getv_local(s, name, default)
if isfield(s,name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end

function summary = build_common_runtime_summary(E, occStd, outBand, pressureDiag, wallDiag, msdAll, msdType0, msdType1)
summary = struct();
summary.Efinal = last_or_nan(E);
summary.occStdFinal = last_or_nan(occStd);
summary.outBandFinal = last_or_nan(outBand);
summary.pMeanFinal = last_or_nan(pressureDiag(:,5));
summary.pBotFinal = last_or_nan(pressureDiag(:,3));
summary.pTopFinal = last_or_nan(pressureDiag(:,4));
summary.wallHitsBottomFinal = last_or_nan(wallDiag(:,1));
summary.wallHitsTopFinal = last_or_nan(wallDiag(:,2));
summary.msdAllFinal = last_or_nan(msdAll);
summary.msdType0Final = last_or_nan(msdType0);
summary.msdType1Final = last_or_nan(msdType1);
end

function val = last_or_nan(x)
val = NaN;
if isempty(x)
    return;
end
x = x(:);
idx = find(isfinite(x), 1, 'last');
if ~isempty(idx)
    val = x(idx);
end
end

function plot_runtime_panels_generic(x, v, type, params, segIter, Ashape, traceDiag, densityWaveDiag, pressureDiag, poiseuilleDiag, viscoDiag, Ux, Uy, Cnt, it)
Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
quiverStrideX = params.quiverStrideX;
quiverStrideY = params.quiverStrideY;
quiverScale = params.quiverScale;
dumpStride = params.dumpStride;
xc = ((0:Nx-1)+0.5) * (Lx/Nx);
yc = ((0:Ny-1)+0.5) * (Ly/Ny);
[Xc, Yc] = meshgrid(xc, yc);
qsx = 1:max(1,quiverStrideX):Nx;
qsy = 1:max(1,quiverStrideY):Ny;
idx = 1:max(1,dumpStride):size(x,1);

figure(1); clf;
subplot(2,2,1);
scatter(x(idx,1), x(idx,2), 4, double(type(idx)), 'filled');
xlim([0 Lx]); ylim([0 Ly]); axis equal tight; colorbar;
title(sprintf('Particules (it=%d)', it));

subplot(2,2,2);
imagesc(xc, yc, Cnt); set(gca,'YDir','normal'); axis equal tight; colorbar;
title(sprintf('Occupation | min=%d max=%d mean=%.2f', min(Cnt(:)), max(Cnt(:)), mean(Cnt(:))));

subplot(2,2,3);
imagesc(xc, yc, Ux); set(gca,'YDir','normal'); hold on;
quiver(Xc(qsy,qsx), Yc(qsy,qsx), Ux(qsy,qsx), Uy(qsy,qsx), quiverScale, 'k');
axis equal tight; colorbar; title('u_x + quiver');

subplot(2,2,4); hold on;
added = false;
if isfield(params,'computeShapeDiagnostics') && params.computeShapeDiagnostics && any(isfinite(Ashape))
    plot(segIter, Ashape, 'k-', 'DisplayName', 'Ashape');
    added = true;
end
if isfield(params,'computeTracerDiagnostics') && params.computeTracerDiagnostics
    if any(isfinite(traceDiag(:,6)))
        plot(segIter, traceDiag(:,6), '-', 'DisplayName', 'R2 traceur');
        added = true;
    end
    if any(isfinite(traceDiag(:,7)))
        plot(segIter, traceDiag(:,7), '--', 'DisplayName', 'Dinst');
        added = true;
    end
end
if isfield(params,'computeDensityWaveDiagnostics') && params.computeDensityWaveDiagnostics
    if any(isfinite(densityWaveDiag(:,2)))
        plot(segIter, densityWaveDiag(:,2), '-', 'DisplayName', 'A1/gamma');
        added = true;
    end
end
if isfield(params,'computePressureSummary') && params.computePressureSummary
    if any(isfinite(pressureDiag(:,5)))
        plot(segIter, pressureDiag(:,5), '-', 'DisplayName', 'pmean');
        added = true;
    end
    if any(isfinite(pressureDiag(:,7)))
        plot(segIter, pressureDiag(:,7), '--', 'DisplayName', 'std(N)/gamma');
        added = true;
    end
end
if isfield(params,'computePoiseuilleDiagnostics') && params.computePoiseuilleDiagnostics
    % Les diagnostics Poiseuille de base sont calcules a chaque pas :
    % ils doivent donc servir de source principale pour l'affichage.
    if any(isfinite(poiseuilleDiag(:,3)))
        plot(segIter, poiseuilleDiag(:,3), '-', 'DisplayName', 'ux moyen');
        added = true;
    end
    if size(poiseuilleDiag,2) >= 11 && any(isfinite(poiseuilleDiag(:,11)))
        plot(segIter, poiseuilleDiag(:,11), '--', 'DisplayName', 'Q');
        added = true;
    end
    if size(poiseuilleDiag,2) >= 1 && any(isfinite(poiseuilleDiag(:,1)))
        plot(segIter, poiseuilleDiag(:,1), ':', 'DisplayName', 'slip bas');
        added = true;
    end
    if size(poiseuilleDiag,2) >= 2 && any(isfinite(poiseuilleDiag(:,2)))
        plot(segIter, poiseuilleDiag(:,2), '-.', 'DisplayName', 'slip haut');
        added = true;
    end
    % Les estimations de viscosite peuvent rester NaN tant que le fit n'est
    % pas etabli : on ne les ajoute que lorsqu'elles existent reellement.
    if any(isfinite(viscoDiag(:,4)))
        plot(segIter, viscoDiag(:,4), 'o-', 'MarkerSize', 3, 'DisplayName', 'nu courbure');
        added = true;
    end
    if any(isfinite(viscoDiag(:,5)))
        plot(segIter, viscoDiag(:,5), 's--', 'MarkerSize', 3, 'DisplayName', 'nu paroi');
        added = true;
    end
end
if added
    grid on; legend('Location','best'); title('Diagnostics actifs');
else
    plot(segIter, zeros(size(segIter)), 'k-');
    grid on; title('Aucun diagnostic spécifique actif');
end
end


function y = movmean_omitnan(x, w)
x = x(:); n = numel(x);
y = nan(size(x));
hw = floor(w/2);
for i = 1:n
    i0 = max(1, i-hw);
    i1 = min(n, i+hw);
    y(i) = mean(x(i0:i1), 'omitnan');
end
end

function [tau, ok] = exp_relaxation_tau(t, y)
mask = isfinite(t) & isfinite(y) & (y > 0);
t = t(mask); y = y(mask);
ok = false; tau = NaN;
if numel(t) < 8
    return
end
y = movmean(y, max(3, floor(numel(y)/30)));
mask2 = y > max(1e-12, 0.05*max(y));
t = t(mask2); y = y(mask2);
if numel(t) < 8
    return
end
p = polyfit(t, log(y), 1);
if p(1) >= 0
    return
end
tau = -1/p(1);
ok = isfinite(tau) && tau > 0;
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


function diagRow = poiseuille_diagnostics_fast(x, v, Ly, Ny, wallInfo, dt, Lx)
aY = Ly/Ny;
botMask = x(:,2) < 1.5*aY;
topMask = x(:,2) > Ly - 1.5*aY;
slipBot = 0; slipTop = 0;
if any(botMask), slipBot = mean(v(botMask,1)); end
if any(topMask), slipTop = mean(v(topMask,1)); end
uxMean = mean(v(:,1));
uxRms = sqrt(mean(v(:,1).^2));
uyRms = sqrt(mean(v(:,2).^2));
tauBot = wallInfo.dPxBot / (Lx*dt);
tauTop = wallInfo.dPxTop / (Lx*dt);
qFlow = mean(v(:,1))*Ly;
diagRow = [slipBot slipTop uxMean uxRms uyRms wallInfo.dEwall wallInfo.nBot wallInfo.nTop tauBot tauTop qFlow];
end

function [uxProf, uyProf, rhoProf] = mean_profiles_y_fast(x, v, Lx, Ly, Nx, Ny)
aY = Ly/Ny;
iy = floor(x(:,2)/aY) + 1;
iy = min(max(iy,1),Ny);
countY = accumarray(iy, 1, [Ny,1], @sum, 0);
sumUxY = accumarray(iy, v(:,1), [Ny,1], @sum, 0);
sumUyY = accumarray(iy, v(:,2), [Ny,1], @sum, 0);
uxProf = zeros(Ny,1); uyProf = zeros(Ny,1);
nz = countY > 0;
uxProf(nz) = sumUxY(nz) ./ countY(nz);
uyProf(nz) = sumUyY(nz) ./ countY(nz);
rhoProf = countY / max(Nx,1);
end

function row = density_wave_diagnostics(Ncell, Nx, Ny, gamma, modeX)
Cnt = reshape(Ncell, [Ny, Nx]);
colMean = mean(Cnt, 1);
j = 0:(Nx-1);
fluc = colMean - mean(colMean);
coef = (2/numel(j)) * sum(fluc .* exp(-1i*2*pi*modeX*j/numel(j)));
ampAbs = abs(coef);
ampRel = ampAbs / max(gamma, eps);
stdColRel = std(colMean) / max(gamma, eps);
minColRel = min(colMean) / max(gamma, eps);
maxColRel = max(colMean) / max(gamma, eps);
rmsColRel = sqrt(mean((fluc/max(gamma,eps)).^2));
phase = angle(coef);
row = [ampAbs, ampRel, mean(colMean), stdColRel, minColRel, maxColRel, rmsColRel, phase];
end

function [row, radial] = tracer_diffusion_diagnostics(x, v, typ, xc0, yc0, r0, nBins, Lx, Ly, dt, iseg)
mask1 = (typ == 1);
x1 = x(mask1,1); y1 = x(mask1,2); N1 = numel(x1);
if N1 == 0
    row = nan(1,13); radial = nan(1,4); return;
end
xcm = mean(x1); ycm = mean(y1);
dx = x1 - xcm; dy = y1 - ycm;
R2x = mean(dx.^2); R2y = mean(dy.^2); R2tot = mean(dx.^2 + dy.^2);
persistent prevR2
if isempty(prevR2) || iseg == 1
    Dinst = NaN;
else
    Dinst = max(0, (R2tot - prevR2)/(4*dt));
end
prevR2 = R2tot;
rr0 = hypot(x1-xc0, y1-yc0);
fracInDisk = mean(rr0 <= r0);
[rBins, cRadial] = tracer_radial_profile(x, typ, xc0, yc0, nBins);
cPeak = max(cRadial); cCenter = cRadial(1);
uxMean = mean(v(:,1)); uxRms = sqrt(mean(v(:,1).^2)); uyRms = sqrt(mean(v(:,2).^2));
row = [N1 xcm ycm R2x R2y R2tot Dinst fracInDisk cPeak cCenter uxMean uxRms uyRms];
targetHalf = 0.5*cCenter;
iHalf = find(cRadial <= targetHalf, 1, 'first');
if isempty(iHalf) || iHalf == 1
    rHalf = NaN;
else
    rHalf = interp1(cRadial(iHalf-1:iHalf), rBins(iHalf-1:iHalf), targetHalf, 'linear', NaN);
end
target95 = 0.05*cCenter;
i95 = find(cRadial <= target95, 1, 'first');
if isempty(i95) || i95 == 1
    r95 = NaN;
else
    r95 = interp1(cRadial(i95-1:i95), rBins(i95-1:i95), target95, 'linear', NaN);
end
radial = [rHalf r95 cCenter cRadial(end)];
end

function [rBins, cRadial] = tracer_radial_profile(x, typ, xc0, yc0, nBins)
rr = hypot(x(:,1)-xc0, x(:,2)-yc0);
rMax = max(rr);
edges = linspace(0, rMax, nBins+1);
idx = discretize(rr, edges);
cRadial = nan(nBins,1);
rBins = 0.5*(edges(1:end-1) + edges(2:end));
for k = 1:nBins
    mask = (idx == k);
    if any(mask)
        cRadial(k) = mean(typ(mask) == 1);
    end
end
cRadial = fillmissing(cRadial, 'nearest');
rBins = rBins(:);
end
function G = compute_liquidclosure_mean_fields(x, v, params)
Nx = params.Nx;
Ny = params.Ny;
Lx = params.Lx;
Ly = params.Ly;
dx = Lx / Nx;
dy = Ly / Ny;
Vc = dx * dy;

ix = floor(x(:,1) / dx) + 1;
iy = floor(x(:,2) / dy) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
ic = sub2ind([Ny, Nx], iy, ix);

Nc = Nx * Ny;
countVec = accumarray(ic, 1, [Nc 1], @sum, 0);
sumVx = accumarray(ic, v(:,1), [Nc 1], @sum, 0);
sumVy = accumarray(ic, v(:,2), [Nc 1], @sum, 0);

UxVec = zeros(Nc,1);
UyVec = zeros(Nc,1);
maskN = (countVec > 0);
UxVec(maskN) = sumVx(maskN) ./ countVec(maskN);
UyVec(maskN) = sumVy(maskN) ./ countVec(maskN);

relx = v(:,1) - UxVec(ic);
rely = v(:,2) - UyVec(ic);
rel2 = relx.^2 + rely.^2;
Erel = 0.5 * accumarray(ic, rel2, [Nc 1], @sum, 0);

TVec = nan(Nc,1);
maskT = (countVec > 1);
TVec(maskT) = 2 * Erel(maskT) ./ (2 * max(countVec(maskT) - 1, 1));
TVec(~isfinite(TVec)) = 0;

rhoVec = countVec ./ Vc;
PkinVec = rhoVec .* TVec;
rhoTarget = closure_target_density_from_params(params);
PvirVec = zeros(Nc,1);
if isfield(params,'useLiquidClosure') && params.useLiquidClosure && isfield(params,'Kvirial') && ~isempty(params.Kvirial)
    PvirVec = params.Kvirial .* (rhoVec - rhoTarget(:));
end
PVec = PkinVec + PvirVec;

xc = ((1:Nx) - 0.5) * dx;
yc = ((1:Ny) - 0.5) * dy;
[Xc, Yc] = meshgrid(xc, yc);

G = struct();
G.dx = dx;
G.dy = dy;
G.xc = xc;
G.yc = yc;
G.Xc = Xc;
G.Yc = Yc;
G.Ux = reshape(UxVec, [Ny, Nx]);
G.Uy = reshape(UyVec, [Ny, Nx]);
G.rho = reshape(rhoVec, [Ny, Nx]);
G.T = reshape(TVec, [Ny, Nx]);
G.Pkin = reshape(PkinVec, [Ny, Nx]);
G.Pvir = reshape(PvirVec, [Ny, Nx]);
G.P = reshape(PVec, [Ny, Nx]);
G.rhoTarget = reshape(rhoTarget, [Ny, Nx]);
end

function rhoTarget = closure_target_density_from_params(params)
Nx = params.Nx;
Ny = params.Ny;
Lx = params.Lx;
Ly = params.Ly;
dx = Lx / Nx;
dy = Ly / Ny;
Vc = dx * dy;

if isfield(params,'gammaCurrent') && ~isempty(params.gammaCurrent) && isfinite(params.gammaCurrent)
    gammaEff = params.gammaCurrent;
elseif isfield(params,'gamma') && ~isempty(params.gamma) && isfinite(params.gamma)
    gammaEff = params.gamma;
else
    gammaEff = params.n / max(Nx*Ny,1);
end

activeFrac = ones(Ny, Nx);
if isfield(params,'activeCellFracCurrent') && ~isempty(params.activeCellFracCurrent)
    activeFrac = params.activeCellFracCurrent;
    if ~isequal(size(activeFrac), [Ny, Nx])
        activeFrac = ones(Ny, Nx);
    end
end

if isfield(params,'redistributionUseActiveVolumeTargets') && params.redistributionUseActiveVolumeTargets
    targetOcc = gammaEff .* activeFrac;
else
    targetOcc = gammaEff * ones(Ny, Nx);
end

rhoTarget = targetOcc ./ Vc;
rhoTarget = rhoTarget(:);
end

function M = finalize_mean_fields(sumFields, nSamples, geom, params, startStep, stride)
M = struct();
M.nSamples = nSamples;
M.startStep = startStep;
M.stride = stride;
M.params = params;
if ~isempty(fieldnames(geom))
    gnames = fieldnames(geom);
    for k = 1:numel(gnames)
        M.(gnames{k}) = geom.(gnames{k});
    end
end
if nSamples <= 0
    M.Ux = []; M.Uy = []; M.rho = []; M.T = []; M.P = [];
    M.Pkin = []; M.Pvir = []; M.rhoTarget = [];
    return;
end
M.Ux = sumFields.Ux / nSamples;
M.Uy = sumFields.Uy / nSamples;
M.rho = sumFields.rho / nSamples;
M.T = sumFields.T / nSamples;
M.P = sumFields.P / nSamples;
M.Pkin = sumFields.Pkin / nSamples;
M.Pvir = sumFields.Pvir / nSamples;
M.rhoTarget = sumFields.rhoTarget / nSamples;
end

function r = range_omitnan(x)
x = x(isfinite(x));
if isempty(x)
    r = NaN;
else
    r = max(x) - min(x);
end
end

function rs = rel_span(x, ref)
x = x(isfinite(x));
if isempty(x) || ~isfinite(ref)
    rs = NaN;
else
    rs = (max(x) - min(x)) / max(abs(ref), 1e-12);
end
end

