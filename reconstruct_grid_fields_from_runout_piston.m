function G = reconstruct_grid_fields_from_runout_piston(runOut, varargin)
%RECONSTRUCT_GRID_FIELDS_FROM_RUNOUT_PISTON
% Post-traitement dedie au cas piston.
%
% Exploite :
%   - l'etat final (ou l'etat moyen temporel si meanFields disponible),
%   - les diagnostics temporels piston enregistres par
%     mpcd_incompressible_simulator_pistondiag,
%   - les diagnostics communs pressureDiag/redistDiag.
%
% Figures produites :
%   1) series temporelles piston : position, hauteur active, compression,
%      gammaGeom, gammaTarget, rho active moyenne, rapport de capacite, pressions bas/haut/moyenne,
%      std(N)/gamma et outBand.
%   2) champs grille de l'etat selectionne (final ou moyen)
%   3) profils verticaux Ux, Uy, rho, T, Pkin et champ cible sous piston
%
% Usage :
%   G = reconstruct_grid_fields_from_runout_piston(runOut);
%   G = reconstruct_grid_fields_from_runout_piston(runOut,'state','mean');
%
% Options :
%   'state' : 'final' (defaut) | 'mean'
%   'plotFigures' : true/false
%
% Sorties :
%   G.final, G.meanFields, G.piston, G.params

p = inputParser;
p.addParameter('state','final',@(s)ischar(s)||isstring(s));
p.addParameter('plotFigures',true,@(x)islogical(x)||isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

G = struct();
G.params = resolve_params(runOut);
G.final = reconstruct_state_from_runout(runOut,'final');
G.meanFields = [];
if has_mean_fields(runOut)
    G.meanFields = reconstruct_state_from_runout(runOut,'mean');
end
G.selected = pick_selected_state(G, lower(string(opt.state)));
G.piston = collect_piston_data(runOut, G.params);

if opt.plotFigures
    if ~isempty(G.piston.time)
        plot_piston_time_series(G.piston);
    end
    plot_state_fields_piston(G.selected);
    plot_state_profiles_piston(G);
end

end

% -------------------------------------------------------------------------
function params = resolve_params(runOut)
params = [];
if isfield(runOut,'simOut') && isfield(runOut.simOut,'params') && ~isempty(runOut.simOut.params)
    params = runOut.simOut.params;
elseif isfield(runOut,'params') && ~isempty(runOut.params)
    params = runOut.params;
elseif isfield(runOut,'finalState') && isfield(runOut.finalState,'params')
    params = runOut.finalState.params;
end
assert(~isempty(params),'Impossible de retrouver params.');
end

function tf = has_mean_fields(runOut)
tf = isfield(runOut,'simOut') && isfield(runOut.simOut,'meanFields') && ~isempty(runOut.simOut.meanFields);
if ~tf
    tf = isfield(runOut,'meanFields') && ~isempty(runOut.meanFields);
end
end

function S = pick_selected_state(G, stateName)
switch stateName
    case "mean"
        if isempty(G.meanFields)
            warning('Aucun champ moyen temporel disponible, bascule vers l''etat final.');
            S = G.final;
        else
            S = G.meanFields;
        end
    otherwise
        S = G.final;
end
end

function Gs = reconstruct_state_from_runout(runOut, stateName)
params = resolve_params(runOut);
if strcmpi(stateName,'mean')
    if isfield(runOut,'simOut') && isfield(runOut.simOut,'meanFields') && ~isempty(runOut.simOut.meanFields)
        MF = runOut.simOut.meanFields;
    else
        MF = runOut.meanFields;
    end
    Gs = struct();
    Gs.params = params;
    Gs.isMean = true;
    Gs.dx = MF.dx; Gs.dy = MF.dy;
    Gs.xc = MF.xc; Gs.yc = MF.yc;
    Gs.Xc = MF.Xc; Gs.Yc = MF.Yc;
    Gs.Ux = MF.Ux; Gs.Uy = MF.Uy;
    Gs.rho = MF.rho; Gs.T = MF.T;
    Gs.Pkin = choose_field(MF,'Pkin',choose_field(MF,'P',MF.rho.*MF.T));
    Gs.Pvir = choose_field(MF,'Pvir',nan(size(Gs.Pkin)));
    Gs.Ptot = choose_field(MF,'P',Gs.Pkin+Gs.Pvir);
    Gs.rhoTarget = choose_field(MF,'rhoTarget',target_map_from_params(params));
    Gs.nSamples = choose_field(MF,'nSamples',NaN);
else
    if isfield(runOut,'finalState') && ~isempty(runOut.finalState)
        S = runOut.finalState;
    elseif isfield(runOut,'simOut')
        S = runOut.simOut;
    else
        S = runOut;
    end
    x = double(S.x); v = double(S.v);
    Nx = params.Nx; Ny = params.Ny; Lx = params.Lx; Ly = params.Ly;
    dx = Lx/Nx; dy = Ly/Ny; Vc = dx*dy;
    [ix,iy,ic] = cell_indices(x, params);
    Nc = Nx*Ny;
    countVec = accumarray(ic,1,[Nc 1],@sum,0);
    sumVx = accumarray(ic,v(:,1),[Nc 1],@sum,0);
    sumVy = accumarray(ic,v(:,2),[Nc 1],@sum,0);
    UxVec = zeros(Nc,1); UyVec = zeros(Nc,1);
    nz = countVec>0;
    UxVec(nz)=sumVx(nz)./countVec(nz);
    UyVec(nz)=sumVy(nz)./countVec(nz);
    uxPart = UxVec(ic); uyPart = UyVec(ic);
    rel2 = (v(:,1)-uxPart).^2 + (v(:,2)-uyPart).^2;
    sumRel2 = accumarray(ic,rel2,[Nc 1],@sum,0);
    dof = 2*max(countVec-1,0);
    Tvec = nan(Nc,1); maskT = dof>0;
    Tvec(maskT)=sumRel2(maskT)./dof(maskT);
    rhoVec = countVec./Vc;
    PkinVec = nan(Nc,1); PkinVec(maskT)=rhoVec(maskT).*Tvec(maskT);

    Gs = struct();
    Gs.params = params;
    Gs.isMean = false;
    Gs.dx = dx; Gs.dy = dy;
    Gs.xc = ((1:Nx)-0.5)*dx; Gs.yc = ((1:Ny)-0.5)*dy;
    [Gs.Xc,Gs.Yc] = meshgrid(Gs.xc,Gs.yc);
    Gs.Ux = reshape(UxVec,[Ny,Nx]); Gs.Uy = reshape(UyVec,[Ny,Nx]);
    Gs.rho = reshape(rhoVec,[Ny,Nx]); Gs.T = reshape(Tvec,[Ny,Nx]);
    Gs.Pkin = reshape(PkinVec,[Ny,Nx]);
    Gs.rhoTarget = target_map_from_params(params);
    Kvirial = getf(params,'Kvirial',NaN);
    Gs.Pvir = Kvirial*(Gs.rho - Gs.rhoTarget);
    Gs.Ptot = Gs.Pkin + Gs.Pvir;
    Gs.Ncell = reshape(countVec,[Ny,Nx]);
end
Gs.activeFrac = active_frac_from_params(params,[size(Gs.rho,1),size(Gs.rho,2)]);
Gs.profiles = build_profiles(Gs);
end

function prof = build_profiles(Gs)
prof = struct();
prof.y = Gs.yc(:);
prof.Ux = mean(Gs.Ux,2,'omitnan');
prof.Uy = mean(Gs.Uy,2,'omitnan');
prof.rho = mean(Gs.rho,2,'omitnan');
prof.T = mean(Gs.T,2,'omitnan');
prof.Pkin = mean(Gs.Pkin,2,'omitnan');
prof.Pvir = mean(Gs.Pvir,2,'omitnan');
prof.Ptot = mean(Gs.Ptot,2,'omitnan');
prof.rhoTarget = mean(Gs.rhoTarget,2,'omitnan');
prof.activeFrac = mean(Gs.activeFrac,2,'omitnan');
end

function P = collect_piston_data(runOut, params)
P = struct();
if isfield(runOut,'simOut') && isfield(runOut.simOut,'pistonDiag') && ~isempty(runOut.simOut.pistonDiag)
    pistonDiag = runOut.simOut.pistonDiag;
elseif isfield(runOut,'pistonDiag') && ~isempty(runOut.pistonDiag)
    pistonDiag = runOut.pistonDiag;
elseif isfield(runOut,'simOut') && isfield(runOut.simOut,'piston') && isfield(runOut.simOut.piston,'pistonDiag')
    pistonDiag = runOut.simOut.piston.pistonDiag;
else
    pistonDiag = [];
end

if isempty(pistonDiag)
    P.time = [];
    return;
end

dt = params.dt;
n = size(pistonDiag,1);
P.time = ((1:n)-1)'*dt;
P.yTop = pistonDiag(:,1);
P.Hactive = pistonDiag(:,2);
P.compression = pistonDiag(:,3);
P.activeVol = pistonDiag(:,4);
P.activeFracMean = pistonDiag(:,5);
P.nActiveEff = pistonDiag(:,6);
P.gammaGeom = pistonDiag(:,7);
if size(pistonDiag,2) >= 8
    P.gammaTarget = pistonDiag(:,8);
else
    P.gammaTarget = getf(params,'gamma',nan(size(P.time)));
end
if size(pistonDiag,2) >= 9
    P.rhoMeanActive = pistonDiag(:,9);
    P.capacityRatio = pistonDiag(:,10);
    P.pBot = pistonDiag(:,11);
    P.pTop = pistonDiag(:,12);
    P.pMean = pistonDiag(:,13);
else
    P.rhoMeanActive = pistonDiag(:,8);
    P.capacityRatio = P.gammaGeom ./ max(P.gammaTarget, eps);
    P.pBot = pistonDiag(:,9);
    P.pTop = pistonDiag(:,10);
    P.pMean = pistonDiag(:,11);
end
if size(pistonDiag,2) >= 17
    P.PkinActive = pistonDiag(:,14);
    P.PvirActive = pistonDiag(:,15);
    P.PtotActive = pistonDiag(:,16);
    P.overfillActive = pistonDiag(:,17);
elseif size(pistonDiag,2) >= 14
    P.PkinActive = pistonDiag(:,12);
    P.PvirActive = pistonDiag(:,13);
    P.PtotActive = pistonDiag(:,14);
    P.overfillActive = nan(size(P.time));
else
    P.PkinActive = nan(size(P.time));
    P.PvirActive = nan(size(P.time));
    P.PtotActive = nan(size(P.time));
    P.overfillActive = nan(size(P.time));
end

if isfield(runOut,'simOut') && isfield(runOut.simOut,'pressureDiag') && ~isempty(runOut.simOut.pressureDiag)
    P.pressureDiag = runOut.simOut.pressureDiag;
end
if isfield(runOut,'simOut') && isfield(runOut.simOut,'redistDiag') && ~isempty(runOut.simOut.redistDiag)
    P.redistDiag = runOut.simOut.redistDiag;
else
    P.redistDiag = [];
end
[P.Beff, P.kEff, P.dPdrhoEff, P.BfromdPdrho] = piston_effective_compressibility_from_series(P.activeVol, P.PtotActive, P.rhoMeanActive);
P.gammaInit = getf(params,'gamma',NaN);
P.Ly = params.Ly;
P.Lx = params.Lx;
end

function plot_piston_time_series(P)
figure('Name','piston - diagnostics temporels','Color','w');
tiledlayout(4,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(P.time, P.yTop,'LineWidth',1.5,'DisplayName','y_{piston}');
plot(P.time, P.Hactive,'--','LineWidth',1.3,'DisplayName','H_{actif}');
grid on; xlabel('t'); ylabel('Longueur'); title('Position piston / hauteur active'); legend('Location','best');

nexttile; hold on;
plot(P.time, P.compression,'LineWidth',1.5,'DisplayName','H_{actif}/L_y');
yline(1,'k:');
grid on; xlabel('t'); ylabel('Compression'); title('Compression geometrique');

nexttile; hold on;
plot(P.time, P.gammaGeom,'LineWidth',1.5,'DisplayName','\gamma_{geom}');
plot(P.time, P.gammaTarget,'--','LineWidth',1.3,'DisplayName','\gamma_{target}');
if isfinite(P.gammaInit), yline(P.gammaInit,'k:','DisplayName','\gamma_{init}'); end
grid on; xlabel('t'); ylabel('\gamma'); title('Occupation geometrique / cible'); legend('Location','best');

nexttile; hold on;
plot(P.time, P.rhoMeanActive,'LineWidth',1.5,'DisplayName','\rho_{moy,actif}');
if any(isfinite(P.capacityRatio))
    yyaxis right;
    plot(P.time, P.capacityRatio,'--','LineWidth',1.2,'DisplayName','\gamma_{geom}/\gamma_{target}');
    ylabel('rapport');
    yyaxis left;
end
grid on; xlabel('t'); ylabel('Densite'); title('Densite active / sur-occupation requise');

nexttile; hold on;
plot(P.time, P.pBot,'LineWidth',1.5,'DisplayName','p_{bot}');
plot(P.time, P.pTop,'--','LineWidth',1.3,'DisplayName','p_{top}');
plot(P.time, P.pMean,'-.','LineWidth',1.3,'DisplayName','p_{mean}');
grid on; xlabel('t'); ylabel('Pression'); title('Pressions piston/paroi'); legend('Location','best');

nexttile; hold on;
plot(P.time, P.PkinActive,'LineWidth',1.5,'DisplayName','P_{kin,actif}');
plot(P.time, P.PvirActive,'--','LineWidth',1.3,'DisplayName','P_{vir,actif}');
plot(P.time, P.PtotActive,'-.','LineWidth',1.3,'DisplayName','P_{tot,actif}');
grid on; xlabel('t'); ylabel('Pression'); title('Pressions actives du modele'); legend('Location','best');

nexttile; hold on;
if ~isempty(P.redistDiag)
    plot(P.time, P.redistDiag(:,1),'LineWidth',1.5,'DisplayName','std(N)/\gamma');
    plot(P.time, P.redistDiag(:,2),'--','LineWidth',1.3,'DisplayName','outBand');
    legend('Location','best');
end
grid on; xlabel('t'); ylabel('Diagnostic'); title('Redistribution');

nexttile; hold on;
plot(P.time, P.Beff,'LineWidth',1.5,'DisplayName','B_{eff}');
if any(isfinite(P.kEff))
    yyaxis right;
    plot(P.time, P.kEff,'--','LineWidth',1.3,'DisplayName','\kappa_{eff}');
    ylabel('\kappa_{eff}');
    yyaxis left;
end
grid on; xlabel('t'); ylabel('B_{eff}'); title('Compressibilite effective');

nexttile; hold on;

mask1 = isfinite(P.dPdrhoEff);
mask2 = isfinite(P.BfromdPdrho);

yyaxis left
if any(mask1)
    plot(P.time(mask1), P.dPdrhoEff(mask1), 'LineWidth', 1.5, ...
        'DisplayName', 'dP_{tot}/d\rho eff');
end
ylabel('dP_{tot}/d\rho');

yyaxis right
if any(mask2)
    plot(P.time(mask2), P.BfromdPdrho(mask2), '--', 'LineWidth', 1.3, ...
        'DisplayName', '\rho\, dP_{tot}/d\rho');
end
ylabel('\rho\, dP_{tot}/d\rho');

grid on;
xlabel('t');
title('Raideur d''EOS effective');
legend('Location','best');

nexttile; axis off;
text(0.0,0.88,sprintf('\gamma_{init}=%.4g', P.gammaInit),'Interpreter','tex');
text(0.0,0.72,sprintf('\gamma_{geom}(t_f)=%.4g', last_or_nan(P.gammaGeom)),'Interpreter','tex');
text(0.0,0.56,sprintf('\gamma_{target}(t_f)=%.4g', last_or_nan(P.gammaTarget)),'Interpreter','tex');
text(0.0,0.40,sprintf('P_{kin,actif}(t_f)=%.4g', last_or_nan(P.PkinActive)),'Interpreter','tex');
text(0.0,0.24,sprintf('P_{vir,actif}(t_f)=%.4g', last_or_nan(P.PvirActive)),'Interpreter','tex');
text(0.0,0.08,sprintf('P_{tot,actif}(t_f)=%.4g', last_or_nan(P.PtotActive)),'Interpreter','tex');
text(0.0,-0.08,sprintf('B_{eff}(t_f)=%.4g', last_or_nan(P.Beff)),'Interpreter','tex');
text(0.0,-0.24,sprintf('\kappa_{eff}(t_f)=%.4g', last_or_nan(P.kEff)),'Interpreter','tex');
text(0.0,-0.40,sprintf('dP_{tot}/d\rho (t_f)=%.4g', last_or_nan(P.dPdrhoEff)),'Interpreter','tex');
text(0.0,-0.56,sprintf('\rho dP_{tot}/d\rho (t_f)=%.4g', last_or_nan(P.BfromdPdrho)),'Interpreter','tex');
set(gca,'XColor','none','YColor','none');
end

function plot_state_fields_piston(Gs)
figure('Name','piston - champs etat selectionne','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(Gs.xc, Gs.yc, Gs.rho); set(gca,'YDir','normal'); axis equal tight; colorbar; title('\rho');
nexttile; imagesc(Gs.xc, Gs.yc, Gs.Ux); set(gca,'YDir','normal'); axis equal tight; colorbar; title('U_x');
nexttile; imagesc(Gs.xc, Gs.yc, Gs.Uy); set(gca,'YDir','normal'); axis equal tight; colorbar; title('U_y');
nexttile; imagesc(Gs.xc, Gs.yc, Gs.T); set(gca,'YDir','normal'); axis equal tight; colorbar; title('kBT');
nexttile; imagesc(Gs.xc, Gs.yc, Gs.Pkin); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{kin}');
nexttile; imagesc(Gs.xc, Gs.yc, Gs.activeFrac); set(gca,'YDir','normal'); axis equal tight; colorbar; title('fraction active');
end

function plot_state_profiles_piston(G)
Gs = G.selected;
y = Gs.profiles.y;

figure('Name','piston - profils verticaux','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(y, Gs.profiles.Ux,'LineWidth',1.6,'DisplayName',label_state(Gs));
if ~isempty(G.meanFields)
    plot(y, G.meanFields.profiles.Ux,'--','LineWidth',1.3,'DisplayName','moyen temporel');
end
grid on; xlabel('y'); ylabel('<U_x>_x'); title('Profil U_x(y)'); legend('Location','best');

nexttile; hold on;
plot(y, Gs.profiles.Uy,'LineWidth',1.6,'DisplayName',label_state(Gs));
if ~isempty(G.meanFields)
    plot(y, G.meanFields.profiles.Uy,'--','LineWidth',1.3,'DisplayName','moyen temporel');
end
grid on; xlabel('y'); ylabel('<U_y>_x'); title('Profil U_y(y)'); legend('Location','best');

nexttile; hold on;
plot(y, Gs.profiles.rho,'LineWidth',1.6,'DisplayName',label_state(Gs));
plot(y, Gs.profiles.rhoTarget,'k--','LineWidth',1.3,'DisplayName','\rho_{target}');
if ~isempty(G.meanFields)
    plot(y, G.meanFields.profiles.rho,':','LineWidth',1.3,'DisplayName','\rho moyen temp.');
end
grid on; xlabel('y'); ylabel('<\rho>_x'); title('Densite'); legend('Location','best');

nexttile; hold on;
plot(y, Gs.profiles.T,'LineWidth',1.6,'DisplayName',label_state(Gs));
if ~isempty(G.meanFields)
    plot(y, G.meanFields.profiles.T,'--','LineWidth',1.3,'DisplayName','moyen temporel');
end
grid on; xlabel('y'); ylabel('<kBT>_x'); title('Temperature'); legend('Location','best');

nexttile; hold on;
plot(y, Gs.profiles.Pkin,'LineWidth',1.6,'DisplayName','P_{kin}');
plot(y, Gs.profiles.Pvir,'--','LineWidth',1.3,'DisplayName','P_{vir}');
plot(y, Gs.profiles.Ptot,'-.','LineWidth',1.3,'DisplayName','P_{tot}');
grid on; xlabel('y'); ylabel('Pression'); title('Pressions'); legend('Location','best');

nexttile; hold on;
plot(y, Gs.profiles.activeFrac,'LineWidth',1.6,'DisplayName','fraction active');
if ~isempty(G.meanFields)
    plot(y, G.meanFields.profiles.activeFrac,'--','LineWidth',1.3,'DisplayName','fraction active moyenne');
end
grid on; xlabel('y'); ylabel('f_{actif}'); title('Profil de domaine actif'); legend('Location','best');
end

function s = label_state(Gs)
if isfield(Gs,'isMean') && Gs.isMean
    s = 'moyen temporel';
else
    s = 'etat final';
end
end

function [Beff, kEff, dPdrhoEff, BfromdPdrho] = ...
    piston_effective_compressibility_from_series(activeVol, PtotActive, rhoMeanActive)

Beff = nan(size(PtotActive));
kEff = nan(size(PtotActive));
dPdrhoEff = nan(size(PtotActive));
BfromdPdrho = nan(size(PtotActive));

if isempty(activeVol) || isempty(PtotActive)
    return;
end

V0 = activeVol(1);
P0 = PtotActive(1);

if ~isfinite(V0) || V0 <= 0 || ~isfinite(P0)
    return;
end

% Compression volumique relative : -Delta V / V0
strain = 1 - activeVol ./ V0;

% Variation de pression totale active
dP = PtotActive - P0;

maskB = isfinite(strain) & isfinite(dP) & (strain > 0);
Beff(maskB) = dP(maskB) ./ strain(maskB);

maskK = maskB & isfinite(Beff) & (abs(Beff) > eps);
kEff(maskK) = 1 ./ Beff(maskK);

if nargin < 3 || isempty(rhoMeanActive)
    return;
end

rho0 = rhoMeanActive(1);
if ~isfinite(rho0)
    return;
end

drho = rhoMeanActive - rho0;

maskR = isfinite(drho) & isfinite(dP) & (abs(drho) > eps);
dPdrhoEff(maskR) = dP(maskR) ./ drho(maskR);

maskB2 = maskR & isfinite(rhoMeanActive);
BfromdPdrho(maskB2) = rhoMeanActive(maskB2) .* dPdrhoEff(maskB2);
end

function val = choose_field(S, name, defaultVal)
if isfield(S,name) && ~isempty(S.(name))
    val = S.(name);
else
    val = defaultVal;
end
end

function rhoTarget = target_map_from_params(params)
Nx = params.Nx; Ny = params.Ny;
dx = params.Lx / Nx; dy = params.Ly / Ny;
Vc = dx*dy;
gammaEff = getf(params,'gammaTargetCurrent',getf(params,'gammaCurrent',getf(params,'gamma',NaN)));
frac = active_frac_from_params(params,[Ny,Nx]);
rhoTarget = gammaEff * frac / Vc;
end

function frac = active_frac_from_params(params, sz)
Ny = sz(1); Nx = sz(2);
if isfield(params,'activeCellFracCurrent') && ~isempty(params.activeCellFracCurrent) && isequal(size(params.activeCellFracCurrent),[Ny,Nx])
    frac = params.activeCellFracCurrent;
    return;
end
frac = ones(Ny,Nx);
if isfield(params,'useMovingPiston') && params.useMovingPiston
    Ly = params.Ly;
    yTop = getf(params,'pistonYCurrent', getf(params,'pistonY0', Ly));
    margin = getf(params,'pistonActiveMargin',0);
    yTopEff = min(max(yTop - margin, 0), Ly);
    dy = Ly/Ny;
    fracRows = zeros(Ny,1);
    for iy = 1:Ny
        y0 = (iy-1)*dy; y1 = iy*dy;
        fracRows(iy) = max(0, min(yTopEff,y1)-y0) / dy;
    end
    frac = repmat(fracRows,1,Nx);
end
end

function [ix,iy,ic] = cell_indices(x, params)
Nx = params.Nx; Ny = params.Ny;
dx = params.Lx/Nx; dy = params.Ly/Ny;
ix = floor(x(:,1)/dx)+1;
iy = floor(x(:,2)/dy)+1;
ix = max(1,min(Nx,ix));
iy = max(1,min(Ny,iy));
ic = sub2ind([Ny,Nx],iy,ix);
end

function v = getf(s, name, defaultVal)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultVal;
end
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