function G = reconstruct_grid_fields_from_runout_liquidclosure_Rh(runOut, varargin)
%RECONSTRUCT_GRID_FIELDS_FROM_RUNOUT_LIQUIDCLOSURE_RH
% Reconstruit des champs grille a partir d'une structure runOut (ou simOut),
% analyse les dumps de liquid closure si disponibles, et exploite aussi les
% champs moyens temporels stockes par le simulateur (simOut.meanFields).
%
% Sorties principales pour l'etat selectionne :
%   G.Ncell   : occupation par cellule
%   G.rho     : densite surfacique = Ncell / (dx*dy)
%   G.Ux, Uy  : vitesse moyenne par cellule
%   G.T       : temperature cinetique locale (ici kBT local)
%   G.P       : pression cinetique locale = rho .* T  (= Pkin)
%   G.Pkin    : alias explicite de G.P
%   G.Pvir    : pression virielle d'etat = Kvirial * (rho - rhoTarget)
%   G.Ptot    : pression d'etat = Pkin + Pvir
%   G.profiles.Rh_target : residu hydrostatique dP/dy - rho_target*g
%   G.profiles.Rh_state  : residu hydrostatique dP/dy - rho*g
%
% Si des dumps de liquid closure sont presents dans runOut/simOut, le champ
% G.liquidClosure contient les memes grandeurs pour chaque etape disponible.
%
% Si le simulateur a accumule des champs moyens temporels dans simOut.meanFields,
% le champ G.meanFields expose :
%   - Ux, Uy, rho, T, Pkin, Pvir, Ptot moyens temporels
%   - les profils verticaux moyens associes
%   - le residu hydrostatique reconstruit a partir de ces champs moyens.
%
% Usage:
%   G = reconstruct_grid_fields_from_runout_liquidclosure_Rh(runOut);
%
% Options:
%   'state'               : 'final' (defaut) | 'simOut' | 'initial' |
%                           'mean' | 'meanfields' |
%                           'liquid_reference' | 'liquid_redistributed' |
%                           'liquid_repaired' | 'liquid_postkick' |
%                           'liquid_prereorientation' |
%                           'liquid_postreorientation'
%   'useType'             : true/false (defaut false)
%   'plotFigures'         : true/false (defaut true)
%   'plotLiquidClosure'   : true/false (defaut true)
%   'plotMeanFields'      : true/false (defaut true)
%   'liquidStages'        : cellstr ou string array. Par defaut : toutes
%                           les etapes presentes dans les dumps.
%   'verticalCoordinate'  : 'y' (defaut) ou 'depth'
%
% Remarque:
%   Les profils verticaux moyens sont toujours des moyennes selon x.
%   Pour les champs moyens temporels, Rh est reconstruit a partir des champs
%   moyens accumules. Ce n'est pas strictement la moyenne temporelle de Rh,
%   mais c'est tres utile pour analyser la tendance hydrostatique moyenne.

p = inputParser;
p.addParameter('state', 'final', @(s) ischar(s) || isstring(s));
p.addParameter('useType', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('plotFigures', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('plotLiquidClosure', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('plotMeanFields', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('plotPoiseuille', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('liquidStages', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('verticalCoordinate', 'y', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

stateName = lower(string(opt.state));
if startsWith(stateName, "liquid_") || startsWith(stateName, "lc_")
    [S, liquidMeta] = resolve_liquid_stage_struct(runOut, char(stateName));
elseif any(strcmp(stateName, ["mean","meanfields","temporalmean","avg","average"]))
    S = resolve_mean_fields_struct(runOut);
    liquidMeta = struct();
else
    S = resolve_state_struct(runOut, char(stateName));
    liquidMeta = struct();
end

if isfield(S,'isMeanFields') && S.isMeanFields
    G = reconstruct_mean_fields_state(S);
else
    assert(isfield(S,'x') && isfield(S,'v'), ...
        'La structure fournie doit contenir x et v.');
    assert(isfield(S,'params') && ~isempty(S.params), ...
        'Impossible de retrouver params dans runOut/simOut.');
    G = reconstruct_single_state(S, logical(opt.useType), liquidMeta);
end

if opt.plotFigures
    plot_single_state_fields(G, opt.verticalCoordinate);
end

lcData = get_liquid_closure_data(runOut);
if ~isempty(lcData.dumps)
    G.liquidClosure = build_liquid_closure_analysis(runOut, lcData, opt);
    if opt.plotFigures && opt.plotLiquidClosure
        plot_liquid_closure_profiles(G.liquidClosure, opt.verticalCoordinate);
    end
else
    G.liquidClosure = struct();
end

meanData = get_mean_fields_data(runOut);
if ~isempty(meanData)
    G.meanFields = reconstruct_mean_fields_state(meanData);
else
    G.meanFields = struct();
end

poiseuilleData = get_poiseuille_data(runOut);
if ~isempty(poiseuilleData)
    poiseuilleRef = extract_poiseuille_reference(poiseuilleData);
    G.poiseuilleRef = poiseuilleRef;
    if isfield(G,'meanFields') && isstruct(G.meanFields) && ~isempty(fieldnames(G.meanFields))
        G.meanFields.poiseuilleRef = poiseuilleRef;
    end
end

if ~isempty(meanData) && opt.plotFigures && opt.plotMeanFields
    plot_mean_fields_analysis(G.meanFields, opt.verticalCoordinate);
end

if ~isempty(poiseuilleData)
    G.poiseuille = build_poiseuille_analysis(G, poiseuilleData);
    if opt.plotFigures && opt.plotPoiseuille
        plot_poiseuille_analysis(G.poiseuille, opt.verticalCoordinate);
    end
else
    G.poiseuille = struct();
end

end

% =====================================================================

function G = reconstruct_single_state(S, useType, liquidMeta)
assert(isfield(S,'x') && isfield(S,'v') && isfield(S,'params'));

x = double(S.x);
v = double(S.v);
params = S.params;

Nx = getf(params, 'Nx', []);
Ny = getf(params, 'Ny', []);
Lx = getf(params, 'Lx', []);
Ly = getf(params, 'Ly', []);
assert(~isempty(Nx) && ~isempty(Ny) && ~isempty(Lx) && ~isempty(Ly), ...
    'params doit contenir Nx, Ny, Lx et Ly.');

dx = Lx / Nx;
dy = Ly / Ny;
Vc = dx * dy;

[ix, iy, ic] = cell_indices_from_positions(x, params);
Nc = Nx * Ny;

countVec = accumarray(ic, 1, [Nc 1], @sum, 0);
sumVx    = accumarray(ic, v(:,1), [Nc 1], @sum, 0);
sumVy    = accumarray(ic, v(:,2), [Nc 1], @sum, 0);

UxVec = zeros(Nc,1);
UyVec = zeros(Nc,1);
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
Tvec(maskT) = sumRel2(maskT) ./ dof(maskT);  % kBT local, coherent avec liquid_closure

rhoVec = countVec ./ Vc;
PkinVec = nan(Nc,1);
PkinVec(maskT) = rhoVec(maskT) .* Tvec(maskT);

Ncell = reshape(countVec, [Ny, Nx]);
Ux    = reshape(UxVec, [Ny, Nx]);
Uy    = reshape(UyVec, [Ny, Nx]);
T     = reshape(Tvec, [Ny, Nx]);
rho   = reshape(rhoVec, [Ny, Nx]);
Pkin  = reshape(PkinVec, [Ny, Nx]);

rhoTarget = [];
if isfield(liquidMeta,'rhoTarget') && ~isempty(liquidMeta.rhoTarget)
    rhoTarget = liquidMeta.rhoTarget;
else
    [~, rhoTarget] = target_maps_from_params(params);
end

Kvirial = getf(params, 'Kvirial', NaN);
if isempty(rhoTarget)
    Pvir = nan(size(Pkin));
    Ptot = nan(size(Pkin));
else
    Pvir = Kvirial .* (rho - rhoTarget);
    Ptot = Pkin + Pvir;
end

xc = ((1:Nx) - 0.5) * dx;
yc = ((1:Ny) - 0.5) * dy;
hc = Ly - yc;
[Xc, Yc] = meshgrid(xc, yc);

G = struct();
G.params = params;
G.dx = dx; G.dy = dy; G.xc = xc; G.yc = yc; G.hc = hc; G.Xc = Xc; G.Yc = Yc;
G.count = Ncell; G.Ncell = Ncell;
G.rho = rho;
G.Ux = Ux; G.Uy = Uy;
G.T = T;
G.P = Pkin;        % compatibilite avec l'ancien script
G.Pkin = Pkin;
G.Pvir = Pvir;
G.Ptot = Ptot;
G.rhoTarget = rhoTarget;
G.Kvirial = Kvirial;
G.isMeanFields = false;

G.nParticles = size(x,1);
G.meanU = [mean(v(:,1)), mean(v(:,2))];
G.meanSpeed2 = mean(sum(v.^2, 2));
G.meanT_nonempty = mean(T(nz_mask_2d(Ncell)), 'omitnan');
G.meanP_nonempty = mean(Pkin(nz_mask_2d(Ncell)), 'omitnan');

G.profiles = compute_vertical_profiles(G);

if isfield(liquidMeta,'PdriveUsed') && ~isempty(liquidMeta.PdriveUsed)
    G.PdriveUsed = liquidMeta.PdriveUsed;
    G.PkinDriveUsed = liquidMeta.PkinDriveUsed;
    G.PvirDriveUsed = liquidMeta.PvirDriveUsed;
    G.duEOSYUsed = liquidMeta.duEOSYUsed;
    G.profiles.PdriveUsed_y = mean(liquidMeta.PdriveUsed, 2, 'omitnan');
    G.profiles.PkinDriveUsed_y = mean(liquidMeta.PkinDriveUsed, 2, 'omitnan');
    G.profiles.PvirDriveUsed_y = mean(liquidMeta.PvirDriveUsed, 2, 'omitnan');
    G.profiles.duEOSYUsed_y = mean(liquidMeta.duEOSYUsed, 2, 'omitnan');
    G.profiles.dPdriveUseddy = dd1_nonperiodic(G.profiles.PdriveUsed_y, G.dy);
    G.profiles.Rh_drive_target = G.profiles.dPdriveUseddy - G.profiles.rhoTarget .* G.profiles.g;
    G.profiles.Rh_drive_state  = G.profiles.dPdriveUseddy - G.profiles.rho .* G.profiles.g;
end

if useType && isfield(S,'type') && ~isempty(S.type)
    typ = double(S.type(:));
    utypes = unique(typ);
    G.byType = struct([]);
    for k = 1:numel(utypes)
        tk = utypes(k);
        mk = (typ == tk);
        if ~any(mk), continue; end
        G.byType(k) = reconstruct_single_type_state(ic, v, mk, Nx, Ny, dx, dy, xc, yc, hc, params); %#ok<AGROW>
        G.byType(k).type = tk;
    end
end
end

function G = reconstruct_mean_fields_state(S)
MF = S.meanFields;
params = S.params;

G = struct();
G.params = params;
G.dx = MF.dx;
G.dy = MF.dy;
G.xc = MF.xc;
G.yc = MF.yc;
G.hc = params.Ly - MF.yc;
G.Xc = MF.Xc;
G.Yc = MF.Yc;
G.Ux = MF.Ux;
G.Uy = MF.Uy;
G.rho = MF.rho;
G.T = MF.T;
G.P = choose_field(MF, 'P', choose_field(MF, 'Pkin', MF.rho .* MF.T));
G.Pkin = choose_field(MF, 'Pkin', G.P);
G.Pvir = choose_field(MF, 'Pvir', G.P - G.Pkin);
G.Ptot = G.Pkin + G.Pvir;
G.rhoTarget = choose_field(MF, 'rhoTarget', nan(size(G.rho)));
G.Kvirial = getf(params, 'Kvirial', NaN);
G.count = [];
G.Ncell = [];
G.isMeanFields = true;
G.nSamples = choose_field(MF, 'nSamples', NaN);
G.meanFieldStartStep = choose_field(MF, 'startStep', NaN);
G.meanFieldStride = choose_field(MF, 'stride', NaN);
G.source = 'simOut.meanFields';
G.profiles = compute_vertical_profiles(G);
end

function val = choose_field(S, name, defaultVal)
if isfield(S,name) && ~isempty(S.(name))
    val = S.(name);
else
    val = defaultVal;
end
end

function Gk = reconstruct_single_type_state(ic, v, mk, Nx, Ny, dx, dy, xc, yc, hc, params)
Nc = Nx * Ny;
ic_k = ic(mk);
vx_k = v(mk,1);
vy_k = v(mk,2);

count_k = accumarray(ic_k, 1, [Nc 1], @sum, 0);
sumVx_k = accumarray(ic_k, vx_k, [Nc 1], @sum, 0);
sumVy_k = accumarray(ic_k, vy_k, [Nc 1], @sum, 0);

Ux_k = zeros(Nc,1); Uy_k = zeros(Nc,1);
nz_k = count_k > 0;
Ux_k(nz_k) = sumVx_k(nz_k) ./ count_k(nz_k);
Uy_k(nz_k) = sumVy_k(nz_k) ./ count_k(nz_k);

uxPart = Ux_k(ic_k);
uyPart = Uy_k(ic_k);
rel2 = (vx_k - uxPart).^2 + (vy_k - uyPart).^2;
sumRel2_k = accumarray(ic_k, rel2, [Nc 1], @sum, 0);
dof_k = 2 * max(count_k - 1, 0);
T_k = nan(Nc,1);
maskT = dof_k > 0;
T_k(maskT) = sumRel2_k(maskT) ./ dof_k(maskT);

N_k = reshape(count_k, [Ny, Nx]);
Ux2 = reshape(Ux_k, [Ny, Nx]);
Uy2 = reshape(Uy_k, [Ny, Nx]);
Tk2 = reshape(T_k, [Ny, Nx]);
rhok = N_k / (dx*dy);
Pk = rhok .* Tk2;
Gk = struct('params',params,'Ncell',N_k,'rho',rhok,'Ux',Ux2,'Uy',Uy2,'T',Tk2,'P',Pk,'Pkin',Pk,'xc',xc,'yc',yc,'hc',hc);
end

function prof = compute_vertical_profiles(G)
prof = struct();
prof.y = G.yc(:);
prof.h = G.hc(:);
prof.rho = mean(G.rho, 2, 'omitnan');
prof.Ux = mean(G.Ux, 2, 'omitnan');
prof.Uy = mean(G.Uy, 2, 'omitnan');
prof.T = mean(G.T, 2, 'omitnan');
prof.Pkin = mean(G.Pkin, 2, 'omitnan');
prof.Pvir = mean(G.Pvir, 2, 'omitnan');
prof.Ptot = mean(G.Ptot, 2, 'omitnan');
if ~isempty(G.rhoTarget)
    prof.rhoTarget = mean(G.rhoTarget, 2, 'omitnan');
else
    prof.rhoTarget = nan(size(prof.y));
end

g = getf(G.params, 'g', NaN);
prof.g = g;
prof.dPtotdy = dd1_nonperiodic(prof.Ptot, G.dy);
prof.Rh_target = prof.dPtotdy - prof.rhoTarget .* g;
prof.Rh_state  = prof.dPtotdy - prof.rho .* g;

if all(isfinite(prof.Pkin))
    prof.dPkindy = dd1_nonperiodic(prof.Pkin, G.dy);
else
    prof.dPkindy = nan(size(prof.y));
end
if all(isfinite(prof.Pvir))
    prof.dPvirdy = dd1_nonperiodic(prof.Pvir, G.dy);
else
    prof.dPvirdy = nan(size(prof.y));
end
end

function plot_single_state_fields(G, verticalCoordinate)
coord = select_vertical_coordinate(G, verticalCoordinate);
[xlab, yplot] = deal(coord.label, coord.values);

nameSuffix = '';
if isfield(G,'isMeanFields') && G.isMeanFields
    nameSuffix = ' - champs moyens';
end

figure('Name',['reconstruct_grid_fields - champs' nameSuffix],'Color','w');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

nexttile; imagesc(G.xc, G.yc, G.rho); set(gca,'YDir','normal'); axis equal tight; colorbar; title('\rho');
nexttile; imagesc(G.xc, G.yc, G.T);   set(gca,'YDir','normal'); axis equal tight; colorbar; title('kBT local');
nexttile; imagesc(G.xc, G.yc, G.Pkin); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{kin} = \rho kBT');
nexttile; imagesc(G.xc, G.yc, G.Pvir); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{vir}');
nexttile; imagesc(G.xc, G.yc, G.Ux);  set(gca,'YDir','normal'); axis equal tight; colorbar; title('U_x');
nexttile; imagesc(G.xc, G.yc, G.Uy);  set(gca,'YDir','normal'); axis equal tight; colorbar; title('U_y');
nexttile; imagesc(G.xc, G.yc, G.Ptot); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{tot} = P_{kin}+P_{vir}');
if ~isempty(G.rhoTarget)
    nexttile; imagesc(G.xc, G.yc, G.rhoTarget); set(gca,'YDir','normal'); axis equal tight; colorbar; title('\rho_{target}');
else
    nexttile; axis off;
end
nexttile; imagesc(G.xc, G.yc, G.Ptot - mean(G.Ptot(:), 'omitnan')); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{tot} - <P_{tot}>');

figure('Name',['reconstruct_grid_fields - profils verticaux' nameSuffix],'Color','w');
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');

nexttile; plot(yplot, G.profiles.Ux, 'LineWidth', 1.5); grid on; xlabel(xlab); ylabel('<U_x>_x'); title('Profil vertical moyen U_x');
nexttile; plot(yplot, G.profiles.Uy, 'LineWidth', 1.5); grid on; xlabel(xlab); ylabel('<U_y>_x'); title('Profil vertical moyen U_y');
nexttile; hold on; plot(yplot, G.profiles.Pkin, 'LineWidth', 1.5, 'DisplayName','P_{kin}');
plot(yplot, G.profiles.Pvir, '--', 'LineWidth', 1.3, 'DisplayName','P_{vir}');
plot(yplot, G.profiles.Ptot, '-.', 'LineWidth', 1.3, 'DisplayName','P_{tot}');
if isfield(G.profiles,'PdriveUsed_y')
    plot(yplot, G.profiles.PdriveUsed_y, 'k:', 'LineWidth', 1.8, 'DisplayName','P_{drive} utilise');
end
grid on; xlabel(xlab); ylabel('Pression'); legend('Location','best'); title('Profils verticaux de pression');
nexttile; hold on; plot(yplot, G.profiles.rho, 'LineWidth', 1.5, 'DisplayName','\rho');
if any(isfinite(G.profiles.rhoTarget))
    plot(yplot, G.profiles.rhoTarget, '--', 'LineWidth', 1.3, 'DisplayName','\rho_{target}');
end
grid on; xlabel(xlab); ylabel('Densite'); legend('Location','best'); title('Profil vertical moyen de densite');
nexttile; plot(yplot, G.profiles.T, 'LineWidth', 1.5); grid on; xlabel(xlab); ylabel('<kBT>_x'); title('Profil vertical moyen de temperature');
nexttile; hold on; plot(yplot, G.profiles.Rh_target, 'LineWidth', 1.5, 'DisplayName','R_h(P_{tot}, \rho_{target})');
plot(yplot, G.profiles.Rh_state, '--', 'LineWidth', 1.3, 'DisplayName','R_h(P_{tot}, \rho)');
if isfield(G.profiles,'Rh_drive_target')
    plot(yplot, G.profiles.Rh_drive_target, 'k:', 'LineWidth', 1.8, 'DisplayName','R_h(P_{drive}, \rho_{target})');
end
yline(0,'k:');
grid on; xlabel(xlab); ylabel('Residuel'); legend('Location','best'); title('Residuel hydrostatique R_h');
nexttile; hold on; plot(yplot, G.profiles.dPtotdy, 'LineWidth', 1.5, 'DisplayName','dP_{tot}/dy');
plot(yplot, G.profiles.rhoTarget .* G.profiles.g, '--', 'LineWidth', 1.3, 'DisplayName','\rho_{target} g');
plot(yplot, G.profiles.rho .* G.profiles.g, '-.', 'LineWidth', 1.3, 'DisplayName','\rho g');
if isfield(G.profiles,'dPdriveUseddy')
    plot(yplot, G.profiles.dPdriveUseddy, 'k:', 'LineWidth', 1.8, 'DisplayName','dP_{drive}/dy');
end
grid on; xlabel(xlab); ylabel('Gradient / forcage'); legend('Location','best'); title('Test hydrostatique');
nexttile; hold on;
plot(yplot, G.profiles.Ux, 'LineWidth', 1.5, 'DisplayName','<U_x>_x');
if isfield(G,'poiseuilleRef') && isfield(G.poiseuilleRef,'uAna') && ~isempty(G.poiseuilleRef.uAna)
    plot(G.poiseuilleRef.yCenters, G.poiseuilleRef.uAna, 'k--', 'LineWidth', 1.5, 'DisplayName','profil analytique');
end
grid on; xlabel(xlab); ylabel('U_x'); legend('Location','best'); title('Poiseuille : U_x(y)');
end

function lc = build_liquid_closure_analysis(runOut, lcData, opt)
params = lcData.params;
allStageNames = fieldnames(lcData.dumps);
requested = normalize_requested_stages(opt.liquidStages, allStageNames);

lc = struct();
lc.stageOrder = requested(:).';
lc.params = params;
lc.Kvirial = getf(params, 'Kvirial', NaN);
lc.stages = struct();

for k = 1:numel(requested)
    name = requested{k};
    S = lcData.dumps.(name);
    S.params = params;
    liquidMeta = struct();
    if strcmpi(name,'reference')
        liquidMeta = extract_liquid_meta_for_reference(lcData, params);
    else
        [~, liquidMeta.rhoTarget] = target_maps_from_params(params);
    end
    Gs = reconstruct_single_state(S, false, liquidMeta);
    lc.stages.(name) = Gs;
end

if isfield(lc.stages,'reference')
    lc.reference = lc.stages.reference;
end
end

function plot_liquid_closure_profiles(lc, verticalCoordinate)
if ~isfield(lc,'stageOrder') || isempty(lc.stageOrder)
    return;
end

firstStage = lc.stages.(lc.stageOrder{1});
coord = select_vertical_coordinate(firstStage, verticalCoordinate);
yplot = coord.values;
xlab = coord.label;

figure('Name','liquid_closure - profils verticaux','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Uy, 'LineWidth', 1.5, 'DisplayName', stage_label(name));
end
xlabel(xlab); ylabel('<U_y>_x'); title('U_y'); grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Pkin, 'LineWidth', 1.5, 'DisplayName', stage_label(name));
end
xlabel(xlab); ylabel('<P_{kin}>_x'); title('P_{kin}'); grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Pvir, 'LineWidth', 1.5, 'DisplayName', stage_label(name));
end
if isfield(lc.stages,'reference') && isfield(lc.stages.reference.profiles,'PvirDriveUsed_y')
    plot(yplot, lc.stages.reference.profiles.PvirDriveUsed_y, 'k:', 'LineWidth', 1.8, 'DisplayName', 'P_{vir} drive utilise');
end
xlabel(xlab); ylabel('<P_{vir}>_x'); title('P_{vir}'); grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Ptot, 'LineWidth', 1.5, 'DisplayName', stage_label(name));
end
if isfield(lc.stages,'reference') && isfield(lc.stages.reference.profiles,'PdriveUsed_y')
    plot(yplot, lc.stages.reference.profiles.PdriveUsed_y, 'k:', 'LineWidth', 1.8, 'DisplayName', 'P_{drive} utilise');
end
xlabel(xlab); ylabel('<P_{tot}>_x'); title('P_{kin}+P_{vir}'); grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.rho, 'LineWidth', 1.5, 'DisplayName', stage_label(name));
end
if isfield(lc.stages,'reference') && any(isfinite(lc.stages.reference.profiles.rhoTarget))
    plot(yplot, lc.stages.reference.profiles.rhoTarget, 'k--', 'LineWidth', 1.5, 'DisplayName', '\rho_{target}');
end
xlabel(xlab); ylabel('<\rho>_x'); title('Densite'); grid on; legend('Location','best');

nexttile; hold on;
if isfield(lc.stages,'reference') && isfield(lc.stages.reference.profiles,'duEOSYUsed_y')
    plot(yplot, lc.stages.reference.profiles.duEOSYUsed_y, 'LineWidth', 1.8, 'DisplayName', '\Delta u_{y,EOS} utilise');
end
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Uy - lc.stages.(lc.stageOrder{1}).profiles.Uy, '--', 'LineWidth', 1.0, 'DisplayName', ['\Delta U_y ' stage_label(name)]);
end
xlabel(xlab); ylabel('Kick / variation'); title('Kick EOS et variations de U_y'); grid on; legend('Location','best');

figure('Name','liquid_closure - test hydrostatique','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Rh_target, 'LineWidth', 1.5, 'DisplayName', [stage_label(name) ' / \rho_{target}']);
end
if isfield(lc.stages,'reference') && isfield(lc.stages.reference.profiles,'Rh_drive_target')
    plot(yplot, lc.stages.reference.profiles.Rh_drive_target, 'k:', 'LineWidth', 1.8, 'DisplayName', 'R_h(P_{drive},\rho_{target})');
end
yline(0,'k:');
xlabel(xlab); ylabel('R_h'); title('Residuel hydrostatique / \rho_{target}'); grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.Rh_state, 'LineWidth', 1.5, 'DisplayName', [stage_label(name) ' / \rho']);
end
if isfield(lc.stages,'reference') && isfield(lc.stages.reference.profiles,'Rh_drive_state')
    plot(yplot, lc.stages.reference.profiles.Rh_drive_state, 'k:', 'LineWidth', 1.8, 'DisplayName', 'R_h(P_{drive},\rho)');
end
yline(0,'k:');
xlabel(xlab); ylabel('R_h'); title('Residuel hydrostatique / \rho'); grid on; legend('Location','best');

nexttile; hold on;
for k = 1:numel(lc.stageOrder)
    name = lc.stageOrder{k};
    Gs = lc.stages.(name);
    plot(yplot, Gs.profiles.dPtotdy, 'LineWidth', 1.5, 'DisplayName', ['dP/dy ' stage_label(name)]);
end
if isfield(lc.stages,'reference') && isfield(lc.stages.reference.profiles,'dPdriveUseddy')
    plot(yplot, lc.stages.reference.profiles.dPdriveUseddy, 'k:', 'LineWidth', 1.8, 'DisplayName', 'dP_{drive}/dy');
end
if isfield(lc.stages,'reference')
    plot(yplot, lc.stages.reference.profiles.rhoTarget .* lc.stages.reference.profiles.g, 'k--', 'LineWidth', 1.5, 'DisplayName', '\rho_{target} g');
end
xlabel(xlab); ylabel('Gradient / forcage'); title('Test hydrostatique'); grid on; legend('Location','best');
end

function plot_mean_fields_analysis(M, verticalCoordinate)
coord = select_vertical_coordinate(M, verticalCoordinate);
yplot = coord.values;
xlab = coord.label;

figure('Name','meanFields - champs moyens temporels','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(M.xc, M.yc, M.rho); set(gca,'YDir','normal'); axis equal tight; colorbar; title(sprintf('\\rho moyen (n=%d)', M.nSamples));
nexttile; imagesc(M.xc, M.yc, M.T); set(gca,'YDir','normal'); axis equal tight; colorbar; title('kBT moyen');
nexttile; imagesc(M.xc, M.yc, M.Pkin); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{kin} moyen');
nexttile; imagesc(M.xc, M.yc, M.Pvir); set(gca,'YDir','normal'); axis equal tight; colorbar; title('P_{vir} moyen');
nexttile; imagesc(M.xc, M.yc, M.Ux); set(gca,'YDir','normal'); axis equal tight; colorbar; title('U_x moyen');
nexttile; imagesc(M.xc, M.yc, M.Uy); set(gca,'YDir','normal'); axis equal tight; colorbar; title('U_y moyen');

figure('Name','meanFields - profils moyens temporels','Color','w');
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile; plot(yplot, M.profiles.Ux, 'LineWidth', 1.5); grid on; xlabel(xlab); ylabel('<U_x>_x'); title('U_x moyen temporel');
nexttile; plot(yplot, M.profiles.Uy, 'LineWidth', 1.5); grid on; xlabel(xlab); ylabel('<U_y>_x'); title('U_y moyen temporel');
nexttile; hold on; plot(yplot, M.profiles.Pkin, 'LineWidth', 1.5, 'DisplayName','P_{kin}'); plot(yplot, M.profiles.Pvir, '--', 'LineWidth', 1.3, 'DisplayName','P_{vir}'); plot(yplot, M.profiles.Ptot, '-.', 'LineWidth', 1.3, 'DisplayName','P_{tot}'); grid on; xlabel(xlab); ylabel('Pression'); legend('Location','best'); title('Pressions moyennes temporelles');
nexttile; hold on; plot(yplot, M.profiles.rho, 'LineWidth', 1.5, 'DisplayName','\rho'); if any(isfinite(M.profiles.rhoTarget)), plot(yplot, M.profiles.rhoTarget, '--', 'LineWidth', 1.3, 'DisplayName','\rho_{target}'); end; grid on; xlabel(xlab); ylabel('Densite'); legend('Location','best'); title('Densite moyenne temporelle');
nexttile; plot(yplot, M.profiles.T, 'LineWidth', 1.5); grid on; xlabel(xlab); ylabel('<kBT>_x'); title('Temperature moyenne temporelle');
nexttile; hold on; plot(yplot, M.profiles.Rh_target, 'LineWidth', 1.5, 'DisplayName','R_h(P_{tot},\rho_{target})'); plot(yplot, M.profiles.Rh_state, '--', 'LineWidth', 1.3, 'DisplayName','R_h(P_{tot},\rho)'); yline(0,'k:'); grid on; xlabel(xlab); ylabel('R_h'); legend('Location','best'); title('Residuel hydrostatique moyen');
nexttile; hold on; plot(yplot, M.profiles.dPtotdy, 'LineWidth', 1.5, 'DisplayName','dP_{tot}/dy'); plot(yplot, M.profiles.rhoTarget .* M.profiles.g, '--', 'LineWidth', 1.3, 'DisplayName','\rho_{target} g'); plot(yplot, M.profiles.rho .* M.profiles.g, '-.', 'LineWidth', 1.3, 'DisplayName','\rho g'); grid on; xlabel(xlab); ylabel('Gradient / forcage'); legend('Location','best'); title('Test hydrostatique moyen');
nexttile; hold on; plot(yplot, M.profiles.Ux, 'LineWidth', 1.5, 'DisplayName','<U_x>_x'); if isfield(M,'poiseuilleRef') && isfield(M.poiseuilleRef,'uAna') && ~isempty(M.poiseuilleRef.uAna), plot(M.poiseuilleRef.yCenters, M.poiseuilleRef.uAna, 'k--', 'LineWidth', 1.5, 'DisplayName','profil analytique'); end; grid on; xlabel(xlab); ylabel('U_x'); legend('Location','best'); title('Poiseuille : U_x(y)');
end


function R = extract_poiseuille_reference(poiseuilleData)
R = struct();
if isfield(poiseuilleData.poiseuille,'yCenters') && ~isempty(poiseuilleData.poiseuille.yCenters)
    R.yCenters = poiseuilleData.poiseuille.yCenters(:);
else
    R.yCenters = [];
end
if isfield(poiseuilleData.poiseuille,'uAna') && ~isempty(poiseuilleData.poiseuille.uAna)
    R.uAna = poiseuilleData.poiseuille.uAna(:);
else
    R.uAna = [];
end
if isfield(poiseuilleData.poiseuille,'uxProf') && ~isempty(poiseuilleData.poiseuille.uxProf)
    R.uxProf = poiseuilleData.poiseuille.uxProf(:);
else
    R.uxProf = [];
end
end

function P = build_poiseuille_analysis(G, poiseuilleData)
P = struct();
P.params = poiseuilleData.params;
P.sim = poiseuilleData.poiseuille;
P.selectedStateLabel = 'etat selectionne';
P.selected = build_poiseuille_profile_summary(G);
if isfield(P.sim,'uxProf') && ~isempty(P.sim.uxProf)
    P.runtimeProfile = struct();
    P.runtimeProfile.y = P.sim.yCenters(:);
    P.runtimeProfile.Ux = P.sim.uxProf(:);
    if isfield(P.sim,'uAna') && ~isempty(P.sim.uAna)
        P.runtimeProfile.uAna = P.sim.uAna(:);
    else
        P.runtimeProfile.uAna = [];
    end
else
    P.runtimeProfile = struct();
end
if isfield(G,'meanFields') && isstruct(G.meanFields) && ~isempty(fieldnames(G.meanFields))
    P.mean = build_poiseuille_profile_summary(G.meanFields);
else
    P.mean = struct();
end
end

function S = build_poiseuille_profile_summary(G)
S = struct();
S.y = G.profiles.y(:);
S.h = G.profiles.h(:);
S.Ux = G.profiles.Ux(:);
S.Uy = G.profiles.Uy(:);
S.UxCenter = max(S.Ux);
S.UxMean = mean(S.Ux, 'omitnan');
S.qFlowProfile = trapz(S.y, S.Ux);
if numel(S.y) >= 3 && all(isfinite(S.Ux))
    p2 = polyfit(S.y, S.Ux, 2);
    S.poly2 = p2;
    S.UxFit = polyval(p2, S.y);
    ssRes = sum((S.Ux - S.UxFit).^2);
    ssTot = sum((S.Ux - mean(S.Ux)).^2);
    S.R2 = 1 - ssRes/max(ssTot, eps);
else
    S.poly2 = [NaN NaN NaN];
    S.UxFit = nan(size(S.Ux));
    S.R2 = NaN;
end
end

function plot_poiseuille_analysis(P, verticalCoordinate)
Gtmp = struct('yc', P.selected.y(:).', 'hc', P.selected.h(:).');
coord = select_vertical_coordinate(Gtmp, verticalCoordinate);
yplot = coord.values;
xlab = coord.label;

figure('Name','poiseuille - champs de vitesse','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile;
axis off;
text(0.0,0.9,sprintf('nu courbure moyen = %.6g', getf(P.sim,'nuCurve_mean',NaN)),'Units','normalized');
text(0.0,0.75,sprintf('nu paroi moyen = %.6g', getf(P.sim,'nuTau_mean',NaN)),'Units','normalized');
text(0.0,0.60,sprintf('Q moyen = %.6g', getf(P.sim,'qFlow_mean',NaN)),'Units','normalized');
text(0.0,0.45,sprintf('R2 moyen = %.6g', getf(P.sim,'R2_mean',NaN)),'Units','normalized');
title('Resume Poiseuille');
nexttile;
if isfield(P,'mean') && ~isempty(P.mean)
    plot(P.mean.Ux, yplot, 'LineWidth', 1.8); grid on;
    xlabel('U_x'); ylabel(xlab); title('Profil moyen temporel U_x(y)');
else
    plot(P.selected.Ux, yplot, 'LineWidth', 1.8); grid on;
    xlabel('U_x'); ylabel(xlab); title('Profil selectionne U_x(y)');
end
set(gca,'YDir','normal');
nexttile; hold on;
plot(P.selected.Ux, yplot, 'LineWidth', 1.5, 'DisplayName','etat selectionne');
if isfield(P.selected,'UxFit') && ~isempty(P.selected.UxFit)
    plot(P.selected.UxFit, yplot, '--', 'LineWidth', 1.3, 'DisplayName','fit quadratique');
end
if isfield(P,'mean') && ~isempty(P.mean)
    plot(P.mean.Ux, yplot, 'LineWidth', 1.5, 'DisplayName','moyenne temporelle');
    if isfield(P.mean,'UxFit') && ~isempty(P.mean.UxFit)
        plot(P.mean.UxFit, yplot, '--', 'LineWidth', 1.3, 'DisplayName','fit quad moyen');
    end
end
if isfield(P,'runtimeProfile') && isfield(P.runtimeProfile,'Ux') && ~isempty(P.runtimeProfile.Ux)
    plot(P.runtimeProfile.Ux, P.runtimeProfile.y, ':', 'LineWidth', 1.8, 'DisplayName','uxProf simOut');
    if isfield(P.runtimeProfile,'uAna') && ~isempty(P.runtimeProfile.uAna)
        plot(P.runtimeProfile.uAna, P.runtimeProfile.y, 'k-.', 'LineWidth', 1.5, 'DisplayName','profil analytique');
    end
end
grid on; xlabel('U_x'); ylabel(xlab); title('Comparaison des profils U_x(y)'); legend('Location','best');
set(gca,'YDir','normal');
end

function S = get_poiseuille_data(runOut)
S = struct();
P = [];
params = [];
if isfield(runOut,'simOut') && isfield(runOut.simOut,'poiseuille') && ~isempty(runOut.simOut.poiseuille)
    P = runOut.simOut.poiseuille;
    if isfield(runOut,'simOut') && isfield(runOut.simOut,'params')
        params = runOut.simOut.params;
    end
elseif isfield(runOut,'poiseuille') && ~isempty(runOut.poiseuille)
    P = runOut.poiseuille;
    if isfield(runOut,'params')
        params = runOut.params;
    end
end
if isempty(P)
    return;
end
S.poiseuille = P;
S.params = params;
end

function [S, liquidMeta] = resolve_liquid_stage_struct(runOut, stateName)
lcData = get_liquid_closure_data(runOut);
assert(~isempty(lcData.dumps), 'Aucun dump liquid_closure trouve dans runOut/simOut.');
name = erase(string(stateName), "lc_");
name = erase(name, "liquid_");
name = char(name);
name = normalize_stage_name(name);
assert(isfield(lcData.dumps, name), 'Etape liquid_closure absente: %s', name);
S = lcData.dumps.(name);
S.params = lcData.params;
if strcmpi(name,'reference')
    liquidMeta = extract_liquid_meta_for_reference(lcData, lcData.params);
else
    liquidMeta = struct();
    [~, liquidMeta.rhoTarget] = target_maps_from_params(lcData.params);
end
end

function S = resolve_mean_fields_struct(runOut)
meanData = get_mean_fields_data(runOut);
assert(~isempty(meanData), 'Aucun champ moyen temporel trouve dans runOut/simOut.meanFields.');
S = meanData;
end

function lcData = get_liquid_closure_data(runOut)
lcData = struct('dumps', struct(), 'closure', struct(), 'params', []);
if isfield(runOut,'simOut') && isfield(runOut.simOut,'params')
    lcData.params = runOut.simOut.params;
elseif isfield(runOut,'params')
    lcData.params = runOut.params;
end

if isfield(runOut,'simOut') && isfield(runOut.simOut,'lastLiquidClosureDump')
    lcData.dumps = runOut.simOut.lastLiquidClosureDump;
elseif isfield(runOut,'lastLiquidClosureDump')
    lcData.dumps = runOut.lastLiquidClosureDump;
elseif isfield(runOut,'simOut') && isfield(runOut.simOut,'lastLiquidClosure') && isfield(runOut.simOut.lastLiquidClosure,'dumps')
    lcData.dumps = runOut.simOut.lastLiquidClosure.dumps;
end

if isfield(runOut,'simOut') && isfield(runOut.simOut,'lastLiquidClosure')
    lcData.closure = runOut.simOut.lastLiquidClosure;
elseif isfield(runOut,'lastLiquidClosure')
    lcData.closure = runOut.lastLiquidClosure;
end
end

function S = get_mean_fields_data(runOut)
S = struct();
MF = [];
params = [];
if isfield(runOut,'simOut') && isfield(runOut.simOut,'meanFields') && ~isempty(runOut.simOut.meanFields)
    MF = runOut.simOut.meanFields;
    if isfield(runOut,'simOut') && isfield(runOut.simOut,'params')
        params = runOut.simOut.params;
    end
elseif isfield(runOut,'meanFields') && ~isempty(runOut.meanFields)
    MF = runOut.meanFields;
    if isfield(runOut,'params')
        params = runOut.params;
    end
end
if isempty(MF)
    return;
end
if isempty(params) && isfield(MF,'params')
    params = MF.params;
end
S.meanFields = MF;
S.params = params;
S.isMeanFields = true;
end

function liquidMeta = extract_liquid_meta_for_reference(lcData, params)
liquidMeta = struct();
[~, liquidMeta.rhoTarget] = target_maps_from_params(params);
if isfield(lcData,'closure') && ~isempty(lcData.closure)
    C = lcData.closure;
    if isfield(C,'rhoTarget') && ~isempty(C.rhoTarget), liquidMeta.rhoTarget = C.rhoTarget; end
    if isfield(C,'Pdrive') && ~isempty(C.Pdrive), liquidMeta.PdriveUsed = C.Pdrive; end
    if isfield(C,'PkinDrive') && ~isempty(C.PkinDrive), liquidMeta.PkinDriveUsed = C.PkinDrive; end
    if isfield(C,'PvirDrive') && ~isempty(C.PvirDrive), liquidMeta.PvirDriveUsed = C.PvirDrive; end
    if isfield(C,'duEOSY') && ~isempty(C.duEOSY), liquidMeta.duEOSYUsed = C.duEOSY; end
end
end

function requested = normalize_requested_stages(requestedIn, available)
if isempty(requestedIn)
    requested = available(:).';
    requested = reorder_stage_names(requested);
    return;
end
if ischar(requestedIn) || isstring(requestedIn)
    requestedIn = cellstr(requestedIn);
end
requested = cellfun(@(s) normalize_stage_name(char(string(s))), requestedIn, 'UniformOutput', false);
requested = requested(ismember(requested, available));
requested = reorder_stage_names(requested);
end

function names = reorder_stage_names(names)
preferred = {'reference','redistributed','repaired','postKick','preReorientation','postReorientation'};
out = {};
for k = 1:numel(preferred)
    if any(strcmp(names, preferred{k})), out{end+1} = preferred{k}; end %#ok<AGROW>
end
for k = 1:numel(names)
    if ~any(strcmp(out, names{k})), out{end+1} = names{k}; end %#ok<AGROW>
end
names = out;
end

function name = normalize_stage_name(name)
s = lower(strrep(strrep(name,'_',''),'-',''));
switch s
    case {'reference','ref'}
        name = 'reference';
    case {'redistributed','redistribution','redist'}
        name = 'redistributed';
    case {'repaired','repair','repare'}
        name = 'repaired';
    case {'postkick','eos','posteos'}
        name = 'postKick';
    case {'prereorientation','prereorient','prereorientation'}
        name = 'preReorientation';
    case {'postreorientation','postreorient'}
        name = 'postReorientation';
    otherwise
        name = name;
end
end

function label = stage_label(name)
switch name
    case 'reference', label = 'reference';
    case 'redistributed', label = 'redistributed';
    case 'repaired', label = 'repaired';
    case 'postKick', label = 'postKick';
    case 'preReorientation', label = 'preReorientation';
    case 'postReorientation', label = 'postReorientation';
    otherwise, label = name;
end
end

function coord = select_vertical_coordinate(G, verticalCoordinate)
mode = lower(string(verticalCoordinate));
switch mode
    case "depth"
        coord.values = G.hc(:);
        coord.label = 'h = Ly - y';
    otherwise
        coord.values = G.yc(:);
        coord.label = 'y';
end
end

function [targetOcc, rhoTarget] = target_maps_from_params(params)
Nx = params.Nx;
Ny = params.Ny;
aX = params.Lx / params.Nx;
aY = params.Ly / params.Ny;
Vc = aX * aY;

gammaEff = getf(params, 'gammaCurrent', getf(params,'gamma', []));
activeFrac = active_cell_fraction_map_light(params, [Ny, Nx]);
targetOcc = gammaEff .* activeFrac;
rhoTarget = targetOcc ./ Vc;
end

function fracMap = active_cell_fraction_map_light(params, sz)
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
    margin = getf(params, 'pistonActiveMargin', 0.0);
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

function [ix, iy, ic] = cell_indices_from_positions(x, params)
Nx = params.Nx; Ny = params.Ny;
dx = params.Lx / params.Nx;
dy = params.Ly / params.Ny;
ix = floor(x(:,1) / dx) + 1;
iy = floor(x(:,2) / dy) + 1;
ix = max(1, min(Nx, ix));
iy = max(1, min(Ny, iy));
ic = sub2ind([Ny, Nx], iy, ix);
end

function S = resolve_state_struct(runOut, stateName)
stateName = lower(string(stateName));

if isfield(runOut,'x') && isfield(runOut,'v') && isfield(runOut,'params')
    S = runOut;
    return;
end

assert(isstruct(runOut), 'L''entree doit etre une structure.');
assert(isfield(runOut,'simOut') || isfield(runOut,'finalState') || isfield(runOut,'initialState'), ...
    'Structure non reconnue : ni simOut, ni finalState, ni initialState.');

switch stateName
    case "simout"
        assert(isfield(runOut,'simOut') && ~isempty(runOut.simOut), 'runOut.simOut absent.');
        S = runOut.simOut;
    case "final"
        if isfield(runOut,'finalState') && ~isempty(runOut.finalState)
            S = runOut.finalState;
            if ~isfield(S,'params') || isempty(S.params)
                if isfield(runOut,'params') && ~isempty(runOut.params)
                    S.params = runOut.params;
                elseif isfield(runOut,'simOut') && isfield(runOut.simOut,'params')
                    S.params = runOut.simOut.params;
                end
            end
        elseif isfield(runOut,'simOut') && ~isempty(runOut.simOut)
            S = runOut.simOut;
        else
            error('Etat final introuvable.');
        end
    case "initial"
        assert(isfield(runOut,'initialState') && ~isempty(runOut.initialState), 'runOut.initialState absent.');
        S = runOut.initialState;
        if ~isfield(S,'params') || isempty(S.params)
            if isfield(runOut,'params') && ~isempty(runOut.params)
                S.params = runOut.params;
            elseif isfield(runOut,'simOut') && isfield(runOut.simOut,'params')
                S.params = runOut.simOut.params;
            end
        end
    otherwise
        error('Option state inconnue: %s', char(stateName));
end

if ~isfield(S,'params') || isempty(S.params)
    if isfield(runOut,'params') && ~isempty(runOut.params)
        S.params = runOut.params;
    elseif isfield(runOut,'simOut') && isfield(runOut.simOut,'params')
        S.params = runOut.simOut.params;
    else
        error('Impossible de retrouver params.');
    end
end
end

function v = getf(s, name, defaultVal)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultVal;
end
end

function m = nz_mask_2d(A)
m = A > 0;
end

function df = dd1_nonperiodic(f, dy)
% Derivee 1D non periodique :
% - difference centree a l'interieur
% - difference decentree d'ordre 1 aux bords
f = f(:);
n = numel(f);
df = nan(size(f));

if n == 0
    return;
elseif n == 1
    df(:) = 0;
    return;
elseif n == 2
    val = (f(2) - f(1)) / dy;
    df(1) = val;
    df(2) = val;
    return;
end

df(1) = (f(2) - f(1)) / dy;
df(end) = (f(end) - f(end-1)) / dy;
df(2:end-1) = (f(3:end) - f(1:end-2)) / (2*dy);
end
