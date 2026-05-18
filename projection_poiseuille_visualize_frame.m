function keepRunning = projection_poiseuille_visualize_frame(state, G, sampleDiag, outPartial, params, it)
%PROJECTION_POISEUILLE_VISUALIZE_FRAME Live monitoring dashboard for Poiseuille.
%
%   keepRunning = projection_poiseuille_visualize_frame(state, G, sampleDiag, outPartial, params, it)
%
% Displays a 2x3 dashboard:
%   1) particles colored by speed;
%   2) relative density N/gamma - 1;
%   3) Ux field + quiver of (Ux,Uy);
%   4) instantaneous profile, running mean profile, quadratic fit;
%   5) kBT, low-k density and rms(div u) histories;
%   6) center-wall velocity and nu_eff history.
%
% A "Stop run" button is created at the top-right of the figure. Clicking it
% requests a graceful stop; the function then returns keepRunning=false.

if nargin < 6
    it = sampleDiag.step;
end
keepRunning = true;

figId = get_param(params, 'visualFigureId', 510);
figName = get_param(params, 'visualFigureName', 'Poiseuille live monitoring');
fig = figure(figId);
% clf does not reliably clear figure appdata. Reset the stop flag at the
% first visualization of a new run so a previous stopped run cannot abort
% immediately at t=0. Do not reset it at later frames, otherwise clicking
% Stop between refreshes would be ignored.
if it == 0 || ~isappdata(fig, 'stopRequested')
    setappdata(fig, 'stopRequested', false);
end
clf(fig);
set(fig, 'Name', figName, 'Color', 'w');
create_stop_button(fig);

T = outPartial.diagTable;
fitTable = table();
if isfield(outPartial, 'fitTable')
    fitTable = outPartial.fitTable;
end
visc = [];
if isfield(outPartial, 'viscosity')
    visc = outPartial.viscosity;
end
if isempty(visc)
    visc = default_visc(outPartial.yCenters);
end

Nx = size(G.N, 1);
Ny = size(G.N, 2);
xc = ((1:Nx) - 0.5) * G.dx;
yc = ((1:Ny) - 0.5) * G.dy;
[Xc, Yc] = ndgrid(xc, yc);

N = double(G.N);
relN = N ./ max(params.gamma, eps) - 1.0;
Ux = double(G.Ux);
Uy = double(G.Uy);
UxYInst = mean(Ux, 1, 'omitnan').';

methodLabel = char(string(get_param(params, 'visualTitleSuffix', '')));
if isempty(methodLabel)
    methodLabel = upper(char(string(params.method)));
else
    methodLabel = sprintf('%s - %s', upper(char(string(params.method))), methodLabel);
end

statusStr = local_status_string(T, fitTable, params);
[currentNuRaw, currentNuScaled, currentNuDisplay, currentNuRefNy] = current_nu_values(fitTable, visc, params);
currentR2 = last_finite_from_table(fitTable, 'R2');
if ~isfinite(currentR2) && isfield(visc,'R2'), currentR2 = visc.R2; end
acc = poiseuille_acceleration_diagnostics(outPartial, get_param(params, 'visualAccelerationWindowTime', 5.0));
UmaxNow = max(abs(Ux(:)), [], 'omitnan');
if isfinite(currentNuScaled) && isfinite(currentNuRaw) && abs(currentNuScaled-currentNuRaw) > 10*eps(max(1,abs(currentNuRaw)))
    nuText = sprintf('nuRaw=%.4g | nu@Ny%d=%.4g', currentNuRaw, currentNuRefNy, currentNuScaled);
else
    nuText = sprintf('nu=%.4g', currentNuDisplay);
end
sg = sprintf('Poiseuille %s | step %d, t=%.4g | %s | Umax=%.4g | accelRatio=%.3g | kBT=%.4g | lowK=%.3e | %s | R2=%.3f', ...
    methodLabel, it, sampleDiag.t, statusStr, UmaxNow, acc.ratioAccel, sampleDiag.kBTCell, sampleDiag.lowKDensityEnergy, nuText, currentR2);

% --- Panel 1: particles -------------------------------------------------
subplot(2,3,1);
plot_particles(state, params);
title('particles (speed color)');

% --- Panel 2: density ---------------------------------------------------
subplot(2,3,2);
imagesc(xc, yc, relN.');
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
climSym = get_param(params, 'visualDensityCLim', 0.5);
if isfinite(climSym) && climSym > 0
    caxis([-climSym climSym]);
end
title('density N/\gamma - 1');
xlabel('x'); ylabel('y');

% --- Panel 3: velocity field -------------------------------------------
subplot(2,3,3);
imagesc(xc, yc, Ux.');
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
uxCLim = get_param(params, 'visualUxCLim', NaN);
if isfinite(uxCLim) && uxCLim > 0
    caxis([-uxCLim uxCLim]);
end
hold on;
strideX = max(1, round(get_param(params, 'visualQuiverStrideX', 4)));
strideY = max(1, round(get_param(params, 'visualQuiverStrideY', 3)));
ix = 1:strideX:Nx;
iy = 1:strideY:Ny;
quiverScale = get_param(params, 'visualQuiverScale', 1.0);
quiver(Xc(ix,iy), Yc(ix,iy), Ux(ix,iy), Uy(ix,iy), quiverScale, 'k');
title(sprintf('U_x field + quiver, max|U_y|=%.3g', max(abs(Uy(:)), [], 'omitnan')));
xlabel('x'); ylabel('y');

% --- Panel 4: profile and fit ------------------------------------------
subplot(2,3,4);
y = outPartial.yCenters(:);
plot(UxYInst, y, '-', 'LineWidth', 1.2, 'DisplayName', 'instant');
hold on;
if isfield(visc,'UxMean') && ~isempty(visc.UxMean)
    plot(visc.UxMean, y, '-', 'LineWidth', 1.8, 'DisplayName', 'mean fit-window');
end
if isfield(visc,'UxFit') && ~isempty(visc.UxFit)
    plot(visc.UxFit, y, '--', 'LineWidth', 1.8, 'DisplayName', 'quadratic fit');
end
grid on;
ylabel('y'); xlabel('U_x(y)');
title(sprintf('profile | center-wall=%.4g', sampleDiag.centerMinusWall));
legend('Location','best');

% --- Panel 5: thermal / low-k / div histories --------------------------
subplot(2,3,5);
tvec = T.t;
yyaxis left;
plot(tvec, T.kBTCell, '-', 'LineWidth', 1.5, 'DisplayName', 'kBT cell');
hold on;
yline(params.kBT, '--', 'DisplayName', 'kBT target');
ylabel('kBT cell');
grid on;
yyaxis right;
set(gca, 'YScale', 'log');
floorVal = get_param(params, 'visualLowKFloor', 1e-12);
plot(tvec, max(T.lowKDensityEnergy, floorVal), '-', 'LineWidth', 1.5, 'DisplayName', 'low-k dens');
hold on;
plot(tvec, max(T.rmsDivU, floorVal), '--', 'LineWidth', 1.2, 'DisplayName', 'rms divU');
ylabel('low-k / divU (log)');
xlabel('t');
curOutBand = last_finite(T.popOutBand);
curPopStd = last_finite(T.popStd);
title(sprintf('thermal & compressive metrics | outBand=%.3g, popStd=%.3g', curOutBand, curPopStd));
legend('Location','best');

% --- Panel 6: center-wall and viscosity histories ----------------------
subplot(2,3,6);
yyaxis left;
plot(tvec, T.centerMinusWall, '-', 'LineWidth', 1.5, 'DisplayName', 'center-wall');
ylabel('center-wall U_x');
hold on;
grid on;
yyaxis right;
if ~isempty(fitTable)
    fitT = fitTable.t;
    if ismember('nuEffScaledNyRef', fitTable.Properties.VariableNames) && logical(get_param(params, 'visualUseScaledNu', true))
        plot(fitT, fitTable.nuEffScaledNyRef, '-', 'LineWidth', 1.6, 'DisplayName', sprintf('nu@Ny%d', currentNuRefNy));
        hold on;
        if ismember('nuEffRaw', fitTable.Properties.VariableNames)
            plot(fitT, fitTable.nuEffRaw, '--', 'LineWidth', 1.1, 'DisplayName', 'nu raw');
        else
            plot(fitT, fitTable.nuEff, '--', 'LineWidth', 1.1, 'DisplayName', 'nu raw');
        end
    else
        plot(fitT, fitTable.nuEff, '-', 'LineWidth', 1.6, 'DisplayName', 'nu raw');
        hold on;
    end
    if isfinite(get_param(params, 'visualReferenceNu', NaN))
        yline(get_param(params, 'visualReferenceNu', NaN), '--', 'DisplayName', 'nu ref');
    end
end
ylabel('nu fit');
xlabel('t');
curSNR = last_finite_from_table(fitTable, 'signalToNoise');
if isfield(visc,'signalToNoise') && ~isfinite(curSNR), curSNR = visc.signalToNoise; end
title(sprintf('profile evolution | R2=%.3f, SNR=%.3g, accel=%.3g', currentR2, curSNR, acc.ratioAccel));
legend('Location','best');

sgtitle(sg, 'Interpreter', 'none');

drawnow limitrate;
pauseTime = get_param(params, 'visualPause', 0.0);
if isfinite(pauseTime) && pauseTime > 0
    pause(pauseTime);
end

if logical(get_param(params, 'visualSaveFrames', false))
    save_frame(fig, params, it);
end

if isvalid(fig)
    if isappdata(fig, 'stopRequested') && logical(getappdata(fig, 'stopRequested'))
        keepRunning = false;
    end
else
    keepRunning = false;
end
end

function create_stop_button(fig)
if ~isvalid(fig)
    return;
end
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop run', ...
    'Units', 'normalized', 'Position', [0.90 0.955 0.085 0.035], ...
    'FontWeight', 'bold', 'ForegroundColor', [0.8 0 0], ...
    'Callback', @(src,evt) setappdata(ancestor(src,'figure'), 'stopRequested', true));
end

function plot_particles(state, params)
Np = size(state.x, 1);
maxParticles = get_param(params, 'visualMaxParticles', 6000);
if isempty(maxParticles) || ~isfinite(maxParticles) || maxParticles <= 0 || maxParticles >= Np
    idx = 1:Np;
else
    idx = unique(round(linspace(1, Np, maxParticles)));
end
speedP = sqrt(sum(state.v(idx,:).^2, 2));
scatter(state.x(idx,1), state.x(idx,2), get_param(params, 'visualParticleSize', 5), speedP, 'filled');
axis equal tight;
xlim([0 params.Lx]);
ylim([0 params.Ly]);
colorbar;
xlabel('x'); ylabel('y');
end

function visc = default_visc(y)
visc = struct('tMin',NaN,'tMax',NaN,'nProfiles',NaN, ...
    'nuEff',NaN,'nuEffRaw',NaN,'nuEffScaledNyRef',NaN,'nuEffDisplay',NaN, ...
    'nuEffCellY',NaN,'nuScaleFactorNyRef',NaN,'nuReferenceNy',NaN,'dy',NaN, ...
    'R2',NaN,'signalToNoise',NaN,'centerMinusWall',NaN,'UxMean',nan(size(y)), ...
    'UxFit',nan(size(y)),'y',y,'fitMask',true(size(y)));
end

function status = local_status_string(T, fitTable, params)
status = 'evolving';
if isempty(T)
    return;
end
kBT = last_finite(T.kBTCell);
lowK = last_finite(T.lowKDensityEnergy);
rmsDiv = last_finite(T.rmsDivU);
outBand = last_finite(T.popOutBand);
if kBT > get_param(params, 'visualDangerKBT', 2*params.kBT) || ...
        lowK > get_param(params, 'visualDangerLowK', 1e-2) || ...
        rmsDiv > get_param(params, 'visualDangerRmsDivU', 1.0) || ...
        outBand > get_param(params, 'visualDangerOutBand', 0.8)
    status = 'DANGER';
    return;
end
m = max(3, round(get_param(params, 'visualStatusTrendWindow', 5)));
if height(T) >= m
    idx = max(1, height(T)-m+1):height(T);
    dt = max(T.t(idx(end)) - T.t(idx(1)), eps);
    slopeCW = (T.centerMinusWall(idx(end)) - T.centerMinusWall(idx(1))) / dt;
    lowK0 = max(T.lowKDensityEnergy(idx(1)), get_param(params,'visualLowKFloor',1e-12));
    lowK1 = max(T.lowKDensityEnergy(idx(end)), get_param(params,'visualLowKFloor',1e-12));
    growthLowK = lowK1 / lowK0;
    slopeThresh = get_param(params, 'visualSteadySlopeThreshold', 5e-4);
    growthThresh = get_param(params, 'visualSteadyLowKGrowthThreshold', 2.0);
    curR2 = last_finite_from_table(fitTable, 'R2');
    if isfinite(curR2) && curR2 > 0.8 && abs(slopeCW) < slopeThresh && growthLowK < growthThresh
        status = 'quasi-steady';
    elseif abs(slopeCW) < 5*slopeThresh
        status = 'slow-evolving';
    else
        status = 'evolving';
    end
end
end

function val = last_finite(x)
val = NaN;
idx = find(isfinite(x), 1, 'last');
if ~isempty(idx)
    val = x(idx);
end
end


function [nuRaw, nuScaled, nuDisplay, refNy] = current_nu_values(fitTable, visc, params)
nuRaw = last_finite_from_table(fitTable, 'nuEffRaw');
if ~isfinite(nuRaw)
    nuRaw = last_finite_from_table(fitTable, 'nuEff');
end
nuScaled = last_finite_from_table(fitTable, 'nuEffScaledNyRef');
refNy = last_finite_from_table(fitTable, 'nuReferenceNy');
if ~isfinite(nuRaw) && isfield(visc, 'nuEffRaw')
    nuRaw = visc.nuEffRaw;
end
if ~isfinite(nuRaw) && isfield(visc, 'nuEff')
    nuRaw = visc.nuEff;
end
if ~isfinite(nuScaled) && isfield(visc, 'nuEffScaledNyRef')
    nuScaled = visc.nuEffScaledNyRef;
end
if ~isfinite(refNy) && isfield(visc, 'nuReferenceNy')
    refNy = visc.nuReferenceNy;
end
if ~isfinite(refNy)
    refNy = get_param(params, 'poiseuilleNuReferenceNy', NaN);
end
if ~isfinite(refNy)
    refNy = NaN;
end
if logical(get_param(params, 'visualUseScaledNu', true)) && isfinite(nuScaled)
    nuDisplay = nuScaled;
else
    nuDisplay = nuRaw;
end
if ~isfinite(nuScaled)
    nuScaled = nuRaw;
end
end

function val = last_finite_from_table(T, varname)
val = NaN;
if isempty(T) || ~istable(T) || ~ismember(varname, T.Properties.VariableNames)
    return;
end
val = last_finite(T.(varname));
end

function save_frame(fig, params, it)
outDir = get_param(params, 'visualFrameDir', 'poiseuille_live_frames');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
filePrefix = char(string(get_param(params, 'visualFramePrefix', 'poiseuille')));
fileName = fullfile(outDir, sprintf('%s_%06d.png', filePrefix, it));
try
    exportgraphics(fig, fileName, 'Resolution', get_param(params, 'visualFrameResolution', 120));
catch
    saveas(fig, fileName);
end
end

function value = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
