function projection_step_visualize_frame(state, sampleDiag, params, it)
%PROJECTION_STEP_VISUALIZE_FRAME Live visualization for the step-channel case.
%
%   projection_step_visualize_frame(state, sampleDiag, params, it)
%
% Four panels are shown:
%   1) particles colored by speed;
%   2) relative cell population (N/gamma - 1);
%   3) velocity magnitude with quiver overlay;
%   4) vorticity omega.
%
% The rectangular step is drawn on every panel.  This function is intended
% for runtime diagnostics; it is lightweight enough for visualEvery of order
% 25--100 on the 48x24 reference grid.

if nargin < 4
    it = sampleDiag.it;
end

figId = get_param(params, 'visualFigureId', 101);
figName = get_param(params, 'visualFigureName', 'Step-channel live visualization');
fig = figure(figId);
clf(fig);
set(fig, 'Name', figName, 'Color', 'w');

G = sampleDiag.G;
step = sampleDiag.step;
solidMask = false(size(G.N));
if isfield(step, 'solidMask') && ~isempty(step.solidMask)
    solidMask = logical(step.solidMask);
elseif isfield(params, 'solidMask') && ~isempty(params.solidMask)
    solidMask = logical(params.solidMask);
end

Nx = size(G.N, 1);
Ny = size(G.N, 2);
xc = ((1:Nx) - 0.5) * G.dx;
yc = ((1:Ny) - 0.5) * G.dy;
[Xc, Yc] = ndgrid(xc, yc);

N = double(G.N);
relN = N ./ max(params.gamma, eps) - 1.0;
relN(solidMask) = NaN;
Ux = double(G.Ux);
Uy = double(G.Uy);
Ux(solidMask) = NaN;
Uy(solidMask) = NaN;
speed = sqrt(Ux.^2 + Uy.^2);
omega = double(step.omega);
omega(solidMask) = NaN;

methodLabel = char(string(get_param(params, 'visualTitleSuffix', '')));
if ~isempty(methodLabel)
    methodLabel = [' - ' methodLabel];
end
mainTitle = sprintf('step %d, t=%.4g%s', it, sampleDiag.t, methodLabel);

% --- Panel 1: particles -------------------------------------------------
subplot(2,2,1);
plot_particles(state, params);
title(['particles ' mainTitle], 'Interpreter', 'none');

% --- Panel 2: density ---------------------------------------------------
subplot(2,2,2);
imagesc(xc, yc, relN.');
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
climSym = get_param(params, 'visualDensityCLim', 0.5);
if isfinite(climSym) && climSym > 0
    caxis([-climSym climSym]);
end
hold on;
draw_step_patch(params);
title('density N/\gamma - 1');
xlabel('x'); ylabel('y');

% --- Panel 3: velocity --------------------------------------------------
subplot(2,2,3);
imagesc(xc, yc, speed.');
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
velCLim = get_param(params, 'visualSpeedCLim', NaN);
if isfinite(velCLim) && velCLim > 0
    caxis([0 velCLim]);
end
hold on;
strideX = max(1, round(get_param(params, 'visualQuiverStrideX', 4)));
strideY = max(1, round(get_param(params, 'visualQuiverStrideY', 3)));
ix = 1:strideX:Nx;
iy = 1:strideY:Ny;
quiverScale = get_param(params, 'visualQuiverScale', 1.5);
quiver(Xc(ix,iy), Yc(ix,iy), Ux(ix,iy), Uy(ix,iy), quiverScale, 'k');
draw_step_patch(params);
title('|u| + quiver');
xlabel('x'); ylabel('y');

% --- Panel 4: vorticity -------------------------------------------------
subplot(2,2,4);
imagesc(xc, yc, omega.');
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
omegaCLim = get_param(params, 'visualOmegaCLim', NaN);
if ~(isfinite(omegaCLim) && omegaCLim > 0)
    omegaCLim = robust_symmetric_clim(omega, get_param(params, 'visualOmegaClipSigma', 3.0));
end
if isfinite(omegaCLim) && omegaCLim > 0
    caxis([-omegaCLim omegaCLim]);
end
hold on;
draw_step_patch(params);
title('\omega = \partial_x u_y - \partial_y u_x');
xlabel('x'); ylabel('y');

sgtitle(sprintf('%s | meanUx=%.4g, enst=%.4g, recircL=%.4g', ...
    mainTitle, sampleDiag.meanUxFluid, step.meanEnstrophy, step.recirculationLength), ...
    'Interpreter', 'none');

drawnow limitrate;
pauseTime = get_param(params, 'visualPause', 0.0);
if isfinite(pauseTime) && pauseTime > 0
    pause(pauseTime);
end

if logical(get_param(params, 'visualSaveFrames', false))
    save_frame(fig, params, it);
end
end

function plot_particles(state, params)
Np = size(state.x, 1);
maxParticles = get_param(params, 'visualMaxParticles', 3000);
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
hold on;
draw_step_patch(params);
xlabel('x'); ylabel('y');
end

function draw_step_patch(params)
x0 = get_param(params, 'stepX0', 0.15 * params.Lx);
x1 = get_param(params, 'stepX1', 0.35 * params.Lx);
h = get_param(params, 'stepHeight', 0.25 * params.Ly);
patch([x0 x1 x1 x0], [0 0 h h], [0.1 0.1 0.1], ...
    'FaceAlpha', 0.35, 'EdgeColor', 'k', 'LineWidth', 1.0);
end

function c = robust_symmetric_clim(A, clipSigma)
vals = A(isfinite(A));
if isempty(vals)
    c = 1;
    return;
end
rmsVal = sqrt(mean(vals.^2));
stdVal = std(vals);
c = clipSigma * max(stdVal, rmsVal / max(clipSigma, eps));
if ~isfinite(c) || c <= 0
    c = max(abs(vals));
end
if ~isfinite(c) || c <= 0
    c = 1;
end
end

function save_frame(fig, params, it)
outDir = get_param(params, 'visualFrameDir', 'step_visual_frames');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
filePrefix = char(string(get_param(params, 'visualFramePrefix', 'step_frame')));
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
