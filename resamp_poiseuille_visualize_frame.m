function resamp_poiseuille_visualize_frame(state, params, prof, md, varargin)
%RESAMP_POISEUILLE_VISUALIZE_FRAME Live visualisation for weighted Poiseuille.
%
% Two figures are produced: a spatial diagnostic figure and a profile/fit
% figure.  This function does not alter the state and is safe to call at low
% cadence during long runs.

opts = parse_options(varargin{:});
active = resamp_active_mask(state);
if isfield(prof, 'G') && isstruct(prof.G)
    G = prof.G;
else
    G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
        'periodicX', true, 'periodicY', false, 'activeMask', active, 'cellWetMask', get_cell_wet_mask(state, params));
end

x = state.x(active,:);
m = state.m(active);
Nrel = G.N ./ max(params.gamma, eps) - 1;
Mrel = G.M ./ max(params.resampTargetCellMass, eps) - 1;
speed = sqrt(G.Ux.^2 + G.Uy.^2);

figure(opts.figureId); clf;
try
    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
catch
end
nexttile_or_subplot(2,3,1);
if ~isempty(x)
    scatter(x(:,1), x(:,2), opts.particleMarkerSize, m, 'filled');
else
    plot(nan,nan,'.');
end
axis equal tight; xlim([0 params.Lx]); ylim([0 params.Ly]); colorbar;
title('fluid particles colored by m_p'); xlabel('x'); ylabel('y');

nexttile_or_subplot(2,3,2);
imagesc([0 params.Lx], [0 params.Ly], Nrel.'); set(gca,'YDir','normal'); axis tight; colorbar; caxis_symmetric(Nrel, 0.5);
title('population N/\gamma - 1'); xlabel('x'); ylabel('y');

nexttile_or_subplot(2,3,3);
imagesc([0 params.Lx], [0 params.Ly], Mrel.'); set(gca,'YDir','normal'); axis tight; colorbar; caxis_symmetric(Mrel, 0.05);
title('cell mass M/M_{target} - 1'); xlabel('x'); ylabel('y');

nexttile_or_subplot(2,3,4);
imagesc([0 params.Lx], [0 params.Ly], G.Ux.'); set(gca,'YDir','normal'); axis tight; colorbar;
title('U_x field'); xlabel('x'); ylabel('y');

nexttile_or_subplot(2,3,5);
imagesc([0 params.Lx], [0 params.Ly], speed.'); set(gca,'YDir','normal'); axis tight; colorbar; hold on;
qstepX = max(1, round(params.Nx/16)); qstepY = max(1, round(params.Ny/12));
xc = ((0:params.Nx-1)+0.5) * params.Lx / params.Nx;
yc = ((0:params.Ny-1)+0.5) * params.Ly / params.Ny;
[XX,YY] = ndgrid(xc, yc);
quiver(XX(1:qstepX:end,1:qstepY:end), YY(1:qstepX:end,1:qstepY:end), ...
    G.Ux(1:qstepX:end,1:qstepY:end), G.Uy(1:qstepX:end,1:qstepY:end), 'k');
title('|u| + quiver'); xlabel('x'); ylabel('y'); hold off;

nexttile_or_subplot(2,3,6);
classMap = zeros(size(G.N));
classMap(G.N < get_param(params,'resampNMin',0)) = -1;
classMap(G.N > get_param(params,'resampNMax',inf)) = 1;
imagesc([0 params.Lx], [0 params.Ly], classMap.'); set(gca,'YDir','normal'); axis tight; colorbar; caxis([-1 1]);
title('population class: poor/ok/over'); xlabel('x'); ylabel('y');

sgtitle_safe(sprintf('%s step=%d t=%.3g kBT=%.4g Mrel=%.2g mRel=%.3g N=[%g,%g] center-wall=%.4g wall=%.4g VP=%.0f', ...
    opts.caseLabel, opts.step, opts.t, get_field(md,'kBTWeighted',NaN), get_field(md,'MRelRms',NaN), ...
    get_field(md,'mParticleRelStd',NaN), get_field(md,'NMin',NaN), get_field(md,'NMax',NaN), ...
    get_field(prof,'centerMinusWall',NaN), get_field(prof,'wallMeanVelocity',NaN), opts.nVirtualWallParticlesEquivalent));
drawnow limitrate;

if opts.showProfileFigure
    figure(opts.profileFigureId); clf; hold on;
    plot(prof.Ux, prof.yCenters, 'o-', 'DisplayName', 'mean U_x(y)');
    if isstruct(opts.fitInfo) && isfield(opts.fitInfo,'y')
        plot(opts.fitInfo.UxFitSlip, opts.fitInfo.y, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('slip fit: nu=%.4g R^2=%.3f', opts.fitInfo.nuEffSlip, opts.fitInfo.R2Slip));
        plot(opts.fitInfo.UxFitNoSlip, opts.fitInfo.y, '--', 'LineWidth', 1.2, 'DisplayName', sprintf('no-slip fit: nu=%.4g R^2=%.3f', opts.fitInfo.nuEffNoSlip, opts.fitInfo.R2NoSlip));
        if isfield(opts.fitInfo,'fitMask')
            yfit = opts.fitInfo.y(opts.fitInfo.fitMask);
            ufit = opts.fitInfo.UxMean(opts.fitInfo.fitMask);
            plot(ufit, yfit, 'ks', 'MarkerSize', 4, 'DisplayName', 'fit rows');
        end
        subtitle_str = sprintf('uSlip=%.4g slipRatio=%.3g center-wall=%.4g', opts.fitInfo.slipVelocity, opts.fitInfo.slipRatio, opts.fitInfo.centerMinusWall);
    else
        subtitle_str = sprintf('center-wall=%.4g wall=%.4g', get_field(prof,'centerMinusWall',NaN), get_field(prof,'wallMeanVelocity',NaN));
    end
    set(gca,'YDir','normal'); grid on; xlabel('U_x'); ylabel('y');
    title(sprintf('%s Poiseuille profile, step=%d, %s', opts.caseLabel, opts.step, subtitle_str), 'Interpreter','none');
    legend('Location','best'); hold off; drawnow limitrate;
end

if opts.saveFrame
    if ~exist(opts.frameDir,'dir'), mkdir(opts.frameDir); end
    spatialFile = fullfile(opts.frameDir, sprintf('%s_step%07d_spatial.png', sanitize(opts.caseLabel), opts.step));
    profileFile = fullfile(opts.frameDir, sprintf('%s_step%07d_profile.png', sanitize(opts.caseLabel), opts.step));
    try
        figure(opts.figureId); exportgraphics(gcf, spatialFile, 'Resolution', opts.pngResolution);
        if opts.showProfileFigure
            figure(opts.profileFigureId); exportgraphics(gcf, profileFile, 'Resolution', opts.pngResolution);
        end
    catch
        figure(opts.figureId); saveas(gcf, spatialFile);
        if opts.showProfileFigure
            figure(opts.profileFigureId); saveas(gcf, profileFile);
        end
    end
end
end

function nexttile_or_subplot(m,n,k)
try
    nexttile(k);
catch
    subplot(m,n,k);
end
end

function caxis_symmetric(A, fallback)
mx = max(abs(A(:)), [], 'omitnan');
if isempty(mx) || ~isfinite(mx) || mx == 0
    mx = fallback;
end
caxis([-mx mx]);
end

function mask = get_cell_wet_mask(state, params)
[mask,~] = resamp_cell_wet_mask(state, params, 'mode', 'auto');
end

function v = get_param(params, name, defaultValue)
if isfield(params,name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end

function v = get_field(s,name,defaultValue)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end

function sgtitle_safe(str)
try
    sgtitle(str, 'Interpreter','none');
catch
    title(str, 'Interpreter','none');
end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^a-zA-Z0-9_\-]', '_');
end

function opts = parse_options(varargin)
opts = struct();
opts.figureId = 720;
opts.profileFigureId = 721;
opts.showProfileFigure = true;
opts.caseLabel = 'poiseuille';
opts.step = 0;
opts.t = 0;
opts.fitInfo = [];
opts.saveFrame = false;
opts.frameDir = fullfile('runs','resamp_poiseuille_frames');
opts.particleMarkerSize = 3;
opts.nVirtualWallParticlesEquivalent = NaN;
opts.pngResolution = 150;
if mod(numel(varargin),2) ~= 0
    error('Options must be name/value pairs.');
end
for k=1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
        case 'figureid', opts.figureId = val;
        case 'profilefigureid', opts.profileFigureId = val;
        case 'showprofilefigure', opts.showProfileFigure = logical(val);
        case 'caselabel', opts.caseLabel = char(string(val));
        case 'step', opts.step = val;
        case 't', opts.t = val;
        case 'fitinfo', opts.fitInfo = val;
        case 'saveframe', opts.saveFrame = logical(val);
        case 'framedir', opts.frameDir = char(string(val));
        case 'particlemarkersize', opts.particleMarkerSize = val;
        case 'nvirtualwallparticlesequivalent', opts.nVirtualWallParticlesEquivalent = val;
        case 'pngresolution', opts.pngResolution = val;
        otherwise, error('Unknown option: %s', key);
    end
end
end
