function projection_step_replay_fields(out, varargin)
%PROJECTION_STEP_REPLAY_FIELDS Replay stored density/velocity/vorticity maps.
%
%   projection_step_replay_fields(out)
%   projection_step_replay_fields(out, 'pauseTime', 0.03, 'stride', 1)
%
% Requires out.NMaps, out.UxMaps, out.UyMaps and out.omegaMaps.  Particle
% trajectories are not replayed because the benchmark stores only the final
% particle state; use params.visualEnable=true for live particle animation.

opts.pauseTime = 0.03;
opts.stride = 1;
opts.figureId = 120;
opts.densityCLim = 0.5;
opts.omegaCLim = NaN;
opts.quiverStrideX = 4;
opts.quiverStrideY = 3;
opts.quiverScale = 1.5;
opts = parse_opts(opts, varargin{:});

if ~isfield(out, 'NMaps') || isempty(out.NMaps)
    error('out.NMaps is missing. Run with params.storeDensityMaps=true.');
end
if ~isfield(out, 'omegaMaps') || isempty(out.omegaMaps)
    error('out.omegaMaps is missing. Run with params.storeVorticityMaps=true.');
end
if ~isfield(out, 'UxMaps') || isempty(out.UxMaps) || ~isfield(out, 'UyMaps') || isempty(out.UyMaps)
    error('out.UxMaps/out.UyMaps are missing. Run with params.storeVelocityMaps=true.');
end

params = out.params;
solidMask = false(params.Nx, params.Ny);
if isfield(out, 'solidMask') && ~isempty(out.solidMask)
    solidMask = logical(out.solidMask);
end
xc = ((1:params.Nx) - 0.5) * (params.Lx / params.Nx);
yc = ((1:params.Ny) - 0.5) * (params.Ly / params.Ny);
[Xc, Yc] = ndgrid(xc, yc);
ix = 1:max(1, round(opts.quiverStrideX)):params.Nx;
iy = 1:max(1, round(opts.quiverStrideY)):params.Ny;

fig = figure(opts.figureId);
set(fig, 'Name', 'Step-channel field replay', 'Color', 'w');
for k = 1:max(1, round(opts.stride)):numel(out.sampleTimes)
    N = out.NMaps(:,:,k);
    Ux = out.UxMaps(:,:,k);
    Uy = out.UyMaps(:,:,k);
    omega = out.omegaMaps(:,:,k);
    relN = N ./ max(params.gamma, eps) - 1.0;
    relN(solidMask) = NaN;
    Ux(solidMask) = NaN;
    Uy(solidMask) = NaN;
    omega(solidMask) = NaN;
    speed = sqrt(Ux.^2 + Uy.^2);

    clf(fig);
    subplot(1,3,1);
    imagesc(xc, yc, relN.'); set(gca, 'YDir', 'normal'); axis equal tight; colorbar;
    caxis([-opts.densityCLim opts.densityCLim]); hold on; draw_step_patch(params);
    title(sprintf('N/\\gamma - 1, t=%.4g', out.sampleTimes(k)));

    subplot(1,3,2);
    imagesc(xc, yc, speed.'); set(gca, 'YDir', 'normal'); axis equal tight; colorbar; hold on;
    quiver(Xc(ix,iy), Yc(ix,iy), Ux(ix,iy), Uy(ix,iy), opts.quiverScale, 'k');
    draw_step_patch(params); title('|u| + quiver');

    subplot(1,3,3);
    imagesc(xc, yc, omega.'); set(gca, 'YDir', 'normal'); axis equal tight; colorbar; hold on;
    omegaCLim = opts.omegaCLim;
    if ~(isfinite(omegaCLim) && omegaCLim > 0)
        vals = omega(isfinite(omega));
        omegaCLim = max(1e-12, 3 * std(vals));
    end
    caxis([-omegaCLim omegaCLim]); draw_step_patch(params); title('\omega');
    drawnow;
    if opts.pauseTime > 0
        pause(opts.pauseTime);
    end
end
end

function opts = parse_opts(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for i = 1:2:numel(varargin)
    name = varargin{i};
    value = varargin{i+1};
    if ~isfield(opts, name)
        error('Unknown option: %s', name);
    end
    opts.(name) = value;
end
end

function draw_step_patch(params)
x0 = get_param(params, 'stepX0', 0.15 * params.Lx);
x1 = get_param(params, 'stepX1', 0.35 * params.Lx);
h = get_param(params, 'stepHeight', 0.25 * params.Ly);
patch([x0 x1 x1 x0], [0 0 h h], [0.1 0.1 0.1], ...
    'FaceAlpha', 0.35, 'EdgeColor', 'k', 'LineWidth', 1.0);
end

function value = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
