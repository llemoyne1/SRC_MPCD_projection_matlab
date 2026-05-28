function info = resamp_pool_visualize_frame(state, params, step, varargin)
%RESAMP_POOL_VISUALIZE_FRAME Live visualization for weighted resampling pool runs.
%
%   info = resamp_pool_visualize_frame(state, params, step)
%
% The figure is deliberately resampling-oriented.  It separates particle
% population N, weighted cell mass M, particle masses m_p, velocity, vorticity
% and the poor/overpopulated/insertion map.  This is different from the
% historical Taylor-Green visualization, where N was also the density.
%
% Optional name-value arguments:
%   'insertDiag'         insertion diagnostic struct from resamp_insert_...
%   'remapDiag'          remap diagnostic struct
%   'stepDiag'           SRC/Q6 diagnostic struct
%   'summary'            optional row/struct/table with scalar diagnostics
%   'figureId'           default params.visualFigureId or 620
%   'figureName'         default 'Weighted resampling pool'
%   'showDebugFigure'    default false
%   'debugFigureId'      default figureId+1
%   'saveFrame'          default false
%   'frameDir'           default fullfile('frames','resamp_pool')
%   'framePrefix'        default 'resamp_pool'
%   'maxParticles'       default params.visualMaxParticles or 12000
%   'quiverStrideX/Y'    default automatic
%   'titleSuffix'        default ''
%
% The function does not modify the state.

if nargin < 3
    error('Usage: info = resamp_pool_visualize_frame(state, params, step, ...)');
end

opts = parse_options(params, varargin{:});
active = resamp_active_mask(state);
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', true, 'periodicY', true, 'activeMask', active);
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);

Nx = params.Nx;
Ny = params.Ny;
dx = params.Lx / Nx;
dy = params.Ly / Ny;
xc = ((0:Nx-1) + 0.5) * dx;
yc = ((0:Ny-1) + 0.5) * dy;
[X, Y] = ndgrid(xc, yc); %#ok<ASGLU>

gamma = getf(params, 'gamma', max(mean(G.N(:), 'omitnan'), eps));
m0 = getf(params, 'resampParticleMass', getf(params, 'particleMass', 1.0));
Mtarget = getf(params, 'resampTargetCellMass', gamma * m0);
Nmin = getf(params, 'resampNMin', ceil(0.5 * gamma));
Nmax = getf(params, 'resampNMax', ceil(1.5 * gamma));

Nrel = double(G.N) / max(gamma, eps) - 1.0;
Mrel = double(G.M) / max(Mtarget, eps) - 1.0;
speed = sqrt(G.Ux.^2 + G.Uy.^2);
omega = periodic_vorticity(G.Ux, G.Uy, dx, dy);
[classMap, classLabels] = population_class_map(G.N, Nmin, Nmax, opts.insertDiag);

fig = figure(opts.figureId); %#ok<NASGU>
clf;
set(gcf, 'Name', opts.figureName, 'Color', 'w');
method = getf(params, 'method', '');
if isempty(method)
    method = 'resamp';
end
sgtitle(sprintf(['%s step %d, t=%.4g %s | Nact=%d/%d free=%d insertedNow=%s cum=%s last=%s ', ...
                 '| N[min,max]=[%g,%g] MrelRMS=%.2e mRelStd=%.2e kBT=%.3g thermAfter=%.3g'], ...
    char(string(method)), step, step * getf(params,'dt',0), opts.titleSuffix, ...
    md.NpActive, md.Ncapacity, md.Nfree, ...
    format_scalar(getf(opts.insertDiag,'nInsertedParticles',NaN)), ...
    format_scalar(getf(opts.insertDiag,'nInsertedParticlesCumulative',NaN)), ...
    format_scalar(getf(opts.insertDiag,'lastInsertionStep',NaN)), ...
    md.NMin, md.NMax, md.MRelRms, md.mParticleRelStd, md.kBTWeighted, ...
    getf(getf(opts.stepDiag,'thermostatAfterRemap',struct()), 'meanKBTAfter', NaN)));

subplot(2,3,1);
plot_particles_by_mass(state, active, params, opts.maxParticles, m0);
title('active particles colored by m_p');

subplot(2,3,2);
imagesc(xc, yc, Nrel'); axis xy equal tight; colorbar;
title('population N/\gamma - 1'); xlabel('x'); ylabel('y');
set_symmetric_clim(gca, Nrel, getf(params,'visualPopulationRelCLim',NaN));

subplot(2,3,3);
imagesc(xc, yc, Mrel'); axis xy equal tight; colorbar;
title('cell mass M/M_{target} - 1'); xlabel('x'); ylabel('y');
set_symmetric_clim(gca, Mrel, getf(params,'visualMassRelCLim',NaN));

subplot(2,3,4);
imagesc(xc, yc, speed'); axis xy equal tight; hold on; colorbar;
strideX = opts.quiverStrideX;
strideY = opts.quiverStrideY;
ix = 1:strideX:Nx;
iy = 1:strideY:Ny;
quiver(X(ix,iy), Y(ix,iy), G.Ux(ix,iy), G.Uy(ix,iy), getf(params,'visualQuiverScale',1.2), 'k');
hold off; title('|u| + quiver'); xlabel('x'); ylabel('y');

subplot(2,3,5);
imagesc(xc, yc, omega'); axis xy equal tight; colorbar;
title('\omega = \partial_x u_y - \partial_y u_x'); xlabel('x'); ylabel('y');
set_symmetric_clim(gca, omega, getf(params,'visualOmegaCLim',NaN));

subplot(2,3,6);
imagesc(xc, yc, classMap'); axis xy equal tight;
caxis([-2 2]);
cbClass = colorbar;
set(cbClass, 'Ticks', -2:2, 'TickLabels', {'empty','poor','ok','over','insert'});
title(sprintf('population class: %s', classLabels)); xlabel('x'); ylabel('y');

info = struct();
info.G = G;
info.massDiagnostics = md;
info.omega = omega;
info.classMap = classMap;
info.Nrel = Nrel;
info.Mrel = Mrel;
info.speed = speed;
info.Nmin = Nmin;
info.Nmax = Nmax;
info.Mtarget = Mtarget;

drawnow;

if opts.showDebugFigure
    plot_debug_figure(state, active, params, G, md, opts, step, xc, yc, classMap);
end

if opts.saveFrame
    if ~exist(opts.frameDir, 'dir')
        mkdir(opts.frameDir);
    end
    exportgraphics(figure(opts.figureId), fullfile(opts.frameDir, sprintf('%s_%07d.png', opts.framePrefix, step)));
    if opts.showDebugFigure
        exportgraphics(figure(opts.debugFigureId), fullfile(opts.frameDir, sprintf('%s_debug_%07d.png', opts.framePrefix, step)));
    end
end
end

function opts = parse_options(params, varargin)
opts = struct();
opts.insertDiag = struct();
opts.remapDiag = struct();
opts.stepDiag = struct();
opts.summary = struct();
opts.figureId = getf(params, 'visualFigureId', 620);
opts.figureName = getf(params, 'visualFigureName', 'Weighted resampling pool');
opts.debugFigureId = opts.figureId + 1;
opts.showDebugFigure = false;
opts.saveFrame = false;
opts.frameDir = fullfile('frames', 'resamp_pool');
opts.framePrefix = 'resamp_pool';
opts.maxParticles = getf(params, 'visualMaxParticles', 12000);
opts.quiverStrideX = getf(params, 'visualQuiverStrideX', max(1, round(params.Nx / 18)));
opts.quiverStrideY = getf(params, 'visualQuiverStrideY', max(1, round(params.Ny / 12)));
opts.titleSuffix = getf(params, 'visualTitleSuffix', '');

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "insertdiag"
            opts.insertDiag = val;
        case "remapdiag"
            opts.remapDiag = val;
        case "stepdiag"
            opts.stepDiag = val;
        case "summary"
            opts.summary = val;
        case "figureid"
            opts.figureId = val;
        case "figurename"
            opts.figureName = char(string(val));
        case "debugfigureid"
            opts.debugFigureId = val;
        case "showdebugfigure"
            opts.showDebugFigure = logical(val);
        case "saveframe"
            opts.saveFrame = logical(val);
        case "framedir"
            opts.frameDir = char(string(val));
        case "frameprefix"
            opts.framePrefix = char(string(val));
        case "maxparticles"
            opts.maxParticles = val;
        case "quiverstridex"
            opts.quiverStrideX = val;
        case "quiverstridey"
            opts.quiverStrideY = val;
        case "titlesuffix"
            opts.titleSuffix = char(string(val));
        otherwise
            error('Unknown option: %s', string(key));
    end
end
end

function plot_particles_by_mass(state, active, params, maxParticles, m0)
idx = find(active);
if isempty(idx)
    cla; axis equal tight; xlim([0 params.Lx]); ylim([0 params.Ly]);
    title('no active particles'); xlabel('x'); ylabel('y');
    return;
end
if numel(idx) > maxParticles
    idx = idx(randperm(numel(idx), maxParticles));
end
m = state.m(idx);
scatter(state.x(idx,1), state.x(idx,2), 5, m ./ max(m0, eps), 'filled');
axis equal tight; xlim([0 params.Lx]); ylim([0 params.Ly]); colorbar;
xlabel('x'); ylabel('y');
end

function omega = periodic_vorticity(Ux, Uy, dx, dy)
dUyDx = (circshift(Uy, [-1, 0]) - circshift(Uy, [1, 0])) / (2 * dx);
dUxDy = (circshift(Ux, [0, -1]) - circshift(Ux, [0, 1])) / (2 * dy);
omega = dUyDx - dUxDy;
omega(~isfinite(omega)) = 0;
end

function [classMap, label] = population_class_map(N, Nmin, Nmax, insertDiag)
% Codes: -2 empty, -1 poor, 0 nominal, +1 over, +2 inserted this step.
classMap = zeros(size(N));
classMap(N == 0) = -2;
classMap(N > 0 & N < Nmin) = -1;
classMap(N > Nmax) = 1;
if isstruct(insertDiag) && isfield(insertDiag, 'insertedPerCellGrid') && ~isempty(insertDiag.insertedPerCellGrid)
    ins = insertDiag.insertedPerCellGrid;
    if isequal(size(ins), size(N))
        classMap(ins > 0) = 2;
    end
end
label = '-2 empty, -1 poor, 0 ok, 1 over, 2 inserted';
end

function plot_debug_figure(state, active, params, G, md, opts, step, xc, yc, classMap) %#ok<INUSD>
Nx = params.Nx;
Ny = params.Ny;
gamma = getf(params, 'gamma', max(mean(G.N(:), 'omitnan'), eps));
m0 = getf(params, 'resampParticleMass', getf(params, 'particleMass', 1.0));
[meanMassGrid, relStdMassGrid] = cell_mass_stats(state, active, G.cellId, Nx, Ny, m0);
insertGrid = zeros(Nx, Ny);
if isstruct(opts.insertDiag) && isfield(opts.insertDiag, 'insertedPerCellGrid') && ~isempty(opts.insertDiag.insertedPerCellGrid)
    if isequal(size(opts.insertDiag.insertedPerCellGrid), [Nx, Ny])
        insertGrid = opts.insertDiag.insertedPerCellGrid;
    end
end
if isfield(state, 'uMemUx') && isequal(size(state.uMemUx), [Nx, Ny])
    memDiff = sqrt((G.Ux - state.uMemUx).^2 + (G.Uy - state.uMemUy).^2);
    if isfield(state, 'uMemValid') && isequal(size(state.uMemValid), [Nx, Ny])
        memDiff(~state.uMemValid) = NaN;
    end
else
    memDiff = NaN(Nx, Ny);
end

figure(opts.debugFigureId); clf;
set(gcf, 'Name', [opts.figureName ' debug'], 'Color', 'w');
sgtitle(sprintf('pool/resampling debug step %d | active %.3f, free %d, m=[%.3g %.3g], class values shown in main fig', ...
    step, md.activeFraction, md.Nfree, md.mParticleMin, md.mParticleMax));

subplot(2,3,1);
imagesc(xc, yc, meanMassGrid'); axis xy equal tight; colorbar;
title('mean particle mass per cell / m_0'); xlabel('x'); ylabel('y');
set_symmetric_clim(gca, meanMassGrid - 1, NaN);

subplot(2,3,2);
imagesc(xc, yc, relStdMassGrid'); axis xy equal tight; colorbar;
title('relative std(m_p) per cell'); xlabel('x'); ylabel('y');

subplot(2,3,3);
imagesc(xc, yc, insertGrid'); axis xy equal tight; colorbar;
title('particles inserted at current step'); xlabel('x'); ylabel('y');

subplot(2,3,4);
imagesc(xc, yc, memDiff'); axis xy equal tight; colorbar;
title('|u - u_{mem}|'); xlabel('x'); ylabel('y');

subplot(2,3,5);
m = state.m(active);
if isempty(m)
    histogram(0);
else
    histogram(m ./ max(m0, eps), max(10, min(60, round(sqrt(numel(m))))));
end
grid on; xlabel('m_p / m_0'); ylabel('count'); title('active particle mass histogram');

subplot(2,3,6);
bar([md.NpActive, md.Nfree, md.Ncapacity, nnz(G.N(:) < ceil(0.5*gamma)), nnz(G.N(:) > ceil(1.5*gamma))]);
set(gca, 'XTickLabel', {'active','free','cap','poor','over'});
grid on; title('pool/population counters'); ylabel('count');

drawnow;
end

function [meanMassGrid, relStdMassGrid] = cell_mass_stats(state, active, cellIdAll, Nx, Ny, m0)
Nc = Nx * Ny;
meanMass = NaN(Nc, 1);
relStdMass = NaN(Nc, 1);
idsActive = find(active);
if ~isempty(idsActive)
    cActive = cellIdAll(idsActive);
    mActive = state.m(idsActive) ./ max(m0, eps);
    valid = cActive > 0 & isfinite(mActive);
    if any(valid)
        meanMassTmp = accumarray(cActive(valid), mActive(valid), [Nc, 1], @mean, NaN);
        stdMassTmp = accumarray(cActive(valid), mActive(valid), [Nc, 1], @std0, NaN);
        meanMass = meanMassTmp;
        relStdMass = stdMassTmp ./ max(abs(meanMassTmp), eps);
    end
end
meanMassGrid = reshape(meanMass, [Ny, Nx]).';
relStdMassGrid = reshape(relStdMass, [Ny, Nx]).';
end

function s = std0(x)
if numel(x) <= 1
    s = 0;
else
    s = std(x, 0, 'omitnan');
end
end

function set_symmetric_clim(ax, data, requested)
if isfinite(requested) && requested > 0
    caxis(ax, [-requested requested]);
    return;
end
vals = data(:);
vals = vals(isfinite(vals));
if isempty(vals)
    return;
end
c = prctile(abs(vals), 98);
if isfinite(c) && c > 0
    caxis(ax, [-c c]);
end
end

function s = format_scalar(v)
if isempty(v) || ~isfinite(v)
    s = 'n/a';
else
    s = sprintf('%g', v);
end
end

function v = getf(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
