function [vOut, info] = mpcd_srd_collision_channel_virtual_walls(x, v, params)
%MPCD_SRD_COLLISION_CHANNEL_VIRTUAL_WALLS SRD collision in an x-periodic channel.
%
%   [vOut, info] = mpcd_srd_collision_channel_virtual_walls(x, v, params)
%
% Compact MATLAB reference kernel for channel-wall SRC/MPCD collisions.
% Real particles are binned once, optional virtual wall particles contribute
% only aggregate collision mass/momentum, and only real particles are
% rotated/updated.  This is designed to map naturally to OpenMP/GPU kernels:
% bin particles -> reduce cell aggregates -> add wall-face aggregates ->
% rotate real particles.
%
% Main wall virtual particle modes:
%   shifted_solid_fraction  physically preferred.  The y collision grid is
%                           randomly shifted; cells intersecting y<0 or y>Ly
%                           receive virtual wall particles in proportion to
%                           their solid fraction.  A small extended grid in y
%                           is used so particles in shifted ghost cells collide
%                           consistently.
%   fixed_boundary_layers   v1 fallback: add a fixed count in the first/last
%                           physical cell layers.
%   none/off                no virtual particles.
%
% Important parameters:
%   wallVirtualParticlesEnable              default false
%   wallVirtualParticlesGeometryMode        default 'shifted_solid_fraction'
%   wallVirtualParticlesForceRandomShiftY   default true when wall VP enabled
%   wallVirtualParticlesDensityFactor       default 1.0
%   wallVirtualParticlesPerCell             optional explicit count per full cell
%   wallVirtualParticlesThermal             default true
%   wallVirtualParticlesKBT                 default params.kBT
%   wallVirtualParticlesStochasticCount     stochastic rounding of fractional counts

validate_state_arrays(x, v);

Lx = get_param(params, 'Lx', 1.0);
Ly = get_param(params, 'Ly', 1.0);
Nx = get_param(params, 'Nx', 1);
Ny = get_param(params, 'Ny', 1);
alphaDeg = get_param(params, 'alphaDeg', 90);

if Lx <= 0 || Ly <= 0 || Nx <= 0 || Ny <= 0
    error('Invalid grid/domain parameters.');
end

dx = Lx / Nx;
dy = Ly / Ny;

enabled = logical(get_param(params, 'wallVirtualParticlesEnable', false));
geometryMode = lower(char(string(get_param(params, 'wallVirtualParticlesGeometryMode', 'shifted_solid_fraction'))));
if any(strcmp(geometryMode, {'off','none','false'}))
    enabled = false;
end
if ~enabled
    geometryMode = 'none';
end

useRandomGridShiftX = logical(get_param(params, 'useRandomGridShiftX', true));
useRandomGridShiftY = logical(get_param(params, 'useRandomGridShiftY', false));
if enabled && strcmp(geometryMode, 'shifted_solid_fraction') && ...
        logical(get_param(params, 'wallVirtualParticlesForceRandomShiftY', true))
    useRandomGridShiftY = true;
end

if useRandomGridShiftX
    shiftX = (rand() - 0.5) * dx;
else
    shiftX = 0.0;
end
if useRandomGridShiftY
    shiftY = (rand() - 0.5) * dy;
else
    shiftY = 0.0;
end

xs = mod(x(:, 1) + shiftX, Lx);
ix = floor(xs / dx) + 1;
ix = min(max(ix, 1), Nx);

switch geometryMode
    case 'shifted_solid_fraction'
        [cellId, Ncy, iyRaw] = assign_y_cells_extended(x(:,2), shiftY, dy, Ny);
        cellId = cellId + Ncy * (ix - 1);
        Nc = Nx * Ncy;
        [Nvirt, PxVirt, PyVirt, virtInfo] = virtual_wall_moments_shifted_solid_fraction(params, Nx, Ny, Ncy, Nc, dy, Ly, shiftY);
    case {'fixed_boundary_layers','fixed_layers','boundary_layers'}
        ys = min(max(x(:, 2) + shiftY, 0), Ly - eps(Ly));
        iy = floor(ys / dy) + 1;
        iy = min(max(iy, 1), Ny);
        Ncy = Ny;
        iyRaw = iy;
        cellId = iy + Ny * (ix - 1);
        Nc = Nx * Ny;
        [Nvirt, PxVirt, PyVirt, virtInfo] = virtual_wall_moments_fixed_layers(params, Nx, Ny, Nc);
    otherwise
        if enabled
            error('Unknown wallVirtualParticlesGeometryMode: %s', geometryMode);
        end
        ys = min(max(x(:, 2) + shiftY, 0), Ly - eps(Ly));
        iy = floor(ys / dy) + 1;
        iy = min(max(iy, 1), Ny);
        Ncy = Ny;
        iyRaw = iy;
        cellId = iy + Ny * (ix - 1);
        Nc = Nx * Ny;
        Nvirt = zeros(Nc, 1);
        PxVirt = zeros(Nc, 1);
        PyVirt = zeros(Nc, 1);
        virtInfo = empty_virtual_info(enabled, geometryMode);
end

Nreal = accumarray(cellId, 1, [Nc, 1], @sum, 0);
PxReal = accumarray(cellId, v(:, 1), [Nc, 1], @sum, 0);
PyReal = accumarray(cellId, v(:, 2), [Nc, 1], @sum, 0);

Ntot = double(Nreal) + double(Nvirt);
PxTot = PxReal + PxVirt;
PyTot = PyReal + PyVirt;

Ux = zeros(Nc, 1);
Uy = zeros(Nc, 1);
occTot = Ntot > 0;
Ux(occTot) = PxTot(occTot) ./ Ntot(occTot);
Uy(occTot) = PyTot(occTot) ./ Ntot(occTot);

alpha = alphaDeg * pi / 180.0;
signs = 2.0 * (rand(Nc, 1) > 0.5) - 1.0;
angles = signs * alpha;
ca = cos(angles(cellId));
sa = sin(angles(cellId));

uxp = Ux(cellId);
uyp = Uy(cellId);
rvx = v(:, 1) - uxp;
rvy = v(:, 2) - uyp;

vOut = v;
vOut(:, 1) = uxp + ca .* rvx - sa .* rvy;
vOut(:, 2) = uyp + sa .* rvx + ca .* rvy;

info = struct();
info.Nx = Nx;
info.Ny = Ny;
info.NcyCollision = Ncy;
info.Nc = Nc;
info.dx = dx;
info.dy = dy;
info.shiftX = shiftX;
info.shiftY = shiftY;
info.useRandomGridShiftX = useRandomGridShiftX;
info.useRandomGridShiftY = useRandomGridShiftY;
info.geometryMode = geometryMode;
info.cellId = cellId;
info.iyRaw = iyRaw;
info.NcellReal = Nreal;
info.NcellVirtual = Nvirt;
info.NcellTotal = Ntot;
info.PxReal = PxReal;
info.PyReal = PyReal;
info.PxVirtual = PxVirt;
info.PyVirtual = PyVirt;
info.alphaDeg = alphaDeg;
info.wallVirtualParticles = virtInfo;
info.NMeanReal = mean(double(Nreal), 'omitnan');
info.NStdReal = std(double(Nreal), 0, 'omitnan');
info.NMinReal = min(Nreal);
info.NMaxReal = max(Nreal);
info.nEmptyRealCells = nnz(Nreal == 0);
info.NMeanTotal = mean(double(Ntot), 'omitnan');
info.NStdTotal = std(double(Ntot), 0, 'omitnan');
info.nVirtualCells = nnz(Nvirt > 0);
info.nVirtualParticlesTotal = sum(Nvirt, 'omitnan');
info.virtualMomentumXTotal = sum(PxVirt, 'omitnan');
info.virtualMomentumYTotal = sum(PyVirt, 'omitnan');
end

function [iyExt, Ncy, iyRaw] = assign_y_cells_extended(y, shiftY, dy, Ny)
% Raw shifted collision-cell index is allowed to be 0..Ny+1.  The +1 shift
% maps it to MATLAB indices 1..Ny+2.  This keeps particles in shifted ghost
% cells in a real collision cell instead of clamping them into the wrong
% boundary layer.
iyRaw = floor((y + shiftY) / dy) + 1;
iyRaw = min(max(iyRaw, 0), Ny + 1);
Ncy = Ny + 2;
iyExt = iyRaw + 1;
end

function [Nvirt, PxVirt, PyVirt, info] = virtual_wall_moments_shifted_solid_fraction(params, Nx, Ny, Ncy, Nc, dy, Ly, shiftY)
Nvirt = zeros(Nc, 1);
PxVirt = zeros(Nc, 1);
PyVirt = zeros(Nc, 1);

enabled = logical(get_param(params, 'wallVirtualParticlesEnable', false));
info = empty_virtual_info(enabled, 'shifted_solid_fraction');
info.shiftY = shiftY;
if ~enabled
    return;
end

includeBottom = logical(get_param(params, 'wallVirtualParticlesIncludeBottom', true));
includeTop = logical(get_param(params, 'wallVirtualParticlesIncludeTop', true));
gamma = get_param(params, 'gamma', 0);
densityFactor = get_param(params, 'wallVirtualParticlesDensityFactor', 1.0);
countParam = get_param(params, 'wallVirtualParticlesPerCell', []);
if isempty(countParam)
    fullCellCount = gamma * densityFactor;
else
    fullCellCount = double(countParam);
end
fullCellCount = max(0, fullCellCount);
thermal = logical(get_param(params, 'wallVirtualParticlesThermal', true));
kBT = get_param(params, 'wallVirtualParticlesKBT', get_param(params, 'kBT', 1.0));
sigma = sqrt(max(kBT, 0));
stochasticCount = logical(get_param(params, 'wallVirtualParticlesStochasticCount', true));
Ubottom = get_param(params, 'Ubottom', 0.0);
Utop = get_param(params, 'Utop', 0.0);
Vbottom = get_param(params, 'Vbottom', 0.0);
Vtop = get_param(params, 'Vtop', 0.0);

info.fullCellCount = fullCellCount;
info.densityFactor = densityFactor;
info.thermal = thermal;
info.kBT = kBT;
info.stochasticCount = stochasticCount;

for ix = 1:Nx
    for jyRaw = 0:(Ny+1)
        cellStart = (jyRaw - 1) * dy - shiftY;
        cellEnd = jyRaw * dy - shiftY;
        bottomSolidLength = max(0.0, min(cellEnd, 0.0) - cellStart);
        topSolidLength = max(0.0, cellEnd - max(cellStart, Ly));
        bottomFrac = min(max(bottomSolidLength / dy, 0.0), 1.0);
        topFrac = min(max(topSolidLength / dy, 0.0), 1.0);
        id = (jyRaw + 1) + Ncy * (ix - 1);
        if includeBottom && bottomFrac > 0
            nAdd = sample_virtual_count(fullCellCount * bottomFrac, stochasticCount);
            [Nvirt, PxVirt, PyVirt] = add_virtual_to_cell(Nvirt, PxVirt, PyVirt, id, nAdd, Ubottom, Vbottom, sigma, thermal);
            info.bottomSolidVolumeCells = info.bottomSolidVolumeCells + bottomFrac;
            info.nBottomVirtualParticles = info.nBottomVirtualParticles + nAdd;
            if nAdd > 0, info.nBottomCells = info.nBottomCells + 1; end
        end
        if includeTop && topFrac > 0
            nAdd = sample_virtual_count(fullCellCount * topFrac, stochasticCount);
            [Nvirt, PxVirt, PyVirt] = add_virtual_to_cell(Nvirt, PxVirt, PyVirt, id, nAdd, Utop, Vtop, sigma, thermal);
            info.topSolidVolumeCells = info.topSolidVolumeCells + topFrac;
            info.nTopVirtualParticles = info.nTopVirtualParticles + nAdd;
            if nAdd > 0, info.nTopCells = info.nTopCells + 1; end
        end
    end
end

info.nVirtualTotal = sum(Nvirt, 'omitnan');
info.PxVirtualTotal = sum(PxVirt, 'omitnan');
info.PyVirtualTotal = sum(PyVirt, 'omitnan');
end

function [Nvirt, PxVirt, PyVirt, info] = virtual_wall_moments_fixed_layers(params, Nx, Ny, Nc)
Nvirt = zeros(Nc, 1);
PxVirt = zeros(Nc, 1);
PyVirt = zeros(Nc, 1);

enabled = logical(get_param(params, 'wallVirtualParticlesEnable', false));
info = empty_virtual_info(enabled, 'fixed_boundary_layers');
if ~enabled
    return;
end

nLayers = max(0, round(get_param(params, 'wallVirtualParticlesCells', 1)));
gamma = get_param(params, 'gamma', 0);
densityFactor = get_param(params, 'wallVirtualParticlesDensityFactor', 1.0);
countParam = get_param(params, 'wallVirtualParticlesPerCell', []);
if isempty(countParam)
    nPerCell = round(gamma * densityFactor);
else
    nPerCell = round(countParam);
end
nPerCell = max(0, nPerCell);
includeBottom = logical(get_param(params, 'wallVirtualParticlesIncludeBottom', true));
includeTop = logical(get_param(params, 'wallVirtualParticlesIncludeTop', true));
thermal = logical(get_param(params, 'wallVirtualParticlesThermal', true));
kBT = get_param(params, 'wallVirtualParticlesKBT', get_param(params, 'kBT', 1.0));
sigma = sqrt(max(kBT, 0));
Ubottom = get_param(params, 'Ubottom', 0.0);
Utop = get_param(params, 'Utop', 0.0);
Vbottom = get_param(params, 'Vbottom', 0.0);
Vtop = get_param(params, 'Vtop', 0.0);

info.nLayers = nLayers;
info.fullCellCount = nPerCell;
info.densityFactor = densityFactor;
info.thermal = thermal;
info.kBT = kBT;

if nLayers <= 0 || nPerCell <= 0
    return;
end
nLayers = min(nLayers, Ny);
bottomLayers = 1:nLayers;
topLayers = max(1, Ny-nLayers+1):Ny;
if includeBottom
    ids = wall_cell_ids(Nx, Ny, bottomLayers);
    for k = 1:numel(ids)
        [Nvirt, PxVirt, PyVirt] = add_virtual_to_cell(Nvirt, PxVirt, PyVirt, ids(k), nPerCell, Ubottom, Vbottom, sigma, thermal);
    end
    info.nBottomCells = info.nBottomCells + numel(ids);
    info.nBottomVirtualParticles = info.nBottomVirtualParticles + nPerCell * numel(ids);
end
if includeTop
    ids = wall_cell_ids(Nx, Ny, topLayers);
    for k = 1:numel(ids)
        [Nvirt, PxVirt, PyVirt] = add_virtual_to_cell(Nvirt, PxVirt, PyVirt, ids(k), nPerCell, Utop, Vtop, sigma, thermal);
    end
    info.nTopCells = info.nTopCells + numel(ids);
    info.nTopVirtualParticles = info.nTopVirtualParticles + nPerCell * numel(ids);
end
info.nVirtualTotal = sum(Nvirt, 'omitnan');
info.PxVirtualTotal = sum(PxVirt, 'omitnan');
info.PyVirtualTotal = sum(PyVirt, 'omitnan');
end

function n = sample_virtual_count(lambda, stochasticCount)
if lambda <= 0
    n = 0;
elseif stochasticCount
    base = floor(lambda);
    n = base + double(rand() < (lambda - base));
else
    n = round(lambda);
end
end

function [Nvirt, PxVirt, PyVirt] = add_virtual_to_cell(Nvirt, PxVirt, PyVirt, id, nAdd, Uwall, Vwall, sigma, thermal)
if nAdd <= 0
    return;
end
Nvirt(id) = Nvirt(id) + nAdd;
if thermal
    PxAdd = nAdd * Uwall + sqrt(nAdd) * sigma * randn();
    PyAdd = nAdd * Vwall + sqrt(nAdd) * sigma * randn();
else
    PxAdd = nAdd * Uwall;
    PyAdd = nAdd * Vwall;
end
PxVirt(id) = PxVirt(id) + PxAdd;
PyVirt(id) = PyVirt(id) + PyAdd;
end

function ids = wall_cell_ids(Nx, Ny, layers)
ids = [];
for ix = 1:Nx
    ids = [ids; layers(:) + Ny * (ix - 1)]; %#ok<AGROW>
end
end

function info = empty_virtual_info(enabled, geometryMode)
info = struct();
info.enabled = enabled;
info.geometryMode = geometryMode;
info.shiftY = NaN;
info.nLayers = NaN;
info.fullCellCount = NaN;
info.densityFactor = NaN;
info.thermal = NaN;
info.kBT = NaN;
info.stochasticCount = NaN;
info.nBottomCells = 0;
info.nTopCells = 0;
info.nBottomVirtualParticles = 0;
info.nTopVirtualParticles = 0;
info.bottomSolidVolumeCells = 0.0;
info.topSolidVolumeCells = 0.0;
info.nVirtualTotal = 0;
info.PxVirtualTotal = 0.0;
info.PyVirtualTotal = 0.0;
end

function validate_state_arrays(x, v)
if size(x, 2) ~= 2 || size(v, 2) ~= 2 || size(x, 1) ~= size(v, 1)
    error('x and v must be Np-by-2 arrays with matching particle count.');
end
end

function value = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
