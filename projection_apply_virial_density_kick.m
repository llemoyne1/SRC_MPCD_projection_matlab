function [vOut, info] = projection_apply_virial_density_kick(x, v, params, varargin)
%PROJECTION_APPLY_VIRIAL_DENSITY_KICK Local liquid-like density correction.
%
%   [vOut, info] = projection_apply_virial_density_kick(x, v, params)
%
% Applies a weak local virial-pressure kick after the pressure projection:
%
%      Pvir = Kvir * (N - gamma)
%      dv   = - strength * dt * grad(Pvir) / max(N, Nmin)
%
% The correction is purely velocity-based; particle positions and therefore the
% instantaneous cell population are not changed. The operation is deliberately
% local on the grid: particle deposit, nearest-neighbor stencil for grad(Pvir),
% and grid-to-particle interpolation. This keeps the method compatible with a
% future OpenMP/MPI/GPU implementation.
%
% Optional name-value arguments:
%   'periodicX'   default true
%   'periodicY'   default false
%   'method'      interpolation method, default params.virialInterpolationMethod
%
% Main parameters:
%   useVirialDensityKick        true/false gate
%   virialDensityKickStrength   dimensionless strength betaVir
%   virialK                     virial pressure coefficient
%   virialSmoothPasses          number of local smoothing passes on Pvir
%   virialMinCellCount          denominator floor for N
%   virialMaxParticleKick       optional per-particle |dv| limiter

if nargin < 3
    error('Usage: [vOut, info] = projection_apply_virial_density_kick(x, v, params)');
end
if size(x, 1) ~= size(v, 1) || size(x, 2) < 2 || size(v, 2) < 2
    error('x and v must be Np-by-2 arrays with the same number of particles.');
end

periodicX = true;
periodicY = false;
method = char(string(get_param(params, 'virialInterpolationMethod', ...
    get_param(params, 'projectionInterpolationMethod', 'nearest'))));
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "method"
            method = char(string(val));
        otherwise
            error('Unknown option: %s', string(key));
    end
end

useKick = logical(get_param(params, 'useVirialDensityKick', false));
strength = get_param(params, 'virialDensityKickStrength', 0.0);
Kvir = get_param(params, 'virialK', 1.0);
smoothPasses = round(get_param(params, 'virialSmoothPasses', 0));
minCellCount = get_param(params, 'virialMinCellCount', 1.0);
maxParticleKick = get_param(params, 'virialMaxParticleKick', Inf);

vOut = v;
baseInfo = empty_info(useKick, strength, Kvir, smoothPasses, minCellCount, maxParticleKick);
if ~useKick || strength == 0 || Kvir == 0
    info = baseInfo;
    return;
end

G = projection_deposit_particles_to_grid(x, v, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minCount', 1);
N = double(G.N);
gamma = get_param(params, 'gamma', mean(N(:)));
if gamma <= 0
    error('params.gamma must be positive.');
end

Pvir = Kvir * (N - gamma);
for ipass = 1:max(0, smoothPasses)
    Pvir = smooth_cell_field(Pvir, periodicX, periodicY);
end

[gradPx, gradPy] = gradient_cell_centered(Pvir, params, periodicX, periodicY);
Nden = max(N, minCellCount);
dUx = -strength * params.dt * gradPx ./ Nden;
dUy = -strength * params.dt * gradPy ./ Nden;

dv = projection_interpolate_grid_delta_to_particles(x, dUx, dUy, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'method', method);

if isfinite(maxParticleKick) && maxParticleKick > 0
    mag = sqrt(sum(dv.^2, 2));
    tooLarge = mag > maxParticleKick;
    if any(tooLarge)
        scale = maxParticleKick ./ max(mag(tooLarge), eps);
        dv(tooLarge, :) = dv(tooLarge, :) .* scale;
    end
elseif maxParticleKick <= 0
    dv(:) = 0;
end

vOut = v + dv;

info = baseInfo;
info.enabled = true;
info.method = method;
info.meanN = mean(N(:));
info.stdN = std(N(:));
info.minN = min(N(:));
info.maxN = max(N(:));
info.pVirMean = mean(Pvir(:));
info.pVirRms = sqrt(mean(Pvir(:).^2));
info.pVirMaxAbs = max(abs(Pvir(:)));
info.gradPvirRms = sqrt(mean(gradPx(:).^2 + gradPy(:).^2));
info.gradPvirMaxAbs = max(sqrt(gradPx(:).^2 + gradPy(:).^2));
info.gridDVRms = sqrt(mean(dUx(:).^2 + dUy(:).^2));
info.gridDVMaxAbs = max(sqrt(dUx(:).^2 + dUy(:).^2));
info.particleDVRms = sqrt(mean(sum(dv.^2, 2)));
info.particleDVMaxAbs = max(sqrt(sum(dv.^2, 2)));
info.particleDVMeanX = mean(dv(:, 1));
info.particleDVMeanY = mean(dv(:, 2));
info.nLimitedParticles = nnz(sqrt(sum(dv.^2, 2)) >= maxParticleKick & isfinite(maxParticleKick) & maxParticleKick > 0);
info.N = N;
info.Pvir = Pvir;
info.gradPvirX = gradPx;
info.gradPvirY = gradPy;
info.dUx = dUx;
info.dUy = dUy;
end

function info = empty_info(useKick, strength, Kvir, smoothPasses, minCellCount, maxParticleKick)
info = struct();
info.enabled = false;
info.requested = logical(useKick);
info.strength = strength;
info.Kvir = Kvir;
info.smoothPasses = smoothPasses;
info.minCellCount = minCellCount;
info.maxParticleKick = maxParticleKick;
info.method = '';
info.meanN = NaN;
info.stdN = NaN;
info.minN = NaN;
info.maxN = NaN;
info.pVirMean = 0.0;
info.pVirRms = 0.0;
info.pVirMaxAbs = 0.0;
info.gradPvirRms = 0.0;
info.gradPvirMaxAbs = 0.0;
info.gridDVRms = 0.0;
info.gridDVMaxAbs = 0.0;
info.particleDVRms = 0.0;
info.particleDVMaxAbs = 0.0;
info.particleDVMeanX = 0.0;
info.particleDVMeanY = 0.0;
info.nLimitedParticles = 0;
info.N = [];
info.Pvir = [];
info.gradPvirX = [];
info.gradPvirY = [];
info.dUx = [];
info.dUy = [];
end

function F = smooth_cell_field(F, periodicX, periodicY)
[Nx, Ny] = size(F);
ixp = [2:Nx, 1];
ixm = [Nx, 1:Nx-1];
if ~periodicX
    ixp(end) = Nx;
    ixm(1) = 1;
end

Fxp = F(ixp, :);
Fxm = F(ixm, :);

Fyp = F;
Fym = F;
if Ny > 1
    Fyp(:, 1:Ny-1) = F(:, 2:Ny);
    Fym(:, 2:Ny) = F(:, 1:Ny-1);
    if periodicY
        Fyp(:, Ny) = F(:, 1);
        Fym(:, 1) = F(:, Ny);
    end
end

% Conservative compact smoothing: center weight 1/2, four neighbors total 1/2.
F = 0.5*F + 0.125*(Fxp + Fxm + Fyp + Fym);
end

function [dFx, dFy] = gradient_cell_centered(F, params, periodicX, periodicY)
[Nx, Ny] = size(F);
dx = params.Lx / Nx;
dy = params.Ly / Ny;

ixp = [2:Nx, 1];
ixm = [Nx, 1:Nx-1];
if periodicX
    dFx = (F(ixp, :) - F(ixm, :)) / (2*dx);
else
    dFx = zeros(Nx, Ny);
    if Nx > 1
        dFx(1, :) = (F(2, :) - F(1, :)) / dx;
        dFx(Nx, :) = (F(Nx, :) - F(Nx-1, :)) / dx;
    end
    if Nx > 2
        dFx(2:Nx-1, :) = (F(3:Nx, :) - F(1:Nx-2, :)) / (2*dx);
    end
end

dFy = zeros(Nx, Ny);
if Ny > 1
    if periodicY
        iyp = [2:Ny, 1];
        iym = [Ny, 1:Ny-1];
        dFy = (F(:, iyp) - F(:, iym)) / (2*dy);
    else
        dFy(:, 1) = (F(:, 2) - F(:, 1)) / dy;
        dFy(:, Ny) = (F(:, Ny) - F(:, Ny-1)) / dy;
        if Ny > 2
            dFy(:, 2:Ny-1) = (F(:, 3:Ny) - F(:, 1:Ny-2)) / (2*dy);
        end
    end
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
