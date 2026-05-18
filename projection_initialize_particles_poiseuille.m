function [state, info] = projection_initialize_particles_poiseuille(params, varargin)
%PROJECTION_INITIALIZE_PARTICLES_POISEUILLE Initialize particles for channel runs.
%
%   [state, info] = projection_initialize_particles_poiseuille(params)
%
% Channel-specific initializer. It supports homogeneous thermal starts and,
% optionally, a deterministic Poiseuille velocity profile used to reduce the
% transient when calibrating viscosity:
%
%   params.initialPoiseuilleProfileEnable = true;
%   params.initialPoiseuilleNuGuess       = 0.03;
%
% The added profile is
%
%   u_x(y) = Ubottom + (Utop-Ubottom)*y/Ly + bodyForceX/(2*nuGuess)*y*(Ly-y),
%
% multiplied by params.initialPoiseuilleScale. It is added after the thermal
% random velocities have been zero-meaned, so the deterministic hydrodynamic
% profile is preserved.

if nargin < 1 || isempty(params)
    error('params is required.');
end

mode = get_param(params, 'initialPopulationMode', 'exact_per_cell');
zeroGlobalMean = logical(get_param(params, 'initialVelocityZeroGlobalMean', true));
initialMeanVelocityX = get_param(params, 'initialMeanVelocityX', 0.0);
initialMeanVelocityY = get_param(params, 'initialMeanVelocityY', 0.0);
poiseuilleProfileEnable = logical(get_param(params, 'initialPoiseuilleProfileEnable', false));
poiseuilleNuGuess = get_param(params, 'initialPoiseuilleNuGuess', 0.03);
poiseuilleScale = get_param(params, 'initialPoiseuilleScale', 1.0);
poiseuilleSubtractMean = logical(get_param(params, 'initialPoiseuilleSubtractMean', false));

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "initialpopulationmode"
            mode = char(val);
        case "initialvelocityzeroglobalmean"
            zeroGlobalMean = logical(val);
        case "initialmeanvelocityx"
            initialMeanVelocityX = val;
        case "initialmeanvelocityy"
            initialMeanVelocityY = val;
        case "initialpoiseuilleprofileenable"
            poiseuilleProfileEnable = logical(val);
        case "initialpoiseuillenuguess"
            poiseuilleNuGuess = val;
        case "initialpoiseuillescale"
            poiseuilleScale = val;
        case "initialpoiseuillesubtractmean"
            poiseuilleSubtractMean = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

mode = lower(strrep(char(mode), '-', '_'));
Nx = params.Nx;
Ny = params.Ny;
Lx = params.Lx;
Ly = params.Ly;
gamma = params.gamma;
kBT = params.kBT;

switch mode
    case {'random_uniform','random','uniform'}
        Np = round(gamma * Nx * Ny);
        x = [Lx * rand(Np, 1), Ly * rand(Np, 1)];
        canonicalMode = 'random_uniform';

    case {'exact_per_cell','uniform_by_cell','homogeneous_per_cell','cell_exact'}
        gammaInt = round(gamma);
        if abs(gamma - gammaInt) > 1e-12
            error('initialPopulationMode=exact_per_cell requires integer gamma. Received gamma=%.16g.', gamma);
        end
        if gammaInt < 0
            error('gamma must be non-negative.');
        end
        Np = gammaInt * Nx * Ny;
        x = zeros(Np, 2);
        dx = Lx / Nx;
        dy = Ly / Ny;
        cursor = 0;
        for ix = 0:Nx-1
            for iy = 0:Ny-1
                ids = cursor + (1:gammaInt);
                if gammaInt > 0
                    x(ids, 1) = (ix + rand(gammaInt, 1)) * dx;
                    x(ids, 2) = (iy + rand(gammaInt, 1)) * dy;
                end
                cursor = cursor + gammaInt;
            end
        end
        canonicalMode = 'exact_per_cell';

    otherwise
        error('Unknown initialPopulationMode: %s', mode);
end

v = sqrt(kBT) * randn(Np, 2);
if zeroGlobalMean && Np > 0
    v(:, 1) = v(:, 1) - mean(v(:, 1), 'omitnan');
    v(:, 2) = v(:, 2) - mean(v(:, 2), 'omitnan');
end
v(:, 1) = v(:, 1) + initialMeanVelocityX;
v(:, 2) = v(:, 2) + initialMeanVelocityY;

profileInfo = struct();
profileInfo.enabled = poiseuilleProfileEnable;
profileInfo.nuGuess = poiseuilleNuGuess;
profileInfo.scale = poiseuilleScale;
profileInfo.subtractMean = poiseuilleSubtractMean;
profileInfo.UmaxNoSlipFormula = NaN;
profileInfo.meanProfileAdded = 0.0;
profileInfo.centerMinusWallProfile = 0.0;

if poiseuilleProfileEnable && Np > 0
    if ~(isfinite(poiseuilleNuGuess) && poiseuilleNuGuess > 0)
        error('initialPoiseuilleNuGuess must be positive and finite. Received %.16g.', poiseuilleNuGuess);
    end
    y = x(:, 2);
    bodyForceX = get_param(params, 'bodyForceX', 0.0);
    Ubottom = get_param(params, 'Ubottom', 0.0);
    Utop = get_param(params, 'Utop', 0.0);
    couette = Ubottom + (Utop - Ubottom) .* (y ./ max(Ly, eps));
    poiseuille = bodyForceX ./ (2.0 * poiseuilleNuGuess) .* y .* (Ly - y);
    uxProfile = poiseuilleScale .* (couette + poiseuille);
    if poiseuilleSubtractMean
        uxProfile = uxProfile - mean(uxProfile, 'omitnan');
    end
    v(:, 1) = v(:, 1) + uxProfile;

    yMid = Ly / 2;
    centerProfile = poiseuilleScale .* (Ubottom + (Utop-Ubottom)*0.5 + bodyForceX/(2*poiseuilleNuGuess)*yMid*(Ly-yMid));
    wallProfileMean = poiseuilleScale .* mean([Ubottom, Utop], 'omitnan');
    if poiseuilleSubtractMean
        % Subtracting a constant does not alter center-wall amplitude.
    end
    profileInfo.UmaxNoSlipFormula = bodyForceX * Ly^2 / (8.0 * poiseuilleNuGuess);
    profileInfo.meanProfileAdded = mean(uxProfile, 'omitnan');
    profileInfo.centerMinusWallProfile = centerProfile - wallProfileMean;
end

state = struct();
state.x = x;
state.v = v;

pop = projection_population_diagnostics(state.x, params, 'periodicX', true, 'periodicY', false);
therm = projection_thermal_diagnostics(state.x, state.v, params, 'periodicX', true, 'periodicY', false, 'minCount', 1);

info = struct();
info.mode = canonicalMode;
info.Np = Np;
info.gamma = gamma;
info.zeroGlobalMean = zeroGlobalMean;
info.initialMeanVelocityX = initialMeanVelocityX;
info.initialMeanVelocityY = initialMeanVelocityY;
info.initialPoiseuilleProfile = profileInfo;
info.population = pop;
info.thermal = therm;
info.initialStdN = pop.stdN;
info.initialMeanN = pop.meanN;
info.initialMinN = pop.minN;
info.initialMaxN = pop.maxN;
info.initialEmptyCells = pop.nEmptyCells;
info.initialOutBand20 = pop.outBandFraction;
info.initialMaxAbsDeltaFromGamma = max(abs(double(pop.N(:)) - gamma));
info.initialRmsRel = sqrt(mean((double(pop.N(:))/gamma - 1).^2, 'omitnan'));
info.initialKBTCell = therm.kBTCellRelative;
end

function val = get_param(params, name, defaultVal)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultVal;
end
end
