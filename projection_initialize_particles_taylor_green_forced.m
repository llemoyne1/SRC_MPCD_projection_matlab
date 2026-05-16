function [state, info] = projection_initialize_particles_taylor_green_forced(params)
%PROJECTION_INITIALIZE_PARTICLES_TAYLOR_GREEN_FORCED Initialize periodic forced Taylor-Green particles.
%
% Particle positions are initialized either exact_per_cell or random, and
% velocities are set to a Taylor-Green mode plus optional thermal noise.

params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'taylorGreenInitialAmplitude', get_param(params, 'taylorGreenAmplitude', 0.10));
params = set_default(params, 'taylorGreenModeX', 1);
params = set_default(params, 'taylorGreenModeY', 1);
params = set_default(params, 'taylorGreenThermalNoise', true);
params = set_default(params, 'kBT', 0.01);

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
gamma = params.gamma;
dx = Lx / Nx;
dy = Ly / Ny;
mode = lower(strrep(char(string(params.initialPopulationMode)), '-', '_'));

switch mode
    case {'exact_per_cell','exact'}
        Np = Nx * Ny * gamma;
        x = zeros(Np, 2);
        ip = 0;
        for ix = 1:Nx
            for iy = 1:Ny
                for k = 1:gamma
                    ip = ip + 1;
                    x(ip,1) = (ix - 1 + rand()) * dx;
                    x(ip,2) = (iy - 1 + rand()) * dy;
                end
            end
        end
    case {'random','poisson'}
        Np = round(Nx * Ny * gamma);
        x = [Lx * rand(Np,1), Ly * rand(Np,1)];
    otherwise
        error('Unknown initialPopulationMode for TG: %s', params.initialPopulationMode);
end

[Ux, Uy] = projection_taylor_green_mode_at_points(x, params, params.taylorGreenInitialAmplitude);
v = [Ux, Uy];
if logical(params.taylorGreenThermalNoise) && params.kBT > 0
    c = sqrt(params.kBT) * randn(size(v));
    if logical(params.initialVelocityZeroGlobalMean)
        c(:,1) = c(:,1) - mean(c(:,1), 'omitnan');
        c(:,2) = c(:,2) - mean(c(:,2), 'omitnan');
    end
    v = v + c;
end
if logical(params.initialVelocityZeroGlobalMean)
    v(:,1) = v(:,1) - mean(v(:,1), 'omitnan');
    v(:,2) = v(:,2) - mean(v(:,2), 'omitnan');
end

state = struct('x', x, 'v', v);
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
tg = projection_taylor_green_diagnostics(G, params);

info = struct();
info.Np = Np;
info.mode = params.initialPopulationMode;
info.initialTGAmplitudeRequested = params.taylorGreenInitialAmplitude;
info.initialTGAmplitudeMeasured = tg.modeAmplitude;
info.initialTGCoherence = tg.modeCoherence;
info.initialTGModeEnergy = tg.modeEnergy;
info.initialPopulationStd = std(double(G.N(:)));
info.initialPopulationOutBand = mean(abs(double(G.N(:)) - gamma) > 0.2*gamma);
info.initialKBTCell = estimate_kBT(v);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function v = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end

function kBT = estimate_kBT(v)
u = mean(v, 1);
c = v - u;
kBT = 0.5 * mean(sum(c.^2, 2), 'omitnan');
end
