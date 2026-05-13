function [state, info] = projection_initialize_particles_taylor_green(params)
%PROJECTION_INITIALIZE_PARTICLES_TAYLOR_GREEN Initialize a periodic TG vortex.
%
%   [state, info] = projection_initialize_particles_taylor_green(params)
%
% Particle positions are initialized with the standard population modes used
% elsewhere.  Velocities are set to a Taylor-Green vortex plus optional
% Maxwellian thermal fluctuations.

params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'taylorGreenAmplitude', 0.20);
params = set_default(params, 'taylorGreenModeX', 1);
params = set_default(params, 'taylorGreenModeY', 1);
params = set_default(params, 'taylorGreenThermalNoise', true);
params = set_default(params, 'kBT', 0.05);

[state, baseInfo] = projection_initialize_particles(params, ...
    'initialPopulationMode', params.initialPopulationMode, ...
    'initialVelocityZeroGlobalMean', false);

[Ux, Uy] = taylor_green_velocity_at_particles(state.x, params);
state.v = [Ux, Uy];
if logical(params.taylorGreenThermalNoise) && params.kBT > 0
    c = sqrt(params.kBT) * randn(size(state.v));
    if logical(params.initialVelocityZeroGlobalMean)
        c(:,1) = c(:,1) - mean(c(:,1), 'omitnan');
        c(:,2) = c(:,2) - mean(c(:,2), 'omitnan');
    end
    state.v = state.v + c;
end
if logical(params.initialVelocityZeroGlobalMean)
    state.v(:,1) = state.v(:,1) - mean(state.v(:,1), 'omitnan');
    state.v(:,2) = state.v(:,2) - mean(state.v(:,2), 'omitnan');
end

G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);
tg = projection_taylor_green_diagnostics(G, params);

info = baseInfo;
info.mode = [baseInfo.mode '_taylor_green_velocity'];
info.taylorGreenAmplitude = params.taylorGreenAmplitude;
info.taylorGreenModeX = params.taylorGreenModeX;
info.taylorGreenModeY = params.taylorGreenModeY;
info.taylorGreenThermalNoise = logical(params.taylorGreenThermalNoise);
info.initialTG = tg;
info.initialTGAmplitude = tg.modeAmplitude;
info.initialTGCoherence = tg.modeCoherence;
info.initialTGModeEnergy = tg.modeEnergy;
end

function [Ux, Uy] = taylor_green_velocity_at_particles(x, params)
A = params.taylorGreenAmplitude;
mx = params.taylorGreenModeX;
my = params.taylorGreenModeY;
kx = 2*pi*mx / params.Lx;
ky = 2*pi*my / params.Ly;
% The Ly/Lx factor keeps div(u)=0 if Lx and Ly or mode numbers differ:
% kx*A_x = ky*A_y.
Ay = A * kx / max(ky, eps);
Ux = A  * sin(kx*x(:,1)) .* cos(ky*x(:,2));
Uy = -Ay * cos(kx*x(:,1)) .* sin(ky*x(:,2));
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end
