function [state, info] = resamp_initialize_particles_taylor_green_forced(params)
%RESAMP_INITIALIZE_PARTICLES_TAYLOR_GREEN_FORCED Initialize weighted periodic TG particles.
%
% This is the weighted/resampling branch equivalent of
% projection_initialize_particles_taylor_green_forced, but it adds
% state.m.  By default m(:)=1, so the weighted path is exactly compatible
% with the historical unit-mass particle interpretation.
%
% Output state fields:
%   state.x : Np-by-2 particle positions
%   state.v : Np-by-2 particle velocities
%   state.m : Np-by-1 particle masses/weights

params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'taylorGreenInitialAmplitude', get_param(params, 'taylorGreenAmplitude', 0.10));
params = set_default(params, 'taylorGreenModeX', 1);
params = set_default(params, 'taylorGreenModeY', 1);
params = set_default(params, 'taylorGreenThermalNoise', true);
params = set_default(params, 'kBT', 0.01);
params = set_default(params, 'resampParticleMass', get_param(params, 'particleMass', 1.0));

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
        error('Unknown initialPopulationMode for weighted TG: %s', params.initialPopulationMode);
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

m0 = params.resampParticleMass;
if ~(isnumeric(m0) && isscalar(m0) && isfinite(m0) && m0 > 0)
    error('params.resampParticleMass must be a positive finite scalar.');
end
m = m0 * ones(Np, 1);

state = struct('x', x, 'v', v, 'm', m);
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', true, 'periodicY', true, 'minMass', eps);
tg = projection_taylor_green_diagnostics(G, params);
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);

info = struct();
info.Np = Np;
info.mode = params.initialPopulationMode;
info.particleMass = m0;
info.nominalCellMass = gamma * m0;
info.nominalMassDensity = gamma * m0 / (dx * dy);
info.initialTGAmplitudeRequested = params.taylorGreenInitialAmplitude;
info.initialTGAmplitudeMeasured = tg.modeAmplitude;
info.initialTGCoherence = tg.modeCoherence;
info.initialTGModeEnergy = tg.modeEnergy;
info.initialPopulationStd = md.NStd;
info.initialMassStd = md.MStd;
info.initialPopulationOutBand = md.NOutBandFraction;
info.initialMassRelRms = md.MRelRms;
info.initialKBTWeighted = estimate_weighted_kBT(v, m);
info.massDiagnostics = md;
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

function kBT = estimate_weighted_kBT(v, m)
M = sum(m, 'omitnan');
if M <= 0
    kBT = NaN;
    return;
end
u = sum(m .* v, 1, 'omitnan') ./ M;
c = v - u;
kBT = 0.5 * sum(m .* sum(c.^2, 2), 'omitnan') ./ M;
end
