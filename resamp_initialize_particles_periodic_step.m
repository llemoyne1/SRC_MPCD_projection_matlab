function [state, info] = resamp_initialize_particles_periodic_step(params)
%RESAMP_INITIALIZE_PARTICLES_PERIODIC_STEP Initialize weighted particles for a periodic step/shear test.
%
%   [state, info] = resamp_initialize_particles_periodic_step(params)
%
% This initializer is deliberately independent of the historical wall/backstep
% scripts.  It creates an exact-per-cell periodic particle support and assigns
% a discontinuous or smoothed staircase velocity field.  The purpose is to
% stress the weighted-resampling machinery in a fully periodic setting before
% reintroducing walls or immersed solids.
%
% Output state fields:
%   state.x : Np-by-2 particle positions
%   state.v : Np-by-2 particle velocities
%   state.m : Np-by-1 particle masses/weights
%
% Useful params fields:
%   periodicStepULower, periodicStepUUpper
%   periodicStepVLower, periodicStepVUpper
%   periodicStepXFraction
%   periodicStepYLowFraction, periodicStepYHighFraction
%   periodicStepTransitionWidth
%   periodicStepSubtractMeanVelocity

params = set_default(params, 'initialPopulationMode', 'exact_per_cell');
params = set_default(params, 'initialVelocityZeroGlobalMean', true);
params = set_default(params, 'periodicStepULower', 0.20);
params = set_default(params, 'periodicStepUUpper', -0.05);
params = set_default(params, 'periodicStepVLower', 0.00);
params = set_default(params, 'periodicStepVUpper', 0.00);
params = set_default(params, 'periodicStepXFraction', 0.35);
params = set_default(params, 'periodicStepYLowFraction', 0.35);
params = set_default(params, 'periodicStepYHighFraction', 0.65);
params = set_default(params, 'periodicStepTransitionWidth', 0.0);
params = set_default(params, 'periodicStepThermalNoise', true);
params = set_default(params, 'periodicStepSubtractMeanVelocity', true);
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
        error('Unknown initialPopulationMode for periodic step: %s', params.initialPopulationMode);
end

[Ux, Uy, side, yInterface] = resamp_periodic_step_velocity_at_points(x, params);
v = [Ux, Uy];
if logical(params.periodicStepThermalNoise) && params.kBT > 0
    c = sqrt(params.kBT) * randn(size(v));
    if logical(params.initialVelocityZeroGlobalMean)
        c(:,1) = c(:,1) - mean(c(:,1), 'omitnan');
        c(:,2) = c(:,2) - mean(c(:,2), 'omitnan');
    end
    v = v + c;
end
if logical(params.periodicStepSubtractMeanVelocity) || logical(params.initialVelocityZeroGlobalMean)
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
md = resamp_population_mass_diagnostics(state, params, 'periodicX', true, 'periodicY', true);

info = struct();
info.Np = Np;
info.mode = params.initialPopulationMode;
info.particleMass = m0;
info.nominalCellMass = gamma * m0;
info.initialPopulationStd = md.NStd;
info.initialMassStd = md.MStd;
info.initialMassRelRms = md.MRelRms;
info.initialKBTWeighted = md.kBTWeighted;
info.initialTotalMomentum = G.totalMomentum;
info.periodicStepULower = params.periodicStepULower;
info.periodicStepUUpper = params.periodicStepUUpper;
info.periodicStepXFraction = params.periodicStepXFraction;
info.periodicStepYLowFraction = params.periodicStepYLowFraction;
info.periodicStepYHighFraction = params.periodicStepYHighFraction;
info.periodicStepTransitionWidth = params.periodicStepTransitionWidth;
info.initialLowerFraction = mean(side < 0, 'omitnan');
info.yInterfaceAtParticles = yInterface;
info.massDiagnostics = md;
end

function [Ux, Uy, side, yInterface] = resamp_periodic_step_velocity_at_points(x, params)
Lx = params.Lx;
Ly = params.Ly;
xn = mod(x(:,1), Lx) ./ Lx;
y = mod(x(:,2), Ly);
xs = params.periodicStepXFraction;
yLow = params.periodicStepYLowFraction * Ly;
yHigh = params.periodicStepYHighFraction * Ly;
yInterface = yLow * ones(size(y));
yInterface(xn >= xs) = yHigh;
side = y - yInterface;
w = params.periodicStepTransitionWidth;
if w > 0
    % Smooth transition, useful if the discontinuous initialization is too stiff.
    a = 0.5 * (1.0 + tanh(side ./ max(w, eps)));
else
    a = double(side >= 0);
end
Ux = (1 - a) * params.periodicStepULower + a * params.periodicStepUUpper;
Uy = (1 - a) * params.periodicStepVLower + a * params.periodicStepVUpper;
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
