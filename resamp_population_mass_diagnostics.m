function diag = resamp_population_mass_diagnostics(state, params, varargin)
%RESAMP_POPULATION_MASS_DIAGNOSTICS Diagnostics for weighted particles.
%
% Separates particle count N, cell mass M, individual particle mass m_p,
% weighted global momentum and weighted kinetic/thermal measures.

validate_state(state);
periodicX = true;
periodicY = true;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', eps);
N = double(G.N(:));
M = double(G.M(:));
m = state.m(:);
v = state.v;

gamma = get_param(params, 'gamma', mean(N, 'omitnan'));
nominalParticleMass = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', median(m, 'omitnan')));
nominalCellMass = get_param(params, 'resampTargetCellMass', gamma * nominalParticleMass);

Mtot = sum(m, 'omitnan');
Ptot = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];
if Mtot > 0
    Uglobal = Ptot ./ Mtot;
    c = v - Uglobal;
    kBTWeighted = 0.5 * sum(m .* sum(c.^2, 2), 'omitnan') ./ Mtot;
    kineticMeanWeighted = 0.5 * sum(m .* sum(v.^2, 2), 'omitnan') ./ Mtot;
else
    Uglobal = [NaN NaN];
    kBTWeighted = NaN;
    kineticMeanWeighted = NaN;
end

massErr = M - nominalCellMass;
diag = struct();
diag.Np = size(state.x,1);
diag.Nx = params.Nx;
diag.Ny = params.Ny;
diag.NCells = params.Nx * params.Ny;
diag.totalMass = Mtot;
diag.totalMomentum = Ptot;
diag.globalVelocityWeighted = Uglobal;
diag.meanVxWeighted = Uglobal(1);
diag.meanVyWeighted = Uglobal(2);
diag.kBTWeighted = kBTWeighted;
diag.kineticEnergyMeanWeighted = kineticMeanWeighted;

diag.NMean = mean(N, 'omitnan');
diag.NStd = std(N, 0, 'omitnan');
diag.NMin = min(N);
diag.NMax = max(N);
diag.nEmptyCells = nnz(N == 0);
diag.NOutBandFraction = mean(abs(N - gamma) > 0.2 * max(gamma, eps), 'omitnan');

diag.MMean = mean(M, 'omitnan');
diag.MStd = std(M, 0, 'omitnan');
diag.MMin = min(M);
diag.MMax = max(M);
diag.MRelRms = sqrt(mean((massErr ./ max(nominalCellMass, eps)).^2, 'omitnan'));
diag.MOutBandFraction = mean(abs(M - nominalCellMass) > 0.2 * max(nominalCellMass, eps), 'omitnan');

diag.mParticleMean = mean(m, 'omitnan');
diag.mParticleStd = std(m, 0, 'omitnan');
diag.mParticleMin = min(m);
diag.mParticleMax = max(m);
diag.mParticleRelStd = diag.mParticleStd / max(diag.mParticleMean, eps);
diag.nominalParticleMass = nominalParticleMass;
diag.nominalCellMass = nominalCellMass;
diag.G = G;
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must be a struct with fields x, v and m.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
if numel(state.m) ~= size(state.x,1)
    error('state.m must have one entry per particle.');
end
if any(~isfinite(state.m(:))) || any(state.m(:) < 0)
    error('state.m must contain finite non-negative masses.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
