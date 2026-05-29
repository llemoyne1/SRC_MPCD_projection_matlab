function diag = resamp_population_mass_diagnostics(state, params, varargin)
%RESAMP_POPULATION_MASS_DIAGNOSTICS Diagnostics for weighted active particles.

validate_state(state);
periodicX = true;
periodicY = true;
cellWetMask = [];
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case {"cellwetmask", "wetmask", "fluidmask"}
            cellWetMask = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

activeMask = resamp_active_mask(state);
if isempty(cellWetMask)
    [cellWetMask, wetInfo] = resamp_cell_wet_mask(state, params, 'mode', 'auto');
else
    [cellWetMask, wetInfo] = resamp_cell_wet_mask(state, params, 'mode', 'explicit', 'cellWetMask', cellWetMask);
end
G = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', eps, 'activeMask', activeMask, ...
    'cellWetMask', cellWetMask);
Nall = double(G.N(:));
Mall = double(G.M(:));
wetVec = cellWetMask(:);
N = Nall(wetVec);
M = Mall(wetVec);
if isempty(N)
    N = zeros(0,1);
    M = zeros(0,1);
end
m = state.m(activeMask);
v = state.v(activeMask, :);

if isempty(m)
    m = zeros(0,1);
    v = zeros(0,2);
end

gamma = get_param(params, 'gamma', mean(N, 'omitnan'));
nominalParticleMass = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', 1.0));
if isempty(nominalParticleMass) || ~isfinite(nominalParticleMass)
    nominalParticleMass = 1.0;
end
nominalCellMass = get_param(params, 'resampTargetCellMass', gamma * nominalParticleMass);

Mtot = sum(m, 'omitnan');
if ~isempty(m)
    Ptot = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];
else
    Ptot = [0 0];
end
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
pool = resamp_particle_pool_info(state);

diag = struct();
diag.Np = pool.Nactive;
diag.NpActive = pool.Nactive;
diag.Ncapacity = pool.Ncapacity;
diag.Nfree = pool.Nfree;
diag.activeFraction = pool.activeFraction;
diag.Nx = params.Nx;
diag.Ny = params.Ny;
diag.NCells = params.Nx * params.Ny;
diag.nWetCells = wetInfo.nWetCells;
diag.nDryCells = wetInfo.nDryCells;
diag.wetFraction = wetInfo.wetFraction;
diag.totalMass = Mtot;
diag.totalMomentum = Ptot;
diag.globalVelocityWeighted = Uglobal;
diag.meanVxWeighted = Uglobal(1);
diag.meanVyWeighted = Uglobal(2);
diag.kBTWeighted = kBTWeighted;
diag.kineticEnergyMeanWeighted = kineticMeanWeighted;

diag.NMeanAllCells = mean(Nall, 'omitnan');
diag.NStdAllCells = std(Nall, 0, 'omitnan');
diag.NMinAllCells = min(Nall);
diag.NMaxAllCells = max(Nall);
diag.nEmptyCellsAllCells = nnz(Nall == 0);
if isempty(N)
    diag.NMean = NaN; diag.NStd = NaN; diag.NMin = NaN; diag.NMax = NaN;
    diag.nEmptyCells = NaN; diag.NOutBandFraction = NaN;
    diag.MMean = NaN; diag.MStd = NaN; diag.MMin = NaN; diag.MMax = NaN;
    diag.MRelRms = NaN; diag.MOutBandFraction = NaN;
else
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
end
diag.MMeanAllCells = mean(Mall, 'omitnan');
diag.MStdAllCells = std(Mall, 0, 'omitnan');
diag.MMinAllCells = min(Mall);
diag.MMaxAllCells = max(Mall);
diag.cellWetMask = cellWetMask;

if isempty(m)
    diag.mParticleMean = NaN;
    diag.mParticleStd = NaN;
    diag.mParticleMin = NaN;
    diag.mParticleMax = NaN;
    diag.mParticleRelStd = NaN;
else
    diag.mParticleMean = mean(m, 'omitnan');
    diag.mParticleStd = std(m, 0, 'omitnan');
    diag.mParticleMin = min(m);
    diag.mParticleMax = max(m);
    diag.mParticleRelStd = diag.mParticleStd / max(diag.mParticleMean, eps);
end

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
    error('state.m must have one entry per particle slot.');
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
