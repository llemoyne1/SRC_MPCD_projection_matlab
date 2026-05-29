function [stateOut, info] = resamp_apply_population_heterogeneity(state, params, varargin)
%RESAMP_APPLY_POPULATION_HETEROGENEITY Create a distributed population defect.
%
%   [stateOut, info] = resamp_apply_population_heterogeneity(state, params, ...)
%
% This routine creates an initial heterogeneous population field while
% conserving the total number of active particles.  It moves particles from
% cells whose target population is below their current population into cells
% whose target population is above their current population.  It is intended
% as a designed intermediate test for the closed-loop recycling method:
%
%   extraction from overpopulated cells -> pool -> insertion in poor cells
%   -> local mass/momentum remap.
%
% Main options:
%   'populationStd'      target standard deviation in particles/cell
%   'populationMin'      lower clamp for target cell population
%   'populationMax'      upper clamp for target cell population
%   'velocityMode'       'taylor_green_at_new_position' (default) or 'keep'
%   'thermalNoise'       add kBT thermal noise for moved particles (default true)
%   'zeroGlobalMean'     remove final weighted mean velocity (default true)
%
% The routine does not change particle masses or active count.

validate_state(state);
Nx = params.Nx;
Ny = params.Ny;
Nc = Nx * Ny;
gamma = get_param(params, 'gamma', 20);
Lx = params.Lx;
Ly = params.Ly;
dx = Lx / Nx;
dy = Ly / Ny;

populationStd = 4.0;
populationMin = max(0, floor(gamma - 3 * populationStd));
populationMax = ceil(gamma + 3 * populationStd);
velocityMode = 'taylor_green_at_new_position';
thermalNoise = true;
zeroGlobalMean = true;
computeDiagnostics = true;
periodicX = true;
periodicY = true;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {'populationstd','heterogeneitystd','stdparticles','std'}
            populationStd = val;
        case {'populationmin','heterogeneitymin','nmininitial','initialnmin'}
            populationMin = val;
        case {'populationmax','heterogeneitymax','nmaxinitial','initialnmax'}
            populationMax = val;
        case {'velocitymode','heterogeneityvelocitymode'}
            velocityMode = lower(char(string(val)));
        case {'thermalnoise','heterogeneitythermalnoise'}
            thermalNoise = logical(val);
        case {'zeroglobalmean','velocityzeroglobalmean'}
            zeroGlobalMean = logical(val);
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case 'computediagnostics'
            computeDiagnostics = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

if ~(isnumeric(populationStd) && isscalar(populationStd) && isfinite(populationStd) && populationStd >= 0)
    error('populationStd must be a finite non-negative scalar.');
end
populationMin = max(0, round(populationMin));
populationMax = max(populationMin, round(populationMax));
if gamma < populationMin || gamma > populationMax
    error('The nominal gamma must lie within [populationMin,populationMax].');
end

stateOut = state;
activeBefore = resamp_active_mask(stateOut);
Gbefore = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'activeMask', activeBefore);
Nbefore = reshape(Gbefore.N.', [Nc, 1]);
cellId = Gbefore.cellId;
activeIds = find(activeBefore);
Ntotal = nnz(activeBefore);

[targetN, targetRaw] = generate_target_population(Nc, Ntotal, gamma, populationStd, populationMin, populationMax);
if sum(targetN) ~= Ntotal
    error('Internal error: target population does not conserve active particle count.');
end

surplus = max(0, Nbefore - targetN);
deficit = max(0, targetN - Nbefore);
nMove = sum(surplus);
if nMove ~= sum(deficit)
    error('Internal error: surplus and deficit do not match.');
end

movedIds = zeros(nMove, 1);
pos = 0;
for c = find(surplus(:) > 0).'
    ids = activeIds(cellId(activeIds) == c);
    nTake = surplus(c);
    if nTake > numel(ids)
        error('Cannot take requested particles from cell %d.', c);
    end
    ids = ids(randperm(numel(ids), nTake));
    movedIds((pos+1):(pos+nTake)) = ids(:);
    pos = pos + nTake;
end

moveTargetCells = zeros(nMove, 1);
pos = 0;
for c = find(deficit(:) > 0).'
    nPut = deficit(c);
    moveTargetCells((pos+1):(pos+nPut)) = c;
    pos = pos + nPut;
end
if nMove > 1
    moveTargetCells = moveTargetCells(randperm(nMove));
end

for kk = 1:nMove
    c = moveTargetCells(kk);
    [ix, iy] = cell_index_from_id(c, Ny);
    stateOut.x(movedIds(kk),1) = (ix - 1 + rand()) * dx;
    stateOut.x(movedIds(kk),2) = (iy - 1 + rand()) * dy;
end

switch velocityMode
    case {'keep','preserve','preserve_velocity'}
        % Leave moved particle velocities unchanged.
    case {'taylor_green_at_new_position','tg','tg_at_new_position'}
        if nMove > 0
            amp = get_param(params, 'taylorGreenInitialAmplitude', get_param(params, 'taylorGreenAmplitude', 0.1));
            [ux, uy] = projection_taylor_green_mode_at_points(stateOut.x(movedIds,:), params, amp);
            vNew = [ux, uy];
            if thermalNoise && get_param(params, 'kBT', 0.0) > 0
                kBT = get_param(params, 'kBT', 0.0);
                fluct = sqrt(kBT) * randn(nMove, 2);
                fluct(:,1) = fluct(:,1) - mean(fluct(:,1), 'omitnan');
                fluct(:,2) = fluct(:,2) - mean(fluct(:,2), 'omitnan');
                vNew = vNew + fluct;
            end
            stateOut.v(movedIds,:) = vNew;
        end
    otherwise
        error('Unknown heterogeneity velocityMode: %s', velocityMode);
end

if zeroGlobalMean
    activeAfterMove = resamp_active_mask(stateOut);
    M = sum(stateOut.m(activeAfterMove), 'omitnan');
    if M > 0
        uMean = [sum(stateOut.m(activeAfterMove) .* stateOut.v(activeAfterMove,1), 'omitnan'), ...
                 sum(stateOut.m(activeAfterMove) .* stateOut.v(activeAfterMove,2), 'omitnan')] ./ M;
        stateOut.v(activeAfterMove,1) = stateOut.v(activeAfterMove,1) - uMean(1);
        stateOut.v(activeAfterMove,2) = stateOut.v(activeAfterMove,2) - uMean(2);
    end
end

activeAfter = resamp_active_mask(stateOut);
stateOut.Nactive = nnz(activeAfter);
if isfield(stateOut, 'Ncapacity') || isfield(stateOut, 'active')
    stateOut.Ncapacity = size(stateOut.x, 1);
end

Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'activeMask', activeAfter);
Nafter = reshape(Gafter.N.', [Nc, 1]);

info = struct();
info.mode = 'population_heterogeneous';
info.populationStdRequested = populationStd;
info.populationMinClamp = populationMin;
info.populationMaxClamp = populationMax;
info.velocityMode = velocityMode;
info.thermalNoise = thermalNoise;
info.zeroGlobalMean = zeroGlobalMean;
info.NactiveBefore = nnz(activeBefore);
info.NactiveAfter = nnz(activeAfter);
info.Ncapacity = size(stateOut.x, 1);
info.nMovedParticles = nMove;
info.nDonorCells = nnz(surplus > 0);
info.nReceiverCells = nnz(deficit > 0);
info.targetNMin = min(targetN);
info.targetNMax = max(targetN);
info.targetNStd = std(double(targetN), 0, 'omitnan');
info.NMinBefore = min(Nbefore);
info.NMaxBefore = max(Nbefore);
info.NStdBefore = std(double(Nbefore), 0, 'omitnan');
info.NMinAfter = min(Nafter);
info.NMaxAfter = max(Nafter);
info.NStdAfter = std(double(Nafter), 0, 'omitnan');
info.MRelRmsAfter = NaN;
if computeDiagnostics
    mdAfter = resamp_population_mass_diagnostics(stateOut, params, 'periodicX', periodicX, 'periodicY', periodicY);
    info.MRelRmsAfter = mdAfter.MRelRms;
    info.massDiagnosticsAfter = mdAfter;
    info.targetNGrid = reshape(targetN, [Ny, Nx]).';
    info.targetRawGrid = reshape(targetRaw, [Ny, Nx]).';
    info.surplusGrid = reshape(surplus, [Ny, Nx]).';
    info.deficitGrid = reshape(deficit, [Ny, Nx]).';
end
end

function [targetN, raw] = generate_target_population(Nc, Ntotal, gamma, sigma, nMin, nMax)
if sigma == 0
    targetN = gamma * ones(Nc, 1);
    raw = targetN;
    return;
end
raw = gamma + sigma * randn(Nc, 1);
targetN = round(raw);
targetN = min(max(targetN, nMin), nMax);
targetN = adjust_sum_with_bounds(targetN, Ntotal, nMin, nMax);
end

function values = adjust_sum_with_bounds(values, targetSum, lowerBound, upperBound)
values = round(values(:));
maxIterations = numel(values) + abs(sum(values) - targetSum) + 100;
it = 0;
while sum(values) ~= targetSum
    it = it + 1;
    if it > maxIterations
        error('Could not adjust heterogeneous target population to requested total.');
    end
    diff = sum(values) - targetSum;
    if diff > 0
        candidates = find(values > lowerBound);
        if isempty(candidates)
            error('Cannot decrease target population while respecting lower bound.');
        end
        candidates = candidates(randperm(numel(candidates)));
        for jj = 1:numel(candidates)
            c = candidates(jj);
            take = min(diff, values(c) - lowerBound);
            values(c) = values(c) - take;
            diff = diff - take;
            if diff == 0
                break;
            end
        end
    else
        need = -diff;
        candidates = find(values < upperBound);
        if isempty(candidates)
            error('Cannot increase target population while respecting upper bound.');
        end
        candidates = candidates(randperm(numel(candidates)));
        for jj = 1:numel(candidates)
            c = candidates(jj);
            add = min(need, upperBound - values(c));
            values(c) = values(c) + add;
            need = need - add;
            if need == 0
                break;
            end
        end
    end
end
end

function [ix, iy] = cell_index_from_id(c, Ny)
ix = floor((c - 1) / Ny) + 1;
iy = c - Ny * (ix - 1);
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state, 'x') || ~isfield(state, 'v') || ~isfield(state, 'm')
    error('state must contain x, v and m.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state.x and state.v must be Np-by-2 arrays with matching particle count.');
end
if numel(state.m) ~= size(state.x,1)
    error('state.m must have one entry per particle slot.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
