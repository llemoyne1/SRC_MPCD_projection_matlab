function [state, diag] = resamp_regularize_particle_masses_preserve_mue(state, params, varargin)
%RESAMP_REGULARIZE_PARTICLE_MASSES_PRESERVE_MUE Regularize particle masses cell-wise.
%
%   [STATE, DIAG] = RESAMP_REGULARIZE_PARTICLE_MASSES_PRESERVE_MUE(STATE, PARAMS)
%   relaxes the masses of fluid particles in each wet cell toward a uniform
%   cell value while preserving, cell by cell:
%
%       M  = sum_p m_p
%       U  = sum_p m_p v_p / M
%       Eth = 1/2 sum_p m_p |v_p-U|^2
%
%   The velocity correction is a uniform shift followed by an optional
%   fluctuation rescale.  The function only acts on fluid-active particles
%   as returned by RESAMP_ACTIVE_MASK, hence latent/free particles are
%   ignored.  Dry cells are skipped.
%
%   The regularizer is intended as a long-run safety mechanism against a
%   slow broadening of the particle mass distribution.  It should not change
%   the hydrodynamic cell fields when preserveEnergy=true.

validate_state(state);

opts = parse_options(varargin{:});
if isempty(opts.nominalParticleMass)
    opts.nominalParticleMass = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', 1.0));
end
if isempty(opts.nominalParticleMass) || ~isfinite(opts.nominalParticleMass) || opts.nominalParticleMass <= 0
    opts.nominalParticleMass = 1.0;
end
opts.massMin = opts.minFactor * opts.nominalParticleMass;
opts.massMax = opts.maxFactor * opts.nominalParticleMass;

activeMask = opts.activeMask;
if isempty(activeMask)
    activeMask = resamp_active_mask(state);
else
    activeMask = logical(activeMask(:));
end

if isempty(opts.cellWetMask)
    [wetMask, wetInfo] = resamp_cell_wet_mask(state, params, 'mode', 'auto');
else
    [wetMask, wetInfo] = resamp_cell_wet_mask(state, params, 'mode', 'explicit', 'cellWetMask', opts.cellWetMask);
end
wetVec = wetMask(:);

cellId = resamp_cell_ids_periodic(state.x, params, ...
    'periodicX', opts.periodicX, 'periodicY', opts.periodicY, ...
    'activeMask', activeMask);

Nc = params.Nx * params.Ny;
cellRegularized = false(Nc,1);
cellSkippedDry = false(Nc,1);
cellSkippedFew = false(Nc,1);
cellTriggeredBounds = false(Nc,1);
cellTriggeredRelStd = false(Nc,1);
cellEnergyInfeasible = false(Nc,1);

relStdBefore = nan(Nc,1);
relStdAfter = nan(Nc,1);
mMinBefore = nan(Nc,1);
mMaxBefore = nan(Nc,1);
mMinAfter = nan(Nc,1);
mMaxAfter = nan(Nc,1);
energyScale = nan(Nc,1);
velocityShiftNorm = nan(Nc,1);
cellMassResidual = zeros(Nc,1);
cellMomentumResidual = zeros(Nc,2);
cellEnergyResidual = zeros(Nc,1);
cellEnergyResidualRel = nan(Nc,1);

for c = 1:Nc
    if ~wetVec(c)
        cellSkippedDry(c) = true;
        continue;
    end
    ids = find(activeMask & cellId == c);
    n = numel(ids);
    if n < opts.minParticlesPerCell
        cellSkippedFew(c) = true;
        continue;
    end

    mOld = state.m(ids);
    vOld = state.v(ids,:);
    M = sum(mOld, 'omitnan');
    if ~(isfinite(M) && M > opts.minCellMass)
        cellSkippedFew(c) = true;
        continue;
    end

    meanM = M / n;
    relStdBefore(c) = std(mOld,0,'omitnan') / max(mean(mOld,'omitnan'), eps);
    mMinBefore(c) = min(mOld);
    mMaxBefore(c) = max(mOld);
    boundsTrigger = (mMinBefore(c) < opts.massMin - opts.tolerance) || (mMaxBefore(c) > opts.massMax + opts.tolerance);
    relStdTrigger = relStdBefore(c) > opts.relStdTrigger;

    switch opts.triggerMode
        case {'always','all'}
            doRegularize = true;
        case {'bounds','bounds_only'}
            doRegularize = boundsTrigger;
        case {'relstd','relstd_only'}
            doRegularize = relStdTrigger;
        case {'bounds_or_relstd','conditional'}
            doRegularize = boundsTrigger || relStdTrigger;
        otherwise
            error('Unknown triggerMode: %s', opts.triggerMode);
    end

    if ~doRegularize
        continue;
    end

    cellTriggeredBounds(c) = boundsTrigger;
    cellTriggeredRelStd(c) = relStdTrigger;

    Uold = [sum(mOld .* vOld(:,1), 'omitnan'), sum(mOld .* vOld(:,2), 'omitnan')] ./ M;
    fluctOld = vOld - Uold;
    EthOld = 0.5 * sum(mOld .* sum(fluctOld.^2,2), 'omitnan');

    mTargetUniform = meanM * ones(n,1);
    strength = min(max(opts.strength,0),1);
    mRaw = (1-strength) .* mOld + strength .* mTargetUniform;
    mNew = local_project_mass_to_box(mRaw, M, opts.massMin, opts.massMax, opts.tolerance);

    % Restore the weighted mean velocity exactly by a uniform shift.
    UwithNewMass = [sum(mNew .* vOld(:,1), 'omitnan'), sum(mNew .* vOld(:,2), 'omitnan')] ./ M;
    shift = Uold - UwithNewMass;
    vShifted = vOld + shift;
    velocityShiftNorm(c) = norm(shift);

    if opts.preserveEnergy
        fluctShifted = vShifted - Uold;
        EthShifted = 0.5 * sum(mNew .* sum(fluctShifted.^2,2), 'omitnan');
        if EthOld <= opts.energyTolerance
            vNew = repmat(Uold, n, 1);
            energyScale(c) = 0;
        elseif EthShifted > opts.energyTolerance
            a = sqrt(EthOld / EthShifted);
            vNew = Uold + a .* fluctShifted;
            energyScale(c) = a;
        else
            vNew = vShifted;
            energyScale(c) = NaN;
            cellEnergyInfeasible(c) = true;
        end
    else
        vNew = vShifted;
        energyScale(c) = 1;
    end

    state.m(ids) = mNew;
    state.v(ids,:) = vNew;
    cellRegularized(c) = true;

    Mnew = sum(mNew, 'omitnan');
    Unew = [sum(mNew .* vNew(:,1), 'omitnan'), sum(mNew .* vNew(:,2), 'omitnan')] ./ max(Mnew, eps);
    fluctNew = vNew - Unew;
    EthNew = 0.5 * sum(mNew .* sum(fluctNew.^2,2), 'omitnan');

    relStdAfter(c) = std(mNew,0,'omitnan') / max(mean(mNew,'omitnan'), eps);
    mMinAfter(c) = min(mNew);
    mMaxAfter(c) = max(mNew);
    cellMassResidual(c) = Mnew - M;
    cellMomentumResidual(c,:) = Mnew .* Unew - M .* Uold;
    cellEnergyResidual(c) = EthNew - EthOld;
    cellEnergyResidualRel(c) = cellEnergyResidual(c) ./ max(abs(EthOld), opts.energyTolerance);
end

regularized = cellRegularized;
if any(regularized)
    massResidualRelRms = sqrt(mean((cellMassResidual(regularized) ./ max(opts.nominalParticleMass, eps)).^2, 'omitnan'));
    momentumResidualRms = sqrt(mean(sum(cellMomentumResidual(regularized,:).^2,2), 'omitnan'));
    energyResidualRms = sqrt(mean(cellEnergyResidual(regularized).^2, 'omitnan'));
    energyResidualRelRms = sqrt(mean(cellEnergyResidualRel(regularized).^2, 'omitnan'));
else
    massResidualRelRms = NaN;
    momentumResidualRms = NaN;
    energyResidualRms = NaN;
    energyResidualRelRms = NaN;
end

diag = struct();
diag.enabled = true;
diag.nCellsRegularized = nnz(cellRegularized);
regCells = find(cellRegularized);
if isempty(regCells)
    diag.nParticlesRegularized = 0;
else
    diag.nParticlesRegularized = nnz(activeMask & ismember(cellId, regCells));
end
diag.nCellsTriggeredBounds = nnz(cellTriggeredBounds);
diag.nCellsTriggeredRelStd = nnz(cellTriggeredRelStd);
diag.nCellsEnergyInfeasible = nnz(cellEnergyInfeasible);
diag.nCellsSkippedDry = nnz(cellSkippedDry);
diag.nCellsSkippedFew = nnz(cellSkippedFew);
diag.nWetCells = wetInfo.nWetCells;
diag.nDryCells = wetInfo.nDryCells;
diag.strength = opts.strength;
diag.relStdTrigger = opts.relStdTrigger;
diag.minFactor = opts.minFactor;
diag.maxFactor = opts.maxFactor;
diag.massResidualRelRms = massResidualRelRms;
diag.momentumResidualRms = momentumResidualRms;
diag.energyResidualRms = energyResidualRms;
diag.energyResidualRelRms = energyResidualRelRms;
diag.meanRelStdBefore = mean(relStdBefore(regularized), 'omitnan');
diag.meanRelStdAfter = mean(relStdAfter(regularized), 'omitnan');
diag.maxRelStdBefore = max(relStdBefore(regularized), [], 'omitnan');
diag.maxRelStdAfter = max(relStdAfter(regularized), [], 'omitnan');
diag.minMassBefore = min(mMinBefore(regularized), [], 'omitnan');
diag.maxMassBefore = max(mMaxBefore(regularized), [], 'omitnan');
diag.minMassAfter = min(mMinAfter(regularized), [], 'omitnan');
diag.maxMassAfter = max(mMaxAfter(regularized), [], 'omitnan');
diag.meanEnergyScale = mean(energyScale(regularized), 'omitnan');
diag.maxVelocityShiftNorm = max(velocityShiftNorm(regularized), [], 'omitnan');
diag.cellRegularized = reshape(cellRegularized, params.Nx, params.Ny);
diag.cellRelStdBefore = reshape(relStdBefore, params.Nx, params.Ny);
diag.cellRelStdAfter = reshape(relStdAfter, params.Nx, params.Ny);
end

function m = local_project_mass_to_box(m0, Mtarget, mMin, mMax, tol)
m = min(max(m0(:), mMin), mMax);
if numel(m) * mMin - Mtarget > 100*tol || Mtarget - numel(m) * mMax > 100*tol
    % Infeasible box for the requested total mass: keep the exact mass by
    % uniform scaling after clipping.  This should not occur in the nominal
    % resampling bands, but avoids a hard failure in diagnostics.
    m = m .* (Mtarget / max(sum(m), eps));
    return;
end
for it = 1:50
    r = Mtarget - sum(m, 'omitnan');
    if abs(r) <= tol * max(1, abs(Mtarget))
        break;
    end
    if r > 0
        free = m < mMax - tol;
    else
        free = m > mMin + tol;
    end
    if ~any(free)
        break;
    end
    m(free) = m(free) + r / nnz(free);
    m = min(max(m, mMin), mMax);
end
% Last tiny exact correction on unconstrained entries.
r = Mtarget - sum(m, 'omitnan');
if abs(r) > tol * max(1, abs(Mtarget))
    free = (m > mMin + tol) & (m < mMax - tol);
    if any(free)
        m(free) = m(free) + r / nnz(free);
    else
        m = m .* (Mtarget / max(sum(m), eps));
    end
end
end

function opts = parse_options(varargin)
opts = struct();
opts.periodicX = true;
opts.periodicY = true;
opts.activeMask = [];
opts.cellWetMask = [];
opts.strength = 1.0;
opts.triggerMode = 'bounds_or_relstd';
opts.relStdTrigger = 0.20;
opts.minFactor = 0.25;
opts.maxFactor = 4.0;
opts.nominalParticleMass = [];
opts.minParticlesPerCell = 2;
opts.minCellMass = eps;
opts.preserveEnergy = true;
opts.tolerance = 1e-12;
opts.energyTolerance = 1e-14;
if mod(numel(varargin),2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
        case 'periodicx', opts.periodicX = logical(val);
        case 'periodicy', opts.periodicY = logical(val);
        case 'activemask', opts.activeMask = logical(val);
        case {'cellwetmask','wetmask','fluidmask'}, opts.cellWetMask = logical(val);
        case 'strength', opts.strength = val;
        case 'triggermode', opts.triggerMode = lower(char(string(val)));
        case {'relstdtrigger','massrelstdtrigger'}, opts.relStdTrigger = val;
        case {'minfactor','massminfactor'}, opts.minFactor = val;
        case {'maxfactor','massmaxfactor'}, opts.maxFactor = val;
        case {'nominalparticlemass','particlemass'}, opts.nominalParticleMass = val;
        case 'minparticlespercell', opts.minParticlesPerCell = val;
        case 'mincellmass', opts.minCellMass = val;
        case {'preserveenergy','energycorrection'}, opts.preserveEnergy = logical(val);
        case {'tolerance','constrainttolerance'}, opts.tolerance = val;
        case 'energytolerance', opts.energyTolerance = val;
        otherwise, error('Unknown option: %s', key);
    end
end
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state,'x') || ~isfield(state,'v') || ~isfield(state,'m')
    error('state must contain x, v and m.');
end
if size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1) || numel(state.m) ~= size(state.x,1)
    error('state.x, state.v and state.m have inconsistent sizes.');
end
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
