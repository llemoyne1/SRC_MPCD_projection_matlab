function [stateOut, diag] = resamp_local_mass_moment_remap(state, params, varargin)
%RESAMP_LOCAL_MASS_MOMENT_REMAP Local conservative mass/moment remap.
%
% Supports preallocated pool states.  Only active particles are remapped;
% inactive storage slots are ignored and left unchanged.

validate_state(state);

Nx = get_param(params, 'Nx', []);
Ny = get_param(params, 'Ny', []);
gamma = get_param(params, 'gamma', []);
nominalParticleMass = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', median(state.m(resamp_active_mask(state)), 'omitnan')));
if isempty(Nx) || isempty(Ny) || isempty(gamma)
    error('params must contain Nx, Ny and gamma.');
end
if isempty(nominalParticleMass) || ~isfinite(nominalParticleMass) || nominalParticleMass <= 0
    nominalParticleMass = 1.0;
end

targetCellMass = get_param(params, 'resampTargetCellMass', gamma * nominalParticleMass);
targetVelocityMode = 'preserve_cell_velocity';
targetUx = [];
targetUy = [];
method = 'scale_preserve_velocity';
massMin = 0.0;
massMax = Inf;
periodicX = true;
periodicY = true;
minMass = eps;
tol = 1e-10;
computeDiagnostics = true;
cellWetMask = [];
massSafetyEnable = false;
massSafetyMinFactor = 0.25;
massSafetyMaxFactor = 4.0;
massSafetyMode = 'uniform_mass_velocity_shift';

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "targetcellmass"
            targetCellMass = val;
        case "targetvelocitymode"
            targetVelocityMode = lower(char(string(val)));
        case "targetux"
            targetUx = val;
        case "targetuy"
            targetUy = val;
        case "method"
            method = lower(char(string(val)));
        case {"massmin", "mmin"}
            massMin = val;
        case {"massmax", "mmax"}
            massMax = val;
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case {"minmass", "minm"}
            minMass = val;
        case {"constrainttolerance", "tol"}
            tol = val;
        case "computediagnostics"
            computeDiagnostics = logical(val);
        case {"cellwetmask", "wetmask", "fluidmask"}
            cellWetMask = logical(val);
        case {"masssafetyenable", "remapmasssafetyenable"}
            massSafetyEnable = logical(val);
        case {"masssafetyminfactor", "remapmasssafetyminfactor"}
            massSafetyMinFactor = val;
        case {"masssafetymaxfactor", "remapmasssafetymaxfactor"}
            massSafetyMaxFactor = val;
        case {"masssafetymode", "remapmasssafetymode"}
            massSafetyMode = lower(char(string(val)));
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Nc = Nx * Ny;
if isempty(cellWetMask)
    [cellWetMask, wetInfo] = resamp_cell_wet_mask(state, params, 'mode', 'auto');
else
    [cellWetMask, wetInfo] = resamp_cell_wet_mask(state, params, 'mode', 'explicit', 'cellWetMask', cellWetMask);
end
wetVec = reshape(cellWetMask.', [Nc, 1]);

if isscalar(targetCellMass)
    targetMassVec = repmat(targetCellMass, Nc, 1);
else
    if ~isequal(size(targetCellMass), [Nx, Ny])
        error('targetCellMass grid must have size Nx-by-Ny.');
    end
    targetMassVec = reshape(targetCellMass.', [Nc, 1]);
end
targetMassVec(~wetVec) = 0;
if strcmp(targetVelocityMode, 'grid')
    if isempty(targetUx) || isempty(targetUy)
        error('targetUx and targetUy are required when targetVelocityMode=''grid''.');
    end
    if ~isequal(size(targetUx), [Nx, Ny]) || ~isequal(size(targetUy), [Nx, Ny])
        error('targetUx and targetUy must have size Nx-by-Ny.');
    end
    targetUxVec = reshape(targetUx.', [Nc, 1]);
    targetUyVec = reshape(targetUy.', [Nc, 1]);
else
    targetUxVec = [];
    targetUyVec = [];
end

activeMask = resamp_active_mask(state);
Gbefore = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', minMass, 'activeMask', activeMask, 'cellWetMask', cellWetMask);
cellId = Gbefore.cellId;
Nvec = reshape(Gbefore.N.', [Nc, 1]);
Mvec = reshape(Gbefore.M.', [Nc, 1]);
UxVec = reshape(Gbefore.Ux.', [Nc, 1]);
UyVec = reshape(Gbefore.Uy.', [Nc, 1]);

mOldAll = state.m(:);
mNewAll = mOldAll;
vNewAll = state.v;
cellSuccess = false(Nc, 1);
cellSkippedEmpty = false(Nc, 1);
cellSkippedDry = false(Nc, 1);
cellUsedSolver = false(Nc, 1);
cellBounded = false(Nc, 1);
massResidual = zeros(Nc, 1);
momentumResidualX = zeros(Nc, 1);
momentumResidualY = zeros(Nc, 1);
solverResidualRel = NaN(Nc, 1);
deltaMassRmsRelCell = NaN(Nc, 1);
massSafetyApplied = false(Nc, 1);
massSafetyTriggeredLow = false(Nc, 1);
massSafetyTriggeredHigh = false(Nc, 1);
massSafetyInfeasible = false(Nc, 1);
massSafetyShiftNorm = NaN(Nc, 1);
massSafetyCandidateMinFactor = NaN(Nc, 1);
massSafetyCandidateMaxFactor = NaN(Nc, 1);

for c = 1:Nc
    if ~wetVec(c)
        cellSkippedDry(c) = true;
        massResidual(c) = 0;
        momentumResidualX(c) = 0;
        momentumResidualY(c) = 0;
        continue;
    end
    ids = find(activeMask & cellId == c);
    if isempty(ids)
        cellSkippedEmpty(c) = true;
        massResidual(c) = -targetMassVec(c);
        continue;
    end

    Mold = Mvec(c);
    Mtarget = targetMassVec(c);
    if Mold <= minMass
        cellSkippedEmpty(c) = true;
        massResidual(c) = Mold - Mtarget;
        continue;
    end

    switch targetVelocityMode
        case {'preserve_cell_velocity', 'preserve_velocity', 'current'}
            Utarget = [UxVec(c), UyVec(c)];
        case 'grid'
            Utarget = [targetUxVec(c), targetUyVec(c)];
        otherwise
            error('Unknown targetVelocityMode: %s', targetVelocityMode);
    end
    Ptarget = Mtarget * Utarget;

    mOld = mNewAll(ids);
    vCell = state.v(ids, :);
    switch method
        case {'scale_preserve_velocity', 'scale'}
            if strcmp(targetVelocityMode, 'grid')
                [mCandidate, sinfo] = resamp_solve_cell_mass_moment_weights(vCell, mOld, Mtarget, Ptarget, ...
                    'massMin', massMin, 'massMax', massMax, 'constraintTolerance', tol);
                cellUsedSolver(c) = true;
            else
                scale = Mtarget / max(Mold, minMass);
                mCandidate = mOld * scale;
                if any(mCandidate < massMin - tol) || any(mCandidate > massMax + tol)
                    [mCandidate, sinfo] = resamp_solve_cell_mass_moment_weights(vCell, mOld, Mtarget, Ptarget, ...
                        'massMin', massMin, 'massMax', massMax, 'constraintTolerance', tol);
                    cellUsedSolver(c) = true;
                else
                    sinfo = local_constraint_info(vCell, mOld, mCandidate, Mtarget, Ptarget, tol);
                end
            end
        case {'min_change', 'bounded_min_change'}
            [mCandidate, sinfo] = resamp_solve_cell_mass_moment_weights(vCell, mOld, Mtarget, Ptarget, ...
                'massMin', massMin, 'massMax', massMax, 'constraintTolerance', tol);
            cellUsedSolver(c) = true;
        otherwise
            error('Unknown remap method: %s', method);
    end

    candidateMinFactor = min(mCandidate) / max(nominalParticleMass, eps);
    candidateMaxFactor = max(mCandidate) / max(nominalParticleMass, eps);
    massSafetyCandidateMinFactor(c) = candidateMinFactor;
    massSafetyCandidateMaxFactor(c) = candidateMaxFactor;
    doSafety = false;
    switch massSafetyMode
        case {'off','none','disabled'}
            doSafety = false;
        case {'uniform_mass_velocity_shift','uniform_velocity_shift','conditional_uniform_velocity_shift'}
            doSafety = massSafetyEnable && (candidateMinFactor < massSafetyMinFactor - tol || ...
                candidateMaxFactor > massSafetyMaxFactor + tol);
        case {'always_uniform_mass_velocity_shift','always_uniform_velocity_shift'}
            doSafety = massSafetyEnable;
        otherwise
            error('Unknown massSafetyMode: %s', massSafetyMode);
    end

    if doSafety
        [mCandidate, vCandidate, sinfoSafety] = local_uniform_mass_velocity_shift(vCell, Mtarget, Utarget, ...
            nominalParticleMass, massSafetyMinFactor, massSafetyMaxFactor, tol);
        vNewAll(ids,:) = vCandidate;
        sinfo = sinfoSafety;
        massSafetyApplied(c) = true;
        massSafetyTriggeredLow(c) = candidateMinFactor < massSafetyMinFactor - tol;
        massSafetyTriggeredHigh(c) = candidateMaxFactor > massSafetyMaxFactor + tol;
        massSafetyInfeasible(c) = get_sinfo_field(sinfoSafety, 'safetyInfeasible', false);
        massSafetyShiftNorm(c) = get_sinfo_field(sinfoSafety, 'velocityShiftNorm', NaN);
    end

    mNewAll(ids) = mCandidate;
    cellSuccess(c) = sinfo.success;
    cellBounded(c) = isfield(sinfo, 'nActiveLower') && (sinfo.nActiveLower + sinfo.nActiveUpper > 0);
    massResidual(c) = sinfo.massResidual;
    momentumResidualX(c) = sinfo.momentumResidual(1);
    momentumResidualY(c) = sinfo.momentumResidual(2);
    solverResidualRel(c) = sinfo.residualRel;
    deltaMassRmsRelCell(c) = sqrt(mean(((mCandidate - mOld) ./ max(mean(mOld), eps)).^2, 'omitnan'));
end

stateOut = state;
stateOut.m = mNewAll;
stateOut.v = vNewAll;
stateOut.Nactive = nnz(activeMask);
stateOut.Ncapacity = size(stateOut.x,1);
Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', minMass, 'activeMask', activeMask, 'cellWetMask', cellWetMask);

Mafter = reshape(Gafter.M.', [Nc, 1]);
PxAfter = reshape(Gafter.Px.', [Nc, 1]);
PyAfter = reshape(Gafter.Py.', [Nc, 1]);
if strcmp(targetVelocityMode, 'grid')
    Utx = targetUxVec;
    Uty = targetUyVec;
else
    Utx = UxVec;
    Uty = UyVec;
end
Pxtarget = targetMassVec .* Utx;
Pytarget = targetMassVec .* Uty;
finalMassResidual = Mafter - targetMassVec;
finalMomentumResidualX = PxAfter - Pxtarget;
finalMomentumResidualY = PyAfter - Pytarget;
finalMassResidual(~wetVec) = 0;
finalMomentumResidualX(~wetVec) = 0;
finalMomentumResidualY(~wetVec) = 0;

mdAfter = resamp_population_mass_diagnostics(stateOut, params, 'periodicX', periodicX, 'periodicY', periodicY, 'cellWetMask', cellWetMask);
activeAfter = resamp_active_mask(stateOut);

diag = struct();
diag.kind = 'local_mass_moment_remap';
diag.method = method;
diag.targetVelocityMode = targetVelocityMode;
diag.targetCellMassMean = mean(targetMassVec, 'omitnan');
diag.massMin = massMin;
diag.massMax = massMax;
diag.massSafetyEnable = massSafetyEnable;
diag.massSafetyMode = massSafetyMode;
diag.massSafetyMinFactor = massSafetyMinFactor;
diag.massSafetyMaxFactor = massSafetyMaxFactor;
diag.nCells = Nc;
diag.nWetCells = wetInfo.nWetCells;
diag.nDryCells = wetInfo.nDryCells;
diag.nCellsSkippedDry = nnz(cellSkippedDry);
diag.nCellsNonEmpty = nnz(wetVec & ~cellSkippedEmpty);
diag.nCellsEmpty = nnz(wetVec & cellSkippedEmpty);
diag.nCellsSolved = nnz(cellSuccess);
diag.nCellsUnresolved = nnz(wetVec & ~cellSuccess & ~cellSkippedEmpty);
diag.nCellsUsedSolver = nnz(cellUsedSolver);
diag.nCellsBounded = nnz(cellBounded);
diag.nCellsMassSafetyApplied = nnz(massSafetyApplied);
diag.nParticlesMassSafetyApplied = sum(Nvec(massSafetyApplied));
diag.nCellsMassSafetyTriggeredLow = nnz(massSafetyTriggeredLow);
diag.nCellsMassSafetyTriggeredHigh = nnz(massSafetyTriggeredHigh);
diag.nCellsMassSafetyInfeasible = nnz(massSafetyInfeasible);
diag.massSafetyVelocityShiftRms = sqrt(mean(massSafetyShiftNorm(massSafetyApplied).^2, 'omitnan'));
diag.massSafetyVelocityShiftMax = max(massSafetyShiftNorm(massSafetyApplied), [], 'omitnan');
diag.massSafetyCandidateMinFactor = min(massSafetyCandidateMinFactor(wetVec), [], 'omitnan');
diag.massSafetyCandidateMaxFactor = max(massSafetyCandidateMaxFactor(wetVec), [], 'omitnan');
diag.successFractionNonEmpty = diag.nCellsSolved / max(diag.nCellsNonEmpty, 1);
if any(wetVec)
    diag.massResidualRms = sqrt(mean(finalMassResidual(wetVec).^2, 'omitnan'));
    diag.massResidualMaxAbs = max(abs(finalMassResidual(wetVec)));
    diag.massResidualRelRms = sqrt(mean((finalMassResidual(wetVec) ./ max(abs(targetMassVec(wetVec)), eps)).^2, 'omitnan'));
    diag.momentumResidualRms = sqrt(mean(finalMomentumResidualX(wetVec).^2 + finalMomentumResidualY(wetVec).^2, 'omitnan'));
    diag.momentumResidualMaxNorm = max(hypot(finalMomentumResidualX(wetVec), finalMomentumResidualY(wetVec)));
else
    diag.massResidualRms = NaN;
    diag.massResidualMaxAbs = NaN;
    diag.massResidualRelRms = NaN;
    diag.momentumResidualRms = NaN;
    diag.momentumResidualMaxNorm = NaN;
end
diag.solverResidualRelMax = max(solverResidualRel, [], 'omitnan');
diag.deltaMassMean = mean(mNewAll(activeAfter) - mOldAll(activeAfter), 'omitnan');
diag.deltaMassRms = sqrt(mean((mNewAll(activeAfter) - mOldAll(activeAfter)).^2, 'omitnan'));
diag.deltaMassRelRms = sqrt(mean(((mNewAll(activeAfter) - mOldAll(activeAfter)) ./ max(mean(mOldAll(activeAfter)), eps)).^2, 'omitnan'));
diag.deltaMassRelCellRmsMean = mean(deltaMassRmsRelCell, 'omitnan');
diag.mParticleMinAfter = min(mNewAll(activeAfter));
diag.mParticleMaxAfter = max(mNewAll(activeAfter));
diag.mParticleStdAfter = std(mNewAll(activeAfter), 0, 'omitnan');
diag.totalMassBefore = sum(mOldAll(activeMask), 'omitnan');
diag.totalMassAfter = sum(mNewAll(activeAfter), 'omitnan');
diag.totalMassTarget = sum(targetMassVec(wetVec), 'omitnan');
diag.totalMomentumBefore = [sum(mOldAll(activeMask) .* state.v(activeMask,1), 'omitnan'), sum(mOldAll(activeMask) .* state.v(activeMask,2), 'omitnan')];
diag.totalMomentumAfter = [sum(mNewAll(activeAfter) .* stateOut.v(activeAfter,1), 'omitnan'), sum(mNewAll(activeAfter) .* stateOut.v(activeAfter,2), 'omitnan')];
diag.NpActive = nnz(activeAfter);
diag.Ncapacity = size(stateOut.x,1);
diag.Nfree = diag.Ncapacity - diag.NpActive;
diag.massDiagnosticsAfter = mdAfter;
diag.cellSuccess = [];
diag.cellSkippedEmpty = [];
diag.massResidualGrid = [];
diag.momentumResidualXGrid = [];
diag.momentumResidualYGrid = [];
diag.massSafetyAppliedGrid = [];
diag.Gbefore = [];
diag.Gafter = [];
if computeDiagnostics
    diag.cellSuccess = reshape(cellSuccess, [Ny, Nx]).';
    diag.cellSkippedEmpty = reshape(cellSkippedEmpty, [Ny, Nx]).';
    diag.cellSkippedDry = reshape(cellSkippedDry, [Ny, Nx]).';
    diag.cellWetMask = cellWetMask;
    diag.massResidualGrid = reshape(finalMassResidual, [Ny, Nx]).';
    diag.momentumResidualXGrid = reshape(finalMomentumResidualX, [Ny, Nx]).';
    diag.momentumResidualYGrid = reshape(finalMomentumResidualY, [Ny, Nx]).';
    diag.massSafetyAppliedGrid = reshape(massSafetyApplied, [Ny, Nx]).';
    diag.Gbefore = Gbefore;
    diag.Gafter = Gafter;
end
end

function info = local_constraint_info(v, mOld, mNew, Mtarget, Ptarget, tol)
res = [sum(mNew) - Mtarget; sum(mNew .* v(:,1)) - Ptarget(1); sum(mNew .* v(:,2)) - Ptarget(2)];
scale = max([abs(Mtarget), norm(Ptarget), 1]);
info = struct();
info.success = norm(res) <= tol * scale;
info.massResidual = res(1);
info.momentumResidual = res(2:3).';
info.residualNorm = norm(res);
info.residualRel = norm(res) / scale;
info.nActiveLower = 0;
info.nActiveUpper = 0;
info.massOld = sum(mOld);
info.massNew = sum(mNew);
info.massTarget = Mtarget;
info.momentumOld = [sum(mOld .* v(:,1)), sum(mOld .* v(:,2))];
info.momentumNew = [sum(mNew .* v(:,1)), sum(mNew .* v(:,2))];
info.momentumTarget = Ptarget(:).';
info.deltaMassL2 = norm(mNew - mOld);
info.deltaMassRelRms = sqrt(mean(((mNew - mOld) ./ max(mean(mOld), eps)).^2, 'omitnan'));
end


function [mNew, vNew, info] = local_uniform_mass_velocity_shift(v, Mtarget, Utarget, nominalParticleMass, safetyMinFactor, safetyMaxFactor, tol)
%LOCAL_UNIFORM_MASS_VELOCITY_SHIFT Bounded fallback for edited populations.
%
% The fallback keeps the particle support well conditioned by assigning the
% same mass Mtarget/N to all particles in the cell, then applies one uniform
% velocity shift so that the cell momentum is exactly Mtarget*Utarget.  This
% preserves relative thermal fluctuations and avoids using extreme masses to
% satisfy the momentum constraint after extraction/insertion.

n = size(v, 1);
if n <= 0
    mNew = zeros(0,1);
    vNew = v;
    info = local_constraint_info(v, zeros(0,1), zeros(0,1), Mtarget, Mtarget * Utarget, tol);
    info.safetyInfeasible = true;
    info.velocityShift = [NaN, NaN];
    info.velocityShiftNorm = NaN;
    return;
end
mUniform = Mtarget / n;
mNew = repmat(mUniform, n, 1);
Ptarget = Mtarget * Utarget;
Ucurrent = sum(mNew .* v, 1) / max(sum(mNew), eps);
deltaU = Utarget - Ucurrent;
vNew = v + deltaU;
info = local_constraint_info(vNew, repmat(nominalParticleMass, n, 1), mNew, Mtarget, Ptarget, tol);
info.usedMassSafety = true;
info.velocityShift = deltaU;
info.velocityShiftNorm = norm(deltaU);
info.uniformMass = mUniform;
info.uniformMassFactor = mUniform / max(nominalParticleMass, eps);
info.safetyInfeasible = (info.uniformMassFactor < safetyMinFactor - tol) || ...
    (info.uniformMassFactor > safetyMaxFactor + tol);
info.nActiveLower = double(info.uniformMassFactor < safetyMinFactor - tol) * n;
info.nActiveUpper = double(info.uniformMassFactor > safetyMaxFactor + tol) * n;
end

function v = get_sinfo_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
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
