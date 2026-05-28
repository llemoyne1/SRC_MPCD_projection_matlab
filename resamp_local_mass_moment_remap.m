function [stateOut, diag] = resamp_local_mass_moment_remap(state, params, varargin)
%RESAMP_LOCAL_MASS_MOMENT_REMAP Local conservative mass/moment remap.
%
%   [stateOut, diag] = resamp_local_mass_moment_remap(state, params)
%
% This is the second weighted-resampling building block.  It changes only
% particle masses; positions and velocities are left unchanged.  It does
% not create, delete, split or merge particles.
%
% By default the target velocity in each cell is the current weighted cell
% velocity.  Therefore the local target momentum is
%
%       Ptarget_c = Mtarget_c * Uold_c,
%
% which means the remap imposes the target cell mass while preserving the
% hydrodynamic cell velocity.  This is deliberately conservative: no force
% or velocity kick is introduced by the remap itself.
%
% Options:
%   'targetCellMass'       scalar or Nx-by-Ny grid, default gamma*m0
%   'targetVelocityMode'   'preserve_cell_velocity' or 'grid'
%   'targetUx','targetUy'  Nx-by-Ny grids, required for targetVelocityMode='grid'
%   'method'               'scale_preserve_velocity' or 'min_change'
%   'massMin'              scalar lower mass bound, default 0
%   'massMax'              scalar upper mass bound, default Inf
%   'periodicX'            default true
%   'periodicY'            default true
%   'minMass'              default eps
%   'constraintTolerance'  default 1e-10
%   'computeDiagnostics'   default true

validate_state(state);

Nx = get_param(params, 'Nx', []);
Ny = get_param(params, 'Ny', []);
gamma = get_param(params, 'gamma', []);
nominalParticleMass = get_param(params, 'resampParticleMass', get_param(params, 'particleMass', median(state.m(:), 'omitnan')));
if isempty(Nx) || isempty(Ny) || isempty(gamma)
    error('params must contain Nx, Ny and gamma.');
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
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Nc = Nx * Ny;
if isscalar(targetCellMass)
    targetMassVec = repmat(targetCellMass, Nc, 1);
else
    if ~isequal(size(targetCellMass), [Nx, Ny])
        error('targetCellMass grid must have size Nx-by-Ny.');
    end
    targetMassVec = reshape(targetCellMass.', [Nc, 1]);
end
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

Gbefore = resamp_deposit_weighted_to_grid(state.x, state.v, state.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', minMass);
cellId = Gbefore.cellId;
Nvec = reshape(Gbefore.N.', [Nc, 1]);
Mvec = reshape(Gbefore.M.', [Nc, 1]);
UxVec = reshape(Gbefore.Ux.', [Nc, 1]);
UyVec = reshape(Gbefore.Uy.', [Nc, 1]);

mOldAll = state.m(:);
mNewAll = mOldAll;
cellSuccess = false(Nc, 1);
cellSkippedEmpty = false(Nc, 1);
cellUsedSolver = false(Nc, 1);
cellBounded = false(Nc, 1);
massResidual = zeros(Nc, 1);
momentumResidualX = zeros(Nc, 1);
momentumResidualY = zeros(Nc, 1);
solverResidualRel = NaN(Nc, 1);
deltaMassRmsRelCell = NaN(Nc, 1);

for c = 1:Nc
    ids = find(cellId == c);
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
Gafter = resamp_deposit_weighted_to_grid(stateOut.x, stateOut.v, stateOut.m, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minMass', minMass);

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

mdAfter = resamp_population_mass_diagnostics(stateOut, params, 'periodicX', periodicX, 'periodicY', periodicY);

diag = struct();
diag.kind = 'local_mass_moment_remap';
diag.method = method;
diag.targetVelocityMode = targetVelocityMode;
diag.targetCellMassMean = mean(targetMassVec, 'omitnan');
diag.massMin = massMin;
diag.massMax = massMax;
diag.nCells = Nc;
diag.nCellsNonEmpty = nnz(~cellSkippedEmpty);
diag.nCellsEmpty = nnz(cellSkippedEmpty);
diag.nCellsSolved = nnz(cellSuccess);
diag.nCellsUnresolved = nnz(~cellSuccess & ~cellSkippedEmpty);
diag.nCellsUsedSolver = nnz(cellUsedSolver);
diag.nCellsBounded = nnz(cellBounded);
diag.successFractionNonEmpty = diag.nCellsSolved / max(diag.nCellsNonEmpty, 1);
diag.massResidualRms = sqrt(mean(finalMassResidual.^2, 'omitnan'));
diag.massResidualMaxAbs = max(abs(finalMassResidual));
diag.massResidualRelRms = sqrt(mean((finalMassResidual ./ max(abs(targetMassVec), eps)).^2, 'omitnan'));
diag.momentumResidualRms = sqrt(mean(finalMomentumResidualX.^2 + finalMomentumResidualY.^2, 'omitnan'));
diag.momentumResidualMaxNorm = max(hypot(finalMomentumResidualX, finalMomentumResidualY));
diag.solverResidualRelMax = max(solverResidualRel, [], 'omitnan');
diag.deltaMassMean = mean(mNewAll - mOldAll, 'omitnan');
diag.deltaMassRms = sqrt(mean((mNewAll - mOldAll).^2, 'omitnan'));
diag.deltaMassRelRms = sqrt(mean(((mNewAll - mOldAll) ./ max(mean(mOldAll), eps)).^2, 'omitnan'));
diag.deltaMassRelCellRmsMean = mean(deltaMassRmsRelCell, 'omitnan');
diag.mParticleMinAfter = min(mNewAll);
diag.mParticleMaxAfter = max(mNewAll);
diag.mParticleStdAfter = std(mNewAll, 0, 'omitnan');
diag.totalMassBefore = sum(mOldAll, 'omitnan');
diag.totalMassAfter = sum(mNewAll, 'omitnan');
diag.totalMassTarget = sum(targetMassVec, 'omitnan');
diag.totalMomentumBefore = [sum(mOldAll .* state.v(:,1), 'omitnan'), sum(mOldAll .* state.v(:,2), 'omitnan')];
diag.totalMomentumAfter = [sum(mNewAll .* state.v(:,1), 'omitnan'), sum(mNewAll .* state.v(:,2), 'omitnan')];
diag.massDiagnosticsAfter = mdAfter;
diag.cellSuccess = [];
diag.cellSkippedEmpty = [];
diag.massResidualGrid = [];
diag.momentumResidualXGrid = [];
diag.momentumResidualYGrid = [];
diag.Gbefore = [];
diag.Gafter = [];
if computeDiagnostics
    diag.cellSuccess = reshape(cellSuccess, [Ny, Nx]).';
    diag.cellSkippedEmpty = reshape(cellSkippedEmpty, [Ny, Nx]).';
    diag.massResidualGrid = reshape(finalMassResidual, [Ny, Nx]).';
    diag.momentumResidualXGrid = reshape(finalMomentumResidualX, [Ny, Nx]).';
    diag.momentumResidualYGrid = reshape(finalMomentumResidualY, [Ny, Nx]).';
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
