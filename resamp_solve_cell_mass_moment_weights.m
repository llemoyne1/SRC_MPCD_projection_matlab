function [mNew, info] = resamp_solve_cell_mass_moment_weights(v, mOld, Mtarget, Ptarget, varargin)
%RESAMP_SOLVE_CELL_MASS_MOMENT_WEIGHTS Minimum-change local mass remap.
%
%   [mNew, info] = resamp_solve_cell_mass_moment_weights(v, mOld, Mtarget, Ptarget)
%
% Solves, for one cell, the bounded least-change problem
%
%       minimize    ||mNew - mOld||_2
%       subject to  sum(mNew)       = Mtarget
%                   sum(mNew .* vx) = Ptarget(1)
%                   sum(mNew .* vy) = Ptarget(2)
%                   massMin <= mNew <= massMax
%
% using a small active-set projection.  No Optimization Toolbox is needed.
% The routine is intended as the algebraic core for weighted SRC/MPCD
% resampling.  It does not create or delete particles.
%
% Options:
%   'massMin'              default 0
%   'massMax'              default Inf
%   'constraintTolerance'  default 1e-10
%   'maxIterations'        default numel(mOld)+8
%
% The output info.success means that the equality residual is below the
% requested tolerance and all bounds are satisfied.  When this is false,
% mNew is still the best effort returned by the active-set projection.

if nargin < 4
    error('Usage: [mNew, info] = resamp_solve_cell_mass_moment_weights(v, mOld, Mtarget, Ptarget, ...)');
end
if size(v, 2) ~= 2
    error('v must be NpCell-by-2.');
end
mOld = mOld(:);
Np = numel(mOld);
if size(v, 1) ~= Np
    error('mOld and v must contain the same number of particles.');
end
if Np == 0
    error('Cannot solve mass/moment weights for an empty cell.');
end
if ~isscalar(Mtarget) || ~isfinite(Mtarget) || Mtarget < 0
    error('Mtarget must be a finite non-negative scalar.');
end
Ptarget = Ptarget(:);
if numel(Ptarget) ~= 2 || any(~isfinite(Ptarget))
    error('Ptarget must be a finite 2-vector.');
end
if any(~isfinite(mOld)) || any(mOld < 0)
    error('mOld must contain finite non-negative masses.');
end
if any(~isfinite(v(:)))
    error('v must be finite.');
end

massMin = 0.0;
massMax = Inf;
tol = 1e-10;
maxIterations = Np + 8;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case {"massmin", "mmin"}
            massMin = val;
        case {"massmax", "mmax"}
            massMax = val;
        case {"constrainttolerance", "tol"}
            tol = val;
        case "maxiterations"
            maxIterations = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end
if ~isscalar(massMin)
    massMin = massMin(:);
    if numel(massMin) ~= Np
        error('massMin must be scalar or one value per particle.');
    end
else
    massMin = repmat(massMin, Np, 1);
end
if ~isscalar(massMax)
    massMax = massMax(:);
    if numel(massMax) ~= Np
        error('massMax must be scalar or one value per particle.');
    end
else
    massMax = repmat(massMax, Np, 1);
end
if any(~isfinite(massMin)) || any(massMin < 0)
    error('massMin must be finite and non-negative.');
end
if any(massMax < massMin)
    error('massMax must be >= massMin.');
end

A = [ones(1, Np); v(:,1).'; v(:,2).'];
b = [Mtarget; Ptarget(1); Ptarget(2)];

x0 = min(max(mOld, massMin), massMax);
x = x0;
free = true(Np, 1);
activeLower = false(Np, 1);
activeUpper = false(Np, 1);
iterations = 0;

for it = 1:maxIterations
    iterations = it;
    rhs = b - A(:, ~free) * x(~free);
    Af = A(:, free);
    x0f = x0(free);

    if isempty(x0f)
        break;
    end

    % Orthogonal projection of x0f onto Af*x = rhs.  pinv handles the
    % rank-deficient cases that appear in nearly monokinetic cells.
    corr = Af.' * (pinv(Af * Af.') * (rhs - Af * x0f));
    xf = x0f + corr;
    x(free) = xf;

    violLow = free & (x < massMin - tol);
    violHigh = free & (x > massMax + tol);
    if ~any(violLow | violHigh)
        break;
    end

    % Activate the most severe offending bounds first.  Activating all
    % offenders at once can over-constrain tiny cells unnecessarily.
    relLow = zeros(Np, 1);
    relHigh = zeros(Np, 1);
    relLow(violLow) = (massMin(violLow) - x(violLow)) ./ max(abs(massMin(violLow)) + 1, 1);
    relHigh(violHigh) = (x(violHigh) - massMax(violHigh)) ./ max(abs(massMax(violHigh)) + 1, 1);
    [~, idx] = max(max(relLow, relHigh));
    if violLow(idx)
        x(idx) = massMin(idx);
        activeLower(idx) = true;
    else
        x(idx) = massMax(idx);
        activeUpper(idx) = true;
    end
    free(idx) = false;
end

% Final clipping protects against roundoff after a successful projection.
x = min(max(x, massMin), massMax);
residual = A * x - b;
residualNorm = norm(residual);
scale = max([abs(Mtarget), norm(Ptarget), 1]);
residualRel = residualNorm / scale;
boundViolation = max([max(massMin - x), max(x - massMax), 0]);

mNew = x;
info = struct();
info.success = residualNorm <= tol * scale && boundViolation <= 10 * tol;
info.method = 'bounded_min_change_active_set';
info.NpCell = Np;
info.iterations = iterations;
info.rankA = rank(A);
info.nFree = nnz(free);
info.nActiveLower = nnz(activeLower);
info.nActiveUpper = nnz(activeUpper);
info.massResidual = residual(1);
info.momentumResidual = residual(2:3).';
info.residualNorm = residualNorm;
info.residualRel = residualRel;
info.boundViolation = boundViolation;
info.massOld = sum(mOld);
info.massNew = sum(mNew);
info.massTarget = Mtarget;
info.momentumOld = [sum(mOld .* v(:,1)), sum(mOld .* v(:,2))];
info.momentumNew = [sum(mNew .* v(:,1)), sum(mNew .* v(:,2))];
info.momentumTarget = Ptarget(:).';
info.deltaMassL2 = norm(mNew - mOld);
info.deltaMassRelRms = sqrt(mean(((mNew - mOld) ./ max(mean(mOld), eps)).^2, 'omitnan'));
end
