function acc = poiseuille_acceleration_diagnostics(outOrTable, windowTime)
%POISEUILLE_ACCELERATION_DIAGNOSTICS Estimate global acceleration in channel runs.
%
%   acc = poiseuille_acceleration_diagnostics(out, windowTime)
%   acc = poiseuille_acceleration_diagnostics(diagTable, windowTime)
%
% Computes slope(meanUx) over the trailing physical-time window and returns
% ratioAccel = slope/bodyForceX when bodyForceX is available.  This is the
% primary diagnostic for detecting non-stationary Poiseuille runs with too
% little wall momentum sink.

if nargin < 2 || isempty(windowTime)
    windowTime = 5.0;
end

bodyForceX = NaN;
if isstruct(outOrTable) && isfield(outOrTable, 'diagTable')
    T = outOrTable.diagTable;
    if isfield(outOrTable, 'params') && isfield(outOrTable.params, 'bodyForceX')
        bodyForceX = outOrTable.params.bodyForceX;
    end
else
    T = outOrTable;
end

acc = struct('slopeMeanUx', NaN, 'ratioAccel', NaN, 'bodyForceX', bodyForceX, ...
    'effectiveShearForceX', NaN, 't0', NaN, 't1', NaN, 'nSamples', 0);

if isempty(T) || ~istable(T) || ~ismember('t', T.Properties.VariableNames) || ~ismember('meanUx', T.Properties.VariableNames)
    return;
end

t = T.t(:);
u = T.meanUx(:);
valid = isfinite(t) & isfinite(u);
if nnz(valid) < 2
    return;
end
t = t(valid);
u = u(valid);

t1 = max(t);
t0 = max(min(t), t1 - windowTime);
idx = t >= t0;
if nnz(idx) < 2
    idx = true(size(t));
end
p = polyfit(t(idx), u(idx), 1);
acc.slopeMeanUx = p(1);
acc.t0 = min(t(idx));
acc.t1 = max(t(idx));
acc.nSamples = nnz(idx);
if isfinite(bodyForceX) && abs(bodyForceX) > eps
    acc.ratioAccel = acc.slopeMeanUx / bodyForceX;
    acc.effectiveShearForceX = bodyForceX - acc.slopeMeanUx;
end
end
