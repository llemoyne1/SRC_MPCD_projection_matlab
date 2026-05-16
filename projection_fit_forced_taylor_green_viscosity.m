function fit = projection_fit_forced_taylor_green_viscosity(t, A, params, varargin)
%PROJECTION_FIT_FORCED_TAYLOR_GREEN_VISCOSITY Fit nu_eff from forced TG amplitude.
%
% Forced Taylor--Green amplitude model, for one divergence-free Fourier mode:
%
%     dA/dt = F - lambda A,      lambda = nu_eff * (kx^2 + ky^2)
%
% with F = params.taylorGreenForceAmplitude.
%
% The function reports three complementary estimates:
%   1. pointwise/mean balance: nu = (F - dA/dt)/(k2 A)
%   2. derivative least-squares fit with known forcing F
%   3. exponential relaxation fit with known forcing F
% It also reports a free-forcing derivative fit dA/dt = F_fit - lambda A.
%
% Usage:
%   fit = projection_fit_forced_taylor_green_viscosity(S.t, S.tgAmplitude, params)
%   fit = projection_fit_forced_taylor_green_viscosity(..., 'fitFraction', 0.6)
%   fit = projection_fit_forced_taylor_green_viscosity(..., 'fitWindow', [t0 t1])

opt = struct();
opt.fitFraction = 0.60;       % use last 60% of samples by default
opt.fitWindow = [];
opt.minAmplitude = 1e-10;
opt.smoothSamples = 0;        % optional moving average before derivative fit
opt.verbose = false;
opt = parse_opts(opt, varargin{:});

fit = init_fit_struct();
if nargin < 3 || isempty(params)
    params = struct();
end

t = t(:); A = A(:);
validAll = isfinite(t) & isfinite(A) & A > opt.minAmplitude;
if nnz(validAll) < 4
    fit.status = 'insufficient_samples';
    return;
end

t = t(validAll); A = A(validAll);

F = get_param(params, 'taylorGreenForceAmplitude', NaN);
Lx = get_param(params, 'Lx', 1.0); Ly = get_param(params, 'Ly', 1.0);
mx = get_param(params, 'taylorGreenModeX', 1); my = get_param(params, 'taylorGreenModeY', 1);
kx = 2*pi*mx/Lx; ky = 2*pi*my/Ly; k2 = kx^2 + ky^2;
fit.forceKnown = F; fit.kx = kx; fit.ky = ky; fit.k2 = k2;
if ~isfinite(F)
    fit.status = 'missing_force_amplitude';
    return;
end

if opt.smoothSamples > 1
    AfitAll = movmean(A, opt.smoothSamples, 'Endpoints', 'shrink');
else
    AfitAll = A;
end

% Fit window.
if ~isempty(opt.fitWindow)
    fitMask = t >= opt.fitWindow(1) & t <= opt.fitWindow(2);
else
    tMin = min(t); tMax = max(t);
    tStart = tMin + (1 - opt.fitFraction) * (tMax - tMin);
    fitMask = t >= tStart;
end
fitMask = fitMask & isfinite(AfitAll) & AfitAll > opt.minAmplitude;
if nnz(fitMask) < 4
    fit.status = 'insufficient_fit_window';
    return;
end

tw = t(fitMask); Aw = AfitAll(fitMask);
fit.tStart = tw(1); fit.tEnd = tw(end); fit.nFit = numel(tw);
fit.meanAmplitude = mean(Aw, 'omitnan');
fit.finalAmplitude = AfitAll(end);

% Derivative estimate on all samples, then restricted to fit window.
dA = derivative_nonuniform(t, AfitAll);
dAw = dA(fitMask);
validD = isfinite(dAw) & isfinite(Aw) & Aw > opt.minAmplitude;
if nnz(validD) >= 3
    Av = Aw(validD); dAv = dAw(validD);

    % Known forcing: F - dA/dt = lambda A.
    lambdaKnown = sum(Av .* (F - dAv), 'omitnan') / max(sum(Av.^2, 'omitnan'), eps);
    fit.lambdaDerivativeKnownF = lambdaKnown;
    fit.nuDerivativeKnownF = lambdaKnown / k2;
    predD = F - lambdaKnown * Av;
    fit.r2DerivativeKnownF = compute_r2(dAv, predD);

    % Free forcing: dA/dt = F_fit - lambda A.
    X = [ones(numel(Av),1), -Av(:)];
    beta = X \ dAv(:);
    Ffit = beta(1); lambdaFree = beta(2);
    fit.forceDerivativeFit = Ffit;
    fit.lambdaDerivativeFreeF = lambdaFree;
    fit.nuDerivativeFreeF = lambdaFree / k2;
    predFree = X * beta;
    fit.r2DerivativeFreeF = compute_r2(dAv, predFree);

    % Pointwise balance over fit window.
    nuPoint = (F - dAv) ./ max(k2 * Av, eps);
    fit.nuPointwiseMean = mean(nuPoint(isfinite(nuPoint)), 'omitnan');
    fit.nuPointwiseMedian = median(nuPoint(isfinite(nuPoint)), 'omitnan');
    fit.nuPointwiseStd = std(nuPoint(isfinite(nuPoint)), 'omitnan');
end

% Plateau estimate: valid only near stationary amplitude, still useful.
fit.nuPlateau = F / max(k2 * fit.meanAmplitude, eps);
fit.lambdaPlateau = fit.nuPlateau * k2;

% Exponential fit with known F: A = F/lambda + C exp(-lambda*(t-t0)).
try
    [lambdaExp, Cexp, r2Exp, sseExp] = fit_exponential_known_force(tw, Aw, F);
    fit.lambdaExpKnownF = lambdaExp;
    fit.nuExpKnownF = lambdaExp / k2;
    fit.CExpKnownF = Cexp;
    fit.r2ExpKnownF = r2Exp;
    fit.sseExpKnownF = sseExp;
    fit.AInfExpKnownF = F / max(lambdaExp, eps);
catch ME
    fit.expFitMessage = ME.message;
end

% Preferred estimate: exponential if plausible, otherwise derivative known-F.
fit.nuPreferred = fit.nuExpKnownF;
fit.lambdaPreferred = fit.lambdaExpKnownF;
fit.preferredMethod = 'exp_knownF';
if ~isfinite(fit.nuPreferred) || fit.nuPreferred <= 0
    fit.nuPreferred = fit.nuDerivativeKnownF;
    fit.lambdaPreferred = fit.lambdaDerivativeKnownF;
    fit.preferredMethod = 'derivative_knownF';
end
if ~isfinite(fit.nuPreferred) || fit.nuPreferred <= 0
    fit.nuPreferred = fit.nuPlateau;
    fit.lambdaPreferred = fit.lambdaPlateau;
    fit.preferredMethod = 'plateau';
end

fit.status = 'ok';
if opt.verbose
    fprintf('TG viscosity fit: nu=%g [%s], window=[%g,%g], n=%d\n', ...
        fit.nuPreferred, fit.preferredMethod, fit.tStart, fit.tEnd, fit.nFit);
end
end

function fit = init_fit_struct()
fit = struct();
fields = {'forceKnown','kx','ky','k2','tStart','tEnd','nFit','meanAmplitude','finalAmplitude', ...
    'lambdaDerivativeKnownF','nuDerivativeKnownF','r2DerivativeKnownF', ...
    'forceDerivativeFit','lambdaDerivativeFreeF','nuDerivativeFreeF','r2DerivativeFreeF', ...
    'nuPointwiseMean','nuPointwiseMedian','nuPointwiseStd', ...
    'nuPlateau','lambdaPlateau','lambdaExpKnownF','nuExpKnownF','CExpKnownF', ...
    'r2ExpKnownF','sseExpKnownF','AInfExpKnownF','nuPreferred','lambdaPreferred'};
for i=1:numel(fields), fit.(fields{i}) = NaN; end
fit.status = 'not_run';
fit.preferredMethod = '';
fit.expFitMessage = '';
end

function d = derivative_nonuniform(t, y)
n = numel(t); d = nan(n,1);
if n < 2, return; end
if n == 2
    d(:) = (y(2)-y(1))/max(t(2)-t(1), eps);
    return;
end
d(1) = (y(2)-y(1))/max(t(2)-t(1), eps);
d(n) = (y(n)-y(n-1))/max(t(n)-t(n-1), eps);
for i=2:n-1
    d(i) = (y(i+1)-y(i-1))/max(t(i+1)-t(i-1), eps);
end
end

function [lambdaBest, Cbest, r2, sseBest] = fit_exponential_known_force(t, A, F)
t = t(:); A = A(:); tau = t - t(1);
Amean = mean(A, 'omitnan');
% Initial lambda guess from plateau. Bound broadly in inverse time units.
lambda0 = abs(F / max(abs(Amean), eps));
if ~isfinite(lambda0) || lambda0 <= 0, lambda0 = 1; end
lo = max(lambda0/1e4, 1e-10);
hi = max(lambda0*1e4, lo*10);
% Also make sure upper bound can represent observed decay/growth over window.
dtSpan = max(tau) - min(tau);
if dtSpan > 0
    hi = max(hi, 50/dtSpan);
end
obj = @(z) sse_for_loglambda(z, tau, A, F);
[zBest, sseBest] = fminbnd(obj, log(lo), log(hi));
lambdaBest = exp(zBest);
base = F / max(lambdaBest, eps);
e = exp(-lambdaBest * tau);
Cbest = sum(e .* (A - base), 'omitnan') / max(sum(e.^2, 'omitnan'), eps);
Apred = base + Cbest * e;
r2 = compute_r2(A, Apred);
end

function sse = sse_for_loglambda(z, tau, A, F)
lambda = exp(z);
base = F / max(lambda, eps);
e = exp(-lambda * tau);
C = sum(e .* (A - base), 'omitnan') / max(sum(e.^2, 'omitnan'), eps);
r = A - (base + C * e);
sse = sum(r.^2, 'omitnan');
if ~isfinite(sse), sse = realmax; end
end

function r2 = compute_r2(y, yhat)
y = y(:); yhat = yhat(:);
valid = isfinite(y) & isfinite(yhat);
if nnz(valid) < 2, r2 = NaN; return; end
y = y(valid); yhat = yhat(valid);
ssRes = sum((y-yhat).^2);
ssTot = sum((y-mean(y)).^2);
if ssTot <= eps, r2 = NaN; else, r2 = 1 - ssRes/ssTot; end
end

function opt = parse_opts(opt, varargin)
if mod(numel(varargin),2) ~= 0
    error('Options must be name/value pairs.');
end
for i=1:2:numel(varargin)
    name = varargin{i}; value = varargin{i+1};
    if ~isfield(opt, name), error('Unknown option: %s', name); end
    opt.(name) = value;
end
end

function v = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end
