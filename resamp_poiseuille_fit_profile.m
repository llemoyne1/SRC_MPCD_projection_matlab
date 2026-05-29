function fit = resamp_poiseuille_fit_profile(y, Ux, bodyForceX, varargin)
%RESAMP_POISEUILLE_FIT_PROFILE Fit weighted Poiseuille profiles robustly.
%
% The canonical forced-channel relation is written as
%
%   u(y) = u_slip + A y (H-y),       A = bodyForceX/(2 nu_eff)
%
% so that the slip-aware viscosity estimate is nu_eff = bodyForceX/(2 A).
% The function also reports a no-slip one-parameter fit and a free quadratic
% fit for diagnostics.  It is intentionally independent from any runner so it
% can be used for live visualisation and post-processing.

y = y(:); Ux = Ux(:);
if numel(y) ~= numel(Ux)
    error('y and Ux must have the same number of entries.');
end

opts = parse_options(varargin{:});
if isempty(opts.channelHeight)
    % For cell-centered grids y is typically (i-0.5)dy, so the channel
    % height is slightly larger than max(y).  If not supplied, use the most
    % conservative geometric estimate.
    dy = median(diff(sort(y(isfinite(y)))));
    if isempty(dy) || ~isfinite(dy) || dy <= 0
        dy = 1;
    end
    H = max(y) + 0.5 * dy;
else
    H = opts.channelHeight;
end

valid = isfinite(y) & isfinite(Ux);
fitMask = valid;
n = numel(y);
nex = max(0, round(opts.excludeWallCells));
if nex > 0 && n > 2*nex
    fitMask(1:nex) = false;
    fitMask(n-nex+1:n) = false;
end
if nnz(fitMask) < 4
    fitMask = valid;
end
if nnz(fitMask) < 3
    error('Too few valid points for Poiseuille fit.');
end

yf = y(fitMask);
uf = Ux(fitMask);
b = yf .* (H - yf);

% Slip-aware fit: u = uslip + A*y*(H-y).
Xslip = [ones(size(b)), b];
coefSlip = Xslip \ uf;
uslip = coefSlip(1);
A = coefSlip(2);
UfitSlip = uslip + A .* y .* (H - y);
statSlip = fit_stats(Ux, UfitSlip, fitMask, opts.profileStd);
nuSlip = force_to_nu(bodyForceX, A);

% No-slip fit: u = A0*y*(H-y).
if sum(b.^2) > 0
    A0 = (b' * uf) / (b' * b);
else
    A0 = NaN;
end
UfitNoSlip = A0 .* y .* (H - y);
statNoSlip = fit_stats(Ux, UfitNoSlip, fitMask, opts.profileStd);
nuNoSlip = force_to_nu(bodyForceX, A0);

% Free quadratic diagnostic: u = a2*y^2 + a1*y + a0.
coefPoly = polyfit(yf, uf, 2);
UfitPoly = polyval(coefPoly, y);
statPoly = fit_stats(Ux, UfitPoly, fitMask, opts.profileStd);
nuPoly = -bodyForceX / max(2 * coefPoly(1), realmin) ;
if coefPoly(1) == 0 || ~isfinite(coefPoly(1)) || bodyForceX == 0
    nuPoly = NaN;
else
    nuPoly = -bodyForceX / (2 * coefPoly(1));
end

centerVelocity = interp1(y, Ux, H/2, 'linear', 'extrap');
wallMeanVelocity = mean([Ux(1), Ux(end)], 'omitnan');
centerMinusWall = centerVelocity - wallMeanVelocity;
slipRatio = abs(uslip) / max(abs(centerVelocity), eps);
wallSlipRatioObserved = abs(wallMeanVelocity) / max(abs(centerMinusWall), eps);

switch lower(opts.fitModel)
    case {'slip','slip_parabolic','offset'}
        chosenName = 'slip'; chosenNu = nuSlip; chosenR2 = statSlip.R2; chosenUfit = UfitSlip; chosenA = A;
    case {'noslip','no_slip','zero_wall'}
        chosenName = 'noslip'; chosenNu = nuNoSlip; chosenR2 = statNoSlip.R2; chosenUfit = UfitNoSlip; chosenA = A0;
    case {'poly','quadratic'}
        chosenName = 'poly'; chosenNu = nuPoly; chosenR2 = statPoly.R2; chosenUfit = UfitPoly; chosenA = -coefPoly(1);
    otherwise
        error('Unknown fitModel: %s', opts.fitModel);
end

fit = struct();
fit.fitModel = chosenName;
fit.y = y;
fit.UxMean = Ux;
fit.fitMask = fitMask;
fit.channelHeight = H;
fit.bodyForceX = bodyForceX;
fit.centerVelocity = centerVelocity;
fit.wallMeanVelocity = wallMeanVelocity;
fit.centerMinusWall = centerMinusWall;
fit.slipVelocity = uslip;
fit.slipRatio = slipRatio;
fit.wallSlipRatioObserved = wallSlipRatioObserved;
fit.nuEff = chosenNu;
fit.R2 = chosenR2;
fit.UxFit = chosenUfit;
fit.A = chosenA;
fit.signalToNoise = statSlip.signalToNoise;
fit.profileNoiseRms = statSlip.profileNoiseRms;
fit.physicalCandidate = isfinite(chosenNu) && chosenNu > 0 && sign_nonzero(chosenA) == sign_nonzero(bodyForceX) && chosenR2 > opts.minR2ForCandidate;

fit.slipFit = struct('uSlip', uslip, 'A', A, 'nuEff', nuSlip, 'R2', statSlip.R2, ...
    'rmse', statSlip.rmse, 'UxFit', UfitSlip, 'physicalCandidate', isfinite(nuSlip) && nuSlip > 0 && sign_nonzero(A)==sign_nonzero(bodyForceX));
fit.noSlipFit = struct('A', A0, 'nuEff', nuNoSlip, 'R2', statNoSlip.R2, ...
    'rmse', statNoSlip.rmse, 'UxFit', UfitNoSlip, 'physicalCandidate', isfinite(nuNoSlip) && nuNoSlip > 0 && sign_nonzero(A0)==sign_nonzero(bodyForceX));
fit.polyFit = struct('coef', coefPoly, 'nuEff', nuPoly, 'R2', statPoly.R2, ...
    'rmse', statPoly.rmse, 'UxFit', UfitPoly);
end

function nu = force_to_nu(g, A)
if ~isfinite(g) || ~isfinite(A) || A == 0 || g == 0
    nu = NaN;
else
    nu = g / (2 * A);
end
end

function s = sign_nonzero(x)
if ~isfinite(x) || abs(x) < eps
    s = 0;
else
    s = sign(x);
end
end

function st = fit_stats(U, Ufit, mask, profileStd)
r = U(mask) - Ufit(mask);
ssRes = sum(r.^2, 'omitnan');
mu = mean(U(mask), 'omitnan');
ssTot = sum((U(mask)-mu).^2, 'omitnan');
R2 = 1 - ssRes / max(ssTot, eps);
rmse = sqrt(mean(r.^2, 'omitnan'));
if isempty(profileStd)
    noise = NaN;
else
    ps = profileStd(:);
    if numel(ps) == numel(U)
        noise = sqrt(mean(ps(mask).^2, 'omitnan'));
    else
        noise = NaN;
    end
end
st = struct('R2', R2, 'rmse', rmse, 'profileNoiseRms', noise, ...
    'signalToNoise', abs(max(U(mask))-min(U(mask))) / max(noise, eps));
end

function opts = parse_options(varargin)
opts = struct();
opts.excludeWallCells = 2;
opts.channelHeight = [];
opts.fitModel = 'slip';
opts.profileStd = [];
opts.minR2ForCandidate = 0.5;
if mod(numel(varargin),2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
        case 'excludewallcells'
            opts.excludeWallCells = val;
        case {'channelheight','h','ly'}
            opts.channelHeight = val;
        case 'fitmodel'
            opts.fitModel = lower(char(string(val)));
        case {'profilestd','uxstd'}
            opts.profileStd = val;
        case 'minr2forcandidate'
            opts.minR2ForCandidate = val;
        otherwise
            error('Unknown option: %s', key);
    end
end
end
