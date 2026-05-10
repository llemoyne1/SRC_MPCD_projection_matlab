function summary = analyze_projection_poiseuille_viscosity(out, varargin)
%ANALYZE_PROJECTION_POISEUILLE_VISCOSITY Fit Poiseuille profile viscosity.
%
%   summary = analyze_projection_poiseuille_viscosity(out)
%
% Expects output from run_projection_poiseuille_demo. Uses the last half of
% stored profiles by default, excludes a few wall-adjacent rows, and fits:
%
%   Ux(y) = a2*y^2 + a1*y + a0,
%   nu_eff = -bodyForceX/(2*a2).

if nargin < 1 || ~isstruct(out)
    error('Usage: summary = analyze_projection_poiseuille_viscosity(out, ...)');
end

excludeWallCells = 2;
tMin = [];
fitStartFraction = [];
requirePositiveCurvature = false;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "excludewallcells"
            excludeWallCells = val;
        case "tmin"
            tMin = val;
        case "fitstartfraction"
            fitStartFraction = val;
        case "requirepositivecurvature"
            requirePositiveCurvature = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

params = out.params;
y = out.yCenters(:);
UxProfiles = out.UxProfiles;
NProfiles = out.NProfiles;
t = out.sampleTimes(:);

if isempty(tMin)
    if isempty(fitStartFraction)
        if isfield(params, 'fitStartFraction') && ~isempty(params.fitStartFraction)
            fitStartFraction = params.fitStartFraction;
        else
            fitStartFraction = 0.5;
        end
    end
    fitStartFraction = min(max(fitStartFraction, 0), 0.95);
    idx0 = max(1, min(numel(t), floor(fitStartFraction * numel(t)) + 1));
    tMin = t(idx0);
end
idx = t >= tMin;
if nnz(idx) < 1
    idx = true(size(t));
end

UxMean = mean(UxProfiles(:, idx), 2, 'omitnan');
UxStd = std(UxProfiles(:, idx), 0, 2, 'omitnan');
NMean = mean(NProfiles(:, idx), 2, 'omitnan');

Ny = numel(y);
fitMask = true(Ny, 1);
fitMask(1:min(excludeWallCells, Ny)) = false;
fitMask(max(1, Ny-excludeWallCells+1):Ny) = false;
fitMask = fitMask & isfinite(UxMean) & isfinite(NMean) & NMean > 0;
if nnz(fitMask) < 5
    fitMask = isfinite(UxMean) & isfinite(NMean) & NMean > 0;
end
if nnz(fitMask) < 3
    error('Too few valid y rows for quadratic fit.');
end

coef = polyfit(y(fitMask), UxMean(fitMask), 2);
UxFit = polyval(coef, y);
a2 = coef(1);
nuEff = -params.bodyForceX / (2*a2);
expectedCurvatureSign = -sign(params.bodyForceX);
if expectedCurvatureSign == 0
    hasCorrectCurvature = true;
else
    hasCorrectCurvature = sign(a2) == expectedCurvatureSign;
end
if requirePositiveCurvature && ~hasCorrectCurvature
    warning('Poiseuille fit curvature has wrong sign: a2=%g, bodyForceX=%g.', a2, params.bodyForceX);
end
res = UxMean(fitMask) - UxFit(fitMask);
SSres = sum(res.^2);
SStot = sum((UxMean(fitMask) - mean(UxMean(fitMask))).^2);
R2 = 1 - SSres / max(SStot, eps);
profileNoiseRms = sqrt(mean(UxStd(fitMask).^2, 'omitnan'));
centerVelocityTmp = interp1(y, UxMean, params.Ly/2, 'linear', 'extrap');
wallMeanVelocityTmp = mean([UxMean(1), UxMean(end)], 'omitnan');
centerMinusWallTmp = centerVelocityTmp - wallMeanVelocityTmp;
signalToNoise = abs(centerMinusWallTmp) / max(profileNoiseRms, eps);
physicalCandidate = hasCorrectCurvature && isfinite(nuEff) && nuEff > 0 && isfinite(R2);

summary = struct();
summary.tMin = tMin;
summary.tMax = max(t(idx));
summary.nProfiles = nnz(idx);
summary.excludeWallCells = excludeWallCells;
summary.a2 = a2;
summary.a1 = coef(2);
summary.a0 = coef(3);
summary.nuEff = nuEff;
summary.hasCorrectCurvature = hasCorrectCurvature;
summary.expectedCurvatureSign = expectedCurvatureSign;
summary.R2 = R2;
summary.UmaxFit = max(UxFit);
summary.Umean = mean(UxMean, 'omitnan');
summary.centerVelocity = centerVelocityTmp;
summary.wallMeanVelocity = wallMeanVelocityTmp;
summary.centerMinusWall = centerMinusWallTmp;
summary.profileNoiseRms = profileNoiseRms;
summary.signalToNoise = signalToNoise;
summary.physicalCandidate = physicalCandidate;
summary.NMean = mean(NMean, 'omitnan');
summary.UxMean = UxMean;
summary.UxStd = UxStd;
summary.NMeanY = NMean;
summary.UxFit = UxFit;
summary.y = y;
summary.fitMask = fitMask;
summary.bodyForceX = params.bodyForceX;
end
