function summary = analyze_projection_poiseuille_viscosity(out, varargin)
%ANALYZE_PROJECTION_POISEUILLE_VISCOSITY Fit Poiseuille profile viscosity.
%
%   summary = analyze_projection_poiseuille_viscosity(out)
%
% Expects output from run_projection_poiseuille_demo or
% run_projection_poiseuille_medium64_demo. Uses a time-window average of
% stored profiles and fits:
%
%   Ux(y) = a2*y^2 + a1*y + a0,
%   nu_eff_raw = -bodyForceX/(2*a2).
%
% Scaling convention added for MPCD/SRC grid comparisons:
%   nu_eff_raw is the viscosity in the coordinate system used by y.
%   If the same physical domain Ly=1 is represented with different Ny, the
%   cell-size-normalized viscosity scales as nu_raw/dy^2.  For easy
%   comparison with a reference grid NyRef, this function also reports
%
%       nuEffScaledNyRef = nuEffRaw * (dyRef/dy)^2,
%       dyRef = physicalLy / NyRef.
%
% For Ly fixed, this is nuEffRaw * (Ny/NyRef)^2.  This is the quantity that
% should be compared across height-convergence runs when the MPCD transport
% coefficients are interpreted in cell units.
%
% Backward-compatible default:
%   fitStartFraction = params.fitStartFraction or 0.5.
%
% Additional transient-safe window options:
%   'fitStartTime', t0       : use t >= t0.
%   'fitWindowTime', DT      : use the trailing physical-time window
%                              t >= t(end)-DT.
%   'fitWindowFraction', f   : use the last fraction f of samples.
%   'fitLastNProfiles', n    : use the last n stored profiles.
%
% Scaling options:
%   'nuScalingEnable', tf          : keep scaled diagnostics enabled.
%   'nuReferenceNy', NyRef         : reference Ny for nuEffScaledNyRef.
%   'physicalLy', LyPhys          : physical/reference height, default params.Ly.
%   'nuScalingMode', mode         : currently 'Ny_over_reference'/'dy2'.
%
% Priority order for selecting fit profiles:
%   fitStartTime > fitLastNProfiles > fitWindowTime > fitWindowFraction >
%   tMin > fitStartFraction.

if nargin < 1 || ~isstruct(out)
    error('Usage: summary = analyze_projection_poiseuille_viscosity(out, ...)');
end

excludeWallCells = 2;
tMin = [];
fitStartFraction = [];
fitStartTime = [];
fitWindowTime = [];
fitWindowFraction = [];
fitLastNProfiles = [];
requirePositiveCurvature = false;
nuScalingEnableOverride = [];
nuReferenceNyOverride = [];
physicalLyOverride = [];
nuScalingModeOverride = [];

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
        case "fitstarttime"
            fitStartTime = val;
        case "fitwindowtime"
            fitWindowTime = val;
        case "fitwindowfraction"
            fitWindowFraction = val;
        case "fitlastnprofiles"
            fitLastNProfiles = val;
        case "requirepositivecurvature"
            requirePositiveCurvature = logical(val);
        case {"nuscalingenable", "poiseuillenu_scalingenable", "poiseuillenuscalingenable"}
            nuScalingEnableOverride = logical(val);
        case {"nureferenceny", "poiseuillenu_referenceny", "poiseuillenureferenceny", "nyref"}
            nuReferenceNyOverride = val;
        case {"physically", "poiseuillephysically", "lyphysical"}
            physicalLyOverride = val;
        case {"nuscalingmode", "poiseuillenu_scalingmode", "poiseuillenuscalingmode"}
            nuScalingModeOverride = char(string(val));
        otherwise
            error('Unknown option: %s', string(key));
    end
end

params = out.params;
y = out.yCenters(:);
UxProfiles = out.UxProfiles;
NProfiles = out.NProfiles;
t = out.sampleTimes(:);

if isempty(t) || isempty(UxProfiles)
    error('No stored Poiseuille profiles available for viscosity fit.');
end

idx = select_fit_indices(t, params, tMin, fitStartFraction, fitStartTime, ...
    fitWindowTime, fitWindowFraction, fitLastNProfiles);
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
nuEffRaw = -params.bodyForceX / (2*a2);
nuEff = nuEffRaw; % Backward-compatible alias: raw fit in y-coordinate units.

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
physicalCandidate = hasCorrectCurvature && isfinite(nuEffRaw) && nuEffRaw > 0 && isfinite(R2);

scale = poiseuille_nu_scaling(params, y, Ny, nuEffRaw, nuScalingEnableOverride, ...
    nuReferenceNyOverride, physicalLyOverride, nuScalingModeOverride);

summary = struct();
summary.tMin = min(t(idx));
summary.tMax = max(t(idx));
summary.nProfiles = nnz(idx);
summary.excludeWallCells = excludeWallCells;
summary.a2 = a2;
summary.a1 = coef(2);
summary.a0 = coef(3);
summary.nuEff = nuEff;
summary.nuEffRaw = nuEffRaw;
summary.nuEffCellY = scale.nuEffCellY;
summary.nuEffScaledNyRef = scale.nuEffScaledNyRef;
summary.nuEffDisplay = scale.nuEffDisplay;
summary.nuScaleFactorNyRef = scale.nuScaleFactorNyRef;
summary.nuReferenceNy = scale.nuReferenceNy;
summary.dy = scale.dy;
summary.dyRef = scale.dyRef;
summary.physicalLy = scale.physicalLy;
summary.nuScalingEnable = scale.nuScalingEnable;
summary.nuScalingMode = scale.nuScalingMode;
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
summary.fitIndexMask = idx;
summary.fitWindowDescription = describe_fit_window(t, idx);
end

function scale = poiseuille_nu_scaling(params, y, Ny, nuEffRaw, enableOverride, refNyOverride, physicalLyOverride, modeOverride)
scale = struct();

if numel(y) >= 2
    dy = median(diff(sort(y(:))), 'omitnan');
else
    dy = get_param(params, 'Ly', 1.0) / max(Ny, 1);
end
if ~isfinite(dy) || dy <= 0
    dy = get_param(params, 'Ly', 1.0) / max(Ny, 1);
end

physicalLy = get_param(params, 'poiseuillePhysicalLy', get_param(params, 'Ly', NaN));
if ~isempty(physicalLyOverride)
    physicalLy = physicalLyOverride;
end
if ~isfinite(physicalLy) || physicalLy <= 0
    physicalLy = max(y(:)) - min(y(:)) + dy;
end

nuReferenceNy = get_param(params, 'poiseuilleNuReferenceNy', 64);
if ~isempty(refNyOverride)
    nuReferenceNy = refNyOverride;
end
nuReferenceNy = max(1, round(nuReferenceNy));

dyRef = physicalLy / nuReferenceNy;
if ~isfinite(dyRef) || dyRef <= 0
    dyRef = dy;
end

nuScalingEnable = logical(get_param(params, 'poiseuilleNuScalingEnable', true));
if ~isempty(enableOverride)
    nuScalingEnable = logical(enableOverride);
end
nuScalingMode = char(string(get_param(params, 'poiseuilleNuScalingMode', 'Ny_over_reference')));
if ~isempty(modeOverride)
    nuScalingMode = char(string(modeOverride));
end

switch lower(strrep(nuScalingMode, '-', '_'))
    case {'ny_over_reference','dy2','dy_squared','cell_y'}
        factor = (dyRef / dy)^2;
    case {'none','raw'}
        factor = 1.0;
    otherwise
        warning('Unknown poiseuilleNuScalingMode=%s; using Ny_over_reference.', nuScalingMode);
        nuScalingMode = 'Ny_over_reference';
        factor = (dyRef / dy)^2;
end

nuEffCellY = nuEffRaw / max(dy^2, eps);
nuEffScaledNyRef = nuEffRaw * factor;
if nuScalingEnable
    nuEffDisplay = nuEffScaledNyRef;
else
    nuEffDisplay = nuEffRaw;
end

scale.nuEffRaw = nuEffRaw;
scale.nuEffCellY = nuEffCellY;
scale.nuEffScaledNyRef = nuEffScaledNyRef;
scale.nuEffDisplay = nuEffDisplay;
scale.nuScaleFactorNyRef = factor;
scale.nuReferenceNy = nuReferenceNy;
scale.dy = dy;
scale.dyRef = dyRef;
scale.physicalLy = physicalLy;
scale.nuScalingEnable = nuScalingEnable;
scale.nuScalingMode = nuScalingMode;
end

function val = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultValue;
end
end

function idx = select_fit_indices(t, params, tMin, fitStartFraction, fitStartTime, fitWindowTime, fitWindowFraction, fitLastNProfiles)
n = numel(t);
idx = false(size(t));

if ~isempty(fitStartTime) && isfinite(fitStartTime)
    idx = t >= fitStartTime;
    return;
end

if ~isempty(fitLastNProfiles) && isfinite(fitLastNProfiles) && fitLastNProfiles > 0
    nLast = max(1, min(n, round(fitLastNProfiles)));
    idx((n-nLast+1):n) = true;
    return;
end

if ~isempty(fitWindowTime) && isfinite(fitWindowTime) && fitWindowTime > 0
    t0 = max(min(t), max(t) - fitWindowTime);
    idx = t >= t0;
    return;
end

if ~isempty(fitWindowFraction) && isfinite(fitWindowFraction) && fitWindowFraction > 0
    f = min(max(fitWindowFraction, 1/n), 1.0);
    nLast = max(1, round(f * n));
    idx((n-nLast+1):n) = true;
    return;
end

if isempty(tMin)
    if isempty(fitStartFraction)
        if isfield(params, 'fitStartFraction') && ~isempty(params.fitStartFraction)
            fitStartFraction = params.fitStartFraction;
        else
            fitStartFraction = 0.5;
        end
    end
    fitStartFraction = min(max(fitStartFraction, 0), 0.95);
    idx0 = max(1, min(n, floor(fitStartFraction * n) + 1));
    tMin = t(idx0);
end
idx = t >= tMin;
end

function txt = describe_fit_window(t, idx)
if nnz(idx) < 1
    txt = 'empty';
else
    txt = sprintf('t=[%.6g, %.6g], n=%d/%d', min(t(idx)), max(t(idx)), nnz(idx), numel(t));
end
end
