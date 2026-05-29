function summary = analyze_resamp_poiseuille_viscosity(out, varargin)
%ANALYZE_RESAMP_POISEUILLE_VISCOSITY Robust profile-based Poiseuille diagnostics.
%
% The returned struct keeps the historical fields nuEff/R2/UxFit, but now
% also reports slip-aware and no-slip fits separately.  The slip-aware fit is
% the default because wallVP or remaining wall slip should not make the
% curvature/viscosity estimate meaningless.

excludeWallCells = 2;
fitStartFraction = 0.5;
fitLastN = [];
fitModel = 'slip';
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
        case 'excludewallcells'
            excludeWallCells = val;
        case 'fitstartfraction'
            fitStartFraction = val;
        case 'fitlastnprofiles'
            fitLastN = val;
        case 'fitmodel'
            fitModel = lower(char(string(val)));
        otherwise
            error('Unknown option: %s', key);
    end
end

params = out.params;
y = out.yCenters(:);
U = out.UxProfiles;
N = out.NProfiles;
t = out.sampleTimes(:);
if isempty(t) || isempty(U)
    error('No profiles in out.');
end
if ~isempty(fitLastN)
    idx = false(size(t));
    idx(max(1, numel(t)-fitLastN+1):end) = true;
else
    t0 = min(t) + (1 - fitStartFraction) * (max(t) - min(t));
    idx = t >= t0;
end
if nnz(idx) < 1
    idx(:) = true;
end
UxMean = mean(U(:,idx), 2, 'omitnan');
UxStd = std(U(:,idx), 0, 2, 'omitnan');
NMean = mean(N(:,idx), 2, 'omitnan');

fit = resamp_poiseuille_fit_profile(y, UxMean, params.bodyForceX, ...
    'excludeWallCells', excludeWallCells, ...
    'channelHeight', params.Ly, ...
    'fitModel', fitModel, ...
    'profileStd', UxStd);

summary = fit;
summary.tMin = min(t(idx));
summary.tMax = max(t(idx));
summary.nProfiles = nnz(idx);
summary.excludeWallCells = excludeWallCells;
summary.NMeanY = NMean;
summary.UxStd = UxStd;
summary.hasCorrectCurvature = sign_nonzero(summary.A) == sign_nonzero(params.bodyForceX) || params.bodyForceX == 0;
summary.physicalCandidate = summary.physicalCandidate && summary.hasCorrectCurvature && isfinite(summary.nuEff) && summary.nuEff > 0;

% Historical aliases used by older summaries.
summary.a2 = -summary.A;
summary.a1 = NaN;
summary.a0 = summary.slipVelocity;
summary.nuEffSlip = summary.slipFit.nuEff;
summary.R2Slip = summary.slipFit.R2;
summary.nuEffNoSlip = summary.noSlipFit.nuEff;
summary.R2NoSlip = summary.noSlipFit.R2;
summary.UxFitSlip = summary.slipFit.UxFit;
summary.UxFitNoSlip = summary.noSlipFit.UxFit;
summary.UxFitPoly = summary.polyFit.UxFit;
end

function s = sign_nonzero(x)
if ~isfinite(x) || abs(x) < eps
    s = 0;
else
    s = sign(x);
end
end
