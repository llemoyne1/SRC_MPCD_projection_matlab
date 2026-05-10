function metrics = projection_density_homogeneity_metrics(NMaps, params, varargin)
%PROJECTION_DENSITY_HOMOGENEITY_METRICS Cell-density homogeneity diagnostics.
%
%   metrics = projection_density_homogeneity_metrics(NMaps, params)
%
% NMaps must be an Nx-by-Ny-by-Nt array of cell populations. The function
% quantifies instantaneous and time-averaged density homogeneity, including
% cell-wise relative fluctuations, out-of-band fractions, wall-excluded
% variants and low-wavenumber density energy.
%
% Name-value options:
%   'sampleTimes'       default 0:Nt-1
%   'targetGamma'       default params.gamma or mean(NMaps(:))
%   'bandFraction'      default 0.20
%   'excludeWallCells'  default 0
%   'lowKMaxIndex'      default 2
%
% The low-k metric is based on FFT2 of rhoRel = N/targetGamma - 1, excluding
% the zero mode. It is a compact measure of coherent large-scale density
% inhomogeneity.

if nargin < 2
    error('Usage: metrics = projection_density_homogeneity_metrics(NMaps, params, ...)');
end
if isempty(NMaps) || ndims(NMaps) ~= 3
    error('NMaps must be an Nx-by-Ny-by-Nt array. Make sure params.storeDensityMaps=true.');
end

[Nx, Ny, Nt] = size(NMaps);
sampleTimes = (0:Nt-1).';
targetGamma = [];
bandFraction = 0.20;
excludeWallCells = 0;
lowKMaxIndex = 2;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "sampletimes"
            sampleTimes = val(:);
        case "targetgamma"
            targetGamma = val;
        case "bandfraction"
            bandFraction = val;
        case "excludewallcells"
            excludeWallCells = val;
        case "lowkmaxindex"
            lowKMaxIndex = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end
if isempty(targetGamma)
    if isfield(params, 'gamma') && ~isempty(params.gamma)
        targetGamma = params.gamma;
    else
        targetGamma = mean(NMaps(:), 'omitnan');
    end
end
if numel(sampleTimes) ~= Nt
    error('sampleTimes must have Nt=%d entries.', Nt);
end

lo = targetGamma * (1 - bandFraction);
hi = targetGamma * (1 + bandFraction);
validMask = true(Nx, Ny);
excludeWallCells = max(0, round(excludeWallCells));
if excludeWallCells > 0 && 2*excludeWallCells < Ny
    validMask(:, 1:excludeWallCells) = false;
    validMask(:, Ny-excludeWallCells+1:Ny) = false;
end
validIdx = validMask(:);

meanN = nan(Nt, 1);
stdN = nan(Nt, 1);
cvN = nan(Nt, 1);
rmsRel = nan(Nt, 1);
maxAbsRel = nan(Nt, 1);
meanAbsRel = nan(Nt, 1);
p05 = nan(Nt, 1);
p95 = nan(Nt, 1);
outBandFraction = nan(Nt, 1);
emptyFraction = nan(Nt, 1);
meanNExcl = nan(Nt, 1);
stdNExcl = nan(Nt, 1);
cvNExcl = nan(Nt, 1);
rmsRelExcl = nan(Nt, 1);
outBandFractionExcl = nan(Nt, 1);
lowKEnergy = nan(Nt, 1);
totalSpectralEnergy = nan(Nt, 1);
lowKFraction = nan(Nt, 1);

lowKMask = make_low_k_mask(Nx, Ny, lowKMaxIndex);
lowKMask(1, 1) = false;
nonzeroMask = true(Nx, Ny);
nonzeroMask(1, 1) = false;

for it = 1:Nt
    N = double(NMaps(:, :, it));
    n = N(:);
    rel = N / targetGamma - 1.0;
    r = rel(:);

    meanN(it) = mean(n, 'omitnan');
    stdN(it) = std(n, 'omitnan');
    cvN(it) = stdN(it) / max(meanN(it), eps);
    rmsRel(it) = sqrt(mean(r.^2, 'omitnan'));
    maxAbsRel(it) = max(abs(r));
    meanAbsRel(it) = mean(abs(r), 'omitnan');
    p05(it) = percentile_simple(n, 5);
    p95(it) = percentile_simple(n, 95);
    outBandFraction(it) = nnz(n < lo | n > hi) / numel(n);
    emptyFraction(it) = nnz(n == 0) / numel(n);

    nEx = n(validIdx);
    rEx = r(validIdx);
    meanNExcl(it) = mean(nEx, 'omitnan');
    stdNExcl(it) = std(nEx, 'omitnan');
    cvNExcl(it) = stdNExcl(it) / max(meanNExcl(it), eps);
    rmsRelExcl(it) = sqrt(mean(rEx.^2, 'omitnan'));
    outBandFractionExcl(it) = nnz(nEx < lo | nEx > hi) / numel(nEx);

    rhoHat = fft2(rel);
    energy = abs(rhoHat).^2 / numel(rhoHat)^2;
    lowKEnergy(it) = sum(energy(lowKMask));
    totalSpectralEnergy(it) = sum(energy(nonzeroMask));
    lowKFraction(it) = lowKEnergy(it) / max(totalSpectralEnergy(it), eps);
end

timeAvgN = mean(double(NMaps), 3);
finalN = double(NMaps(:, :, end));
initialN = double(NMaps(:, :, 1));
timeAvgRel = timeAvgN / targetGamma - 1.0;
finalRel = finalN / targetGamma - 1.0;

metrics = struct();
metrics.Nx = Nx;
metrics.Ny = Ny;
metrics.Nt = Nt;
metrics.sampleTimes = sampleTimes;
metrics.targetGamma = targetGamma;
metrics.bandFraction = bandFraction;
metrics.excludeWallCells = excludeWallCells;
metrics.lowKMaxIndex = lowKMaxIndex;
metrics.meanN = meanN;
metrics.stdN = stdN;
metrics.cvN = cvN;
metrics.rmsRel = rmsRel;
metrics.maxAbsRel = maxAbsRel;
metrics.meanAbsRel = meanAbsRel;
metrics.p05 = p05;
metrics.p95 = p95;
metrics.outBandFraction = outBandFraction;
metrics.emptyFraction = emptyFraction;
metrics.meanNExcl = meanNExcl;
metrics.stdNExcl = stdNExcl;
metrics.cvNExcl = cvNExcl;
metrics.rmsRelExcl = rmsRelExcl;
metrics.outBandFractionExcl = outBandFractionExcl;
metrics.lowKEnergy = lowKEnergy;
metrics.totalSpectralEnergy = totalSpectralEnergy;
metrics.lowKFraction = lowKFraction;
metrics.initialN = initialN;
metrics.finalN = finalN;
metrics.timeAvgN = timeAvgN;
metrics.timeAvgRel = timeAvgRel;
metrics.finalRel = finalRel;
metrics.validMask = validMask;

summary = struct();
summary.initialStdN = stdN(1);
summary.initialCvN = cvN(1);
summary.initialRmsRel = rmsRel(1);
summary.initialOutBandFraction = outBandFraction(1);
summary.initialEmptyFraction = emptyFraction(1);
summary.initialLowKEnergy = lowKEnergy(1);
summary.meanStdN = mean(stdN, 'omitnan');
summary.finalStdN = stdN(end);
summary.meanCvN = mean(cvN, 'omitnan');
summary.finalCvN = cvN(end);
summary.meanRmsRel = mean(rmsRel, 'omitnan');
summary.finalRmsRel = rmsRel(end);
summary.meanOutBandFraction = mean(outBandFraction, 'omitnan');
summary.finalOutBandFraction = outBandFraction(end);
summary.meanEmptyFraction = mean(emptyFraction, 'omitnan');
summary.finalEmptyFraction = emptyFraction(end);
summary.meanStdNExcl = mean(stdNExcl, 'omitnan');
summary.finalStdNExcl = stdNExcl(end);
summary.meanRmsRelExcl = mean(rmsRelExcl, 'omitnan');
summary.finalRmsRelExcl = rmsRelExcl(end);
summary.meanOutBandFractionExcl = mean(outBandFractionExcl, 'omitnan');
summary.finalOutBandFractionExcl = outBandFractionExcl(end);
summary.meanLowKEnergy = mean(lowKEnergy, 'omitnan');
summary.finalLowKEnergy = lowKEnergy(end);
summary.meanLowKFraction = mean(lowKFraction, 'omitnan');
summary.finalLowKFraction = lowKFraction(end);
summary.finalP05 = p05(end);
summary.finalP95 = p95(end);
summary.timeAvgRelRms = sqrt(mean(timeAvgRel(:).^2, 'omitnan'));
summary.finalRelRms = sqrt(mean(finalRel(:).^2, 'omitnan'));
metrics.summary = summary;
end

function mask = make_low_k_mask(Nx, Ny, kmax)
ix = [0:floor(Nx/2), -ceil(Nx/2)+1:-1];
iy = [0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix) > Nx
    ix = ix(1:Nx);
end
if numel(iy) > Ny
    iy = iy(1:Ny);
end
[KX, KY] = ndgrid(ix, iy);
mask = sqrt(double(KX).^2 + double(KY).^2) <= kmax;
end

function q = percentile_simple(x, p)
x = sort(x(:));
x = x(isfinite(x));
if isempty(x)
    q = NaN;
    return;
end
pos = 1 + (numel(x)-1) * p / 100;
i0 = floor(pos);
i1 = ceil(pos);
if i0 == i1
    q = x(i0);
else
    a = pos - i0;
    q = (1-a) * x(i0) + a * x(i1);
end
end
