function shed = projection_wake_shedding_diagnostics(t, signal, params, Uref)
%PROJECTION_WAKE_SHEDDING_DIAGNOSTICS Estimate shedding frequency and Strouhal.
%
%   shed = projection_wake_shedding_diagnostics(t, signal, params, Uref)
%
% The signal should be an off-center wake probe, typically Uy or omega behind
% the cylinder. The diagnostic discards an initial transient fraction, removes
% the mean and a linear trend, applies a Hann window, then picks the dominant
% positive FFT peak.
%
% Returned fields include frequency, Strouhal number St = f D / Uref,
% dominant period, peak amplitude, and a simple peak-to-median spectral SNR.
% Optionally, the peak search can be restricted to a plausible Strouhal band
% to avoid selecting high-frequency collision noise on long particle runs.

if nargin < 4 || isempty(Uref)
    Uref = NaN;
end

t = t(:);
signal = signal(:);
valid = isfinite(t) & isfinite(signal);
t = t(valid);
signal = signal(valid);

shed = empty_shed();
if numel(t) < 8
    shed.note = 'not enough finite samples';
    return;
end

transientFraction = get_param(params, 'sheddingTransientFraction', 0.40);
transientFraction = min(max(transientFraction, 0), 0.95);
t0 = min(t) + transientFraction * (max(t) - min(t));
keep = t >= t0;
t = t(keep);
signal = signal(keep);
if numel(t) < 8
    shed.note = 'not enough samples after transient cut';
    return;
end

% Require near-uniform sampling; if not exactly uniform, interpolate to a
% uniform grid based on the median time step.
dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    shed.note = 'invalid sample time step';
    return;
end
tu = (t(1):dt:t(end)).';
if numel(tu) < 8
    shed.note = 'uniform resampling too short';
    return;
end
su = interp1(t, signal, tu, 'linear', 'extrap');

su = su - mean(su, 'omitnan');
if numel(su) >= 3
    p = polyfit(tu - tu(1), su, 1);
    su = su - polyval(p, tu - tu(1));
end

N = numel(su);
if N < 8 || all(abs(su) < eps)
    shed.note = 'flat signal';
    return;
end

window = hann_local(N);
sw = su .* window;
Y = fft(sw);
P = abs(Y / N).^2;
f = (0:N-1).' / (N * dt);
positive = 2:floor(N/2); % skip zero frequency
if isempty(positive)
    shed.note = 'no positive FFT bins';
    return;
end

PposFull = P(positive);
fposFull = f(positive);
D = 2 * get_param(params, 'cylinderRadius', NaN);

useBand = logical(get_param(params, 'sheddingUseStrouhalBand', true));
StMin = get_param(params, 'sheddingStrouhalMin', 0.05);
StMax = get_param(params, 'sheddingStrouhalMax', 0.50);
if useBand && isfinite(Uref) && abs(Uref) > eps && isfinite(D) && D > 0
    StFull = fposFull * D / abs(Uref);
    band = isfinite(StFull) & StFull >= StMin & StFull <= StMax;
else
    band = true(size(fposFull));
end

if ~any(band)
    shed.note = 'no FFT bins inside Strouhal search band';
    shed.freq = fposFull;
    shed.power = PposFull;
    shed.D = D;
    shed.Uref = Uref;
    shed.searchStrouhalMin = StMin;
    shed.searchStrouhalMax = StMax;
    shed.usedBandLimitedSearch = useBand;
    return;
end

Ppos = PposFull(band);
fpos = fposFull(band);
[peakPower, idx] = max(Ppos);
fPeak = fpos(idx);

if isfinite(Uref) && abs(Uref) > eps && isfinite(D) && D > 0
    St = fPeak * D / abs(Uref);
else
    St = NaN;
end

medianPower = median(Ppos, 'omitnan');
if ~isfinite(medianPower) || medianPower <= 0
    snr = NaN;
else
    snr = peakPower / medianPower;
end

shed.valid = true;
shed.note = 'ok';
shed.frequency = fPeak;
shed.period = 1 / fPeak;
shed.strouhal = St;
shed.peakPower = peakPower;
shed.spectralSNR = snr;
shed.signalRms = sqrt(mean(su.^2, 'omitnan'));
shed.nSamplesUsed = N;
shed.tStartUsed = tu(1);
shed.tEndUsed = tu(end);
shed.dtSample = dt;
shed.Uref = Uref;
shed.D = D;
shed.freq = fpos;
shed.power = Ppos;
shed.freqFull = fposFull;
shed.powerFull = PposFull;
shed.searchStrouhalMin = StMin;
shed.searchStrouhalMax = StMax;
shed.usedBandLimitedSearch = useBand;
shed.nSearchBins = numel(fpos);
end

function shed = empty_shed()
shed = struct();
shed.valid = false;
shed.note = '';
shed.frequency = NaN;
shed.period = NaN;
shed.strouhal = NaN;
shed.peakPower = NaN;
shed.spectralSNR = NaN;
shed.signalRms = NaN;
shed.nSamplesUsed = 0;
shed.tStartUsed = NaN;
shed.tEndUsed = NaN;
shed.dtSample = NaN;
shed.Uref = NaN;
shed.D = NaN;
shed.freq = [];
shed.power = [];
shed.freqFull = [];
shed.powerFull = [];
shed.searchStrouhalMin = NaN;
shed.searchStrouhalMax = NaN;
shed.usedBandLimitedSearch = false;
shed.nSearchBins = 0;
end

function w = hann_local(N)
if N <= 1
    w = ones(N, 1);
else
    n = (0:N-1).';
    w = 0.5 - 0.5 * cos(2*pi*n/(N-1));
end
end

function value = get_param(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
