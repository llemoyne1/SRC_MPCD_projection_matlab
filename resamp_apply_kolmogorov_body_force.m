function [stateOut, info] = resamp_apply_kolmogorov_body_force(state, params, varargin)
%RESAMP_APPLY_KOLMOGOROV_BODY_FORCE Apply a mean-free periodic Kolmogorov force.
%
%   [stateOut, info] = resamp_apply_kolmogorov_body_force(state, params, ...)
%
% The default forcing is
%
%       a_x(y) = A sin(2*pi*k*y/Ly + phase),  a_y = 0,
%
% applied to active particles as v <- v + dt*a.  Because the particle
% sampling is finite and weighted, the function can subtract the exact
% mass-weighted mean kick so that the force does not accelerate the periodic
% box as a whole.
%
% Options:
%   'amplitude'    : acceleration amplitude A, default 0
%   'waveNumber'   : integer/real k, default 1
%   'phase'        : phase in radians, default 0
%   'direction'    : 'x' or 'y', default 'x'
%   'zeroMeanKick' : subtract mass-weighted mean kick, default true
%   'applyDt'      : multiply acceleration by params.dt, default true
%   'activeMask'   : optional logical active-particle mask

opts = parse_options(varargin{:});
stateOut = state;
active = opts.activeMask;
if isempty(active)
    active = resamp_active_mask(state);
end
active = active(:);
idx = find(active);

info = empty_info(opts);
info.nActive = numel(idx);
if isempty(idx) || opts.amplitude == 0
    info.enabled = false;
    return;
end

if ~isfield(params, 'Ly') || isempty(params.Ly) || ~isfield(params, 'dt') || isempty(params.dt)
    error('params must contain Ly and dt for Kolmogorov forcing.');
end
Ly = params.Ly;
dt = params.dt;

x = state.x(idx,:);
v = state.v(idx,:);
m = state.m(idx);
if any(~isfinite(m)) || any(m < 0)
    error('state.m contains invalid active masses.');
end
mtot = sum(m, 'omitnan');
if ~(isfinite(mtot) && mtot > 0)
    info.enabled = false;
    return;
end

phaseArg = 2*pi*opts.waveNumber*mod(x(:,2), Ly)/Ly + opts.phase;
a = opts.amplitude * sin(phaseArg);
if opts.applyDt
    kickScalar = dt * a;
else
    kickScalar = a;
end

kick = zeros(numel(idx), 2);
switch lower(char(string(opts.direction)))
    case 'x'
        kick(:,1) = kickScalar;
    case 'y'
        kick(:,2) = kickScalar;
    otherwise
        error('Kolmogorov force direction must be ''x'' or ''y'', got %s.', opts.direction);
end

rawMeanKick = [sum(m .* kick(:,1), 'omitnan'), sum(m .* kick(:,2), 'omitnan')] / mtot;
if opts.zeroMeanKick
    kick(:,1) = kick(:,1) - rawMeanKick(1);
    kick(:,2) = kick(:,2) - rawMeanKick(2);
end
meanKick = [sum(m .* kick(:,1), 'omitnan'), sum(m .* kick(:,2), 'omitnan')] / mtot;

momentumBefore = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];
v = v + kick;
momentumAfter = [sum(m .* v(:,1), 'omitnan'), sum(m .* v(:,2), 'omitnan')];

stateOut.v(idx,:) = v;

kickMag = sqrt(sum(kick.^2, 2));
info.enabled = true;
info.mode = 'kolmogorov_x';
if strcmpi(opts.direction, 'y')
    info.mode = 'kolmogorov_y';
end
info.amplitude = opts.amplitude;
info.waveNumber = opts.waveNumber;
info.phase = opts.phase;
info.direction = char(string(opts.direction));
info.zeroMeanKick = opts.zeroMeanKick;
info.applyDt = opts.applyDt;
info.rawMeanKick = rawMeanKick;
info.meanKick = meanKick;
info.kickRms = sqrt(mean(kickMag.^2, 'omitnan'));
info.kickMax = max(kickMag, [], 'omitnan');
info.momentumBefore = momentumBefore;
info.momentumAfter = momentumAfter;
info.momentumDelta = momentumAfter - momentumBefore;
info.momentumDeltaNorm = norm(info.momentumDelta);
info.accelerationRms = info.kickRms / max(dt, eps);
info.accelerationMax = info.kickMax / max(dt, eps);
end

function opts = parse_options(varargin)
opts = struct();
opts.amplitude = 0.0;
opts.waveNumber = 1;
opts.phase = 0.0;
opts.direction = 'x';
opts.zeroMeanKick = true;
opts.applyDt = true;
opts.activeMask = [];
if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
        case {'amplitude','forceamplitude','kolmogorovforceamplitude'}
            opts.amplitude = val;
        case {'wavenumber','modenumber','k','kolmogorovforcewavenumber'}
            opts.waveNumber = val;
        case {'phase','kolmogorovforcephase'}
            opts.phase = val;
        case {'direction','forcedirection','kolmogorovforcedirection'}
            opts.direction = lower(char(string(val)));
        case {'zeromeankick','forcemeanmomentumcorrection','kolmogorovforcezeromeankick'}
            opts.zeroMeanKick = logical(val);
        case {'applydt','kolmogorovforceapplydt'}
            opts.applyDt = logical(val);
        case 'activemask'
            opts.activeMask = logical(val(:));
        otherwise
            error('Unknown option for resamp_apply_kolmogorov_body_force: %s', key);
    end
end
end

function info = empty_info(opts)
info = struct();
info.enabled = false;
info.mode = 'none';
info.amplitude = opts.amplitude;
info.waveNumber = opts.waveNumber;
info.phase = opts.phase;
info.direction = char(string(opts.direction));
info.zeroMeanKick = opts.zeroMeanKick;
info.applyDt = opts.applyDt;
info.nActive = 0;
info.rawMeanKick = [0 0];
info.meanKick = [0 0];
info.kickRms = 0;
info.kickMax = 0;
info.momentumBefore = [NaN NaN];
info.momentumAfter = [NaN NaN];
info.momentumDelta = [0 0];
info.momentumDeltaNorm = 0;
info.accelerationRms = 0;
info.accelerationMax = 0;
end
