function tr = projection_population_transport_diagnostics(x, vClassic, vProjected, params, varargin)
%PROJECTION_POPULATION_TRANSPORT_DIAGNOSTICS One-step density transport forecast.
%
%   tr = projection_population_transport_diagnostics(x, vClassic, vProjected, params)
%
% The pressure projection modifies velocities, not positions. Therefore the
% instantaneous population is unchanged by the projection substep. However,
% the corrected velocities change how particles will be transported during the
% next streaming substep. This diagnostic quantifies that effect without
% modifying the simulation state.
%
% Starting from the same particle positions x, it computes two kinematic
% one-step forecasts:
%   xClassicNext   = stream(x, vClassic)
%   xProjectedNext = stream(x, vProjected)
%
% and compares the resulting cell populations. It does not include SRD
% collision randomness; it is a deterministic diagnostic for the density
% transport induced by the velocity correction.
%
% Optional name-value arguments:
%   'includeBodyForce'   default true
%   'periodicX'          default true
%   'periodicY'          default false
%   'bandFraction'       default 0.20
%
% Output fields include popNow, popClassicForecast, popProjectedForecast,
% deltaClassic, deltaProjected, deltaProjectedMinusClassic, and compact scalar
% summaries such as rms/max/std/empty-cell changes.

if nargin < 4
    error('Usage: tr = projection_population_transport_diagnostics(x, vClassic, vProjected, params, ...)');
end
if size(x,2) ~= 2 || size(vClassic,2) ~= 2 || size(vProjected,2) ~= 2
    error('x, vClassic and vProjected must be Np-by-2 arrays.');
end
if size(x,1) ~= size(vClassic,1) || size(x,1) ~= size(vProjected,1)
    error('x, vClassic and vProjected must have the same particle count.');
end

includeBodyForce = true;
periodicX = true;
periodicY = false;
bandFraction = 0.20;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "includebodyforce"
            includeBodyForce = logical(val);
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        case "bandfraction"
            bandFraction = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

popNow = projection_population_diagnostics(x, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'bandFraction', bandFraction);

xClassicNext = forecast_stream_positions(x, vClassic, params, includeBodyForce, periodicX, periodicY);
xProjectedNext = forecast_stream_positions(x, vProjected, params, includeBodyForce, periodicX, periodicY);

popClassic = projection_population_diagnostics(xClassicNext, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'bandFraction', bandFraction);
popProjected = projection_population_diagnostics(xProjectedNext, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'bandFraction', bandFraction);

deltaClassic = double(popClassic.N) - double(popNow.N);
deltaProjected = double(popProjected.N) - double(popNow.N);
deltaProjectedMinusClassic = double(popProjected.N) - double(popClassic.N);

tr = struct();
tr.includeBodyForce = includeBodyForce;
tr.periodicX = periodicX;
tr.periodicY = periodicY;
tr.popNow = popNow;
tr.popClassicForecast = popClassic;
tr.popProjectedForecast = popProjected;
tr.deltaClassic = deltaClassic;
tr.deltaProjected = deltaProjected;
tr.deltaProjectedMinusClassic = deltaProjectedMinusClassic;
tr.classic = summarize_delta(deltaClassic, popNow, popClassic);
tr.projected = summarize_delta(deltaProjected, popNow, popProjected);
tr.projectedMinusClassic = summarize_delta(deltaProjectedMinusClassic, popClassic, popProjected);
tr.massErrorClassic = sum(popClassic.N(:)) - sum(popNow.N(:));
tr.massErrorProjected = sum(popProjected.N(:)) - sum(popNow.N(:));
tr.massErrorProjectedMinusClassic = sum(popProjected.N(:)) - sum(popClassic.N(:));
tr.yProfileNow = popNow.yProfile;
tr.yProfileClassicForecast = popClassic.yProfile;
tr.yProfileProjectedForecast = popProjected.yProfile;
tr.yProfileDeltaClassic = popClassic.yProfile - popNow.yProfile;
tr.yProfileDeltaProjected = popProjected.yProfile - popNow.yProfile;
tr.yProfileProjectedMinusClassic = popProjected.yProfile - popClassic.yProfile;
end

function xNext = forecast_stream_positions(x, v, params, includeBodyForce, periodicX, periodicY)
Lx = params.Lx;
Ly = params.Ly;
dt = params.dt;
bodyForceX = get_param(params, 'bodyForceX', 0.0);
bodyForceY = get_param(params, 'bodyForceY', 0.0);

vf = v;
if includeBodyForce
    vf(:,1) = vf(:,1) + dt * bodyForceX;
    vf(:,2) = vf(:,2) + dt * bodyForceY;
end

xNext = x + dt * vf;
if periodicX
    xNext(:,1) = mod(xNext(:,1), Lx);
else
    xNext(:,1) = min(max(xNext(:,1), 0), Lx - eps(Lx));
end

if periodicY
    xNext(:,2) = mod(xNext(:,2), Ly);
else
    % Specular position reflection only. For population forecasts we need
    % positions, not wall momentum exchange. The triangular-wave map handles
    % occasional displacements larger than Ly robustly.
    y = mod(xNext(:,2), 2*Ly);
    high = y > Ly;
    y(high) = 2*Ly - y(high);
    xNext(:,2) = min(max(y, 0), Ly - eps(Ly));
end
end

function s = summarize_delta(deltaN, popA, popB)
d = double(deltaN(:));
s = struct();
s.rms = sqrt(mean(d.^2));
s.maxAbs = max(abs(d));
s.meanAbs = mean(abs(d));
s.std = std(d);
s.sum = sum(d);
s.rmsOverGamma = s.rms / max(popA.targetGamma, eps);
s.maxAbsOverGamma = s.maxAbs / max(popA.targetGamma, eps);
s.stdBefore = popA.stdN;
s.stdAfter = popB.stdN;
s.stdDelta = popB.stdN - popA.stdN;
s.emptyBefore = popA.nEmptyCells;
s.emptyAfter = popB.nEmptyCells;
s.emptyDelta = popB.nEmptyCells - popA.nEmptyCells;
s.outBandBefore = popA.outBandFraction;
s.outBandAfter = popB.outBandFraction;
s.outBandDelta = popB.outBandFraction - popA.outBandFraction;
s.cvBefore = popA.coefficientOfVariation;
s.cvAfter = popB.coefficientOfVariation;
s.cvDelta = popB.coefficientOfVariation - popA.coefficientOfVariation;
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
