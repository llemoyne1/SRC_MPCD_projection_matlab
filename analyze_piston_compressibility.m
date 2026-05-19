function analysis = analyze_piston_compressibility(out, varargin)
%ANALYZE_PISTON_COMPRESSIBILITY Fit effective compressibility diagnostics.
%
%   analysis = analyze_piston_compressibility(out)
%
% Uses the piston diagnostic table to fit pressure-like quantities versus
% rhoPhysicalMean.  For a linear fit P = a*rho + b, the reported effective
% bulk modulus is K_eff = rhoFitMean * a.  With the active thermostat this is
% primarily an isothermal numerical compressibility diagnostic.
%
% Name-value options:
%   fitStartFraction  default 0.2
%   fitEndFraction    default 1.0
%   pressureFields    default {'PkinMean','PkinTopLayerMean','PvirMean','PtotMean','PdriveMean','pressureTopWallTotal'}

opts = parse_inputs(varargin{:});
opts.fitStartFraction = getf(opts, 'fitStartFraction', 0.2);
opts.fitEndFraction = getf(opts, 'fitEndFraction', 1.0);
opts.pressureFields = getf(opts, 'pressureFields', {'PkinMean','PkinTopLayerMean','PvirMean','PtotMean','PdriveMean','pressureTopWallTotal'});

analysis = struct();
analysis.ok = false;
analysis.errorMessage = '';
analysis.fitStartFraction = opts.fitStartFraction;
analysis.fitEndFraction = opts.fitEndFraction;

if ~isstruct(out) || ~isfield(out, 'diagTable') || isempty(out.diagTable) || height(out.diagTable) < 3
    analysis.errorMessage = 'No usable diagTable.';
    analysis = fill_defaults(analysis);
    return;
end
T = out.diagTable;
if ~ismember('rhoPhysicalMean', T.Properties.VariableNames)
    analysis.errorMessage = 'diagTable lacks rhoPhysicalMean.';
    analysis = fill_defaults(analysis);
    return;
end

n = height(T);
i0 = max(1, min(n, floor(opts.fitStartFraction*n) + 1));
i1 = max(i0, min(n, ceil(opts.fitEndFraction*n)));
idx = false(n,1);
idx(i0:i1) = true;
rho = T.rhoPhysicalMean;
idx = idx & isfinite(rho);
analysis.fitIndexStart = i0;
analysis.fitIndexEnd = i1;
analysis.tFitStart = safe_get(T, 't', i0);
analysis.tFitEnd = safe_get(T, 't', i1);
analysis.nFitSamples = nnz(idx);
analysis.rhoFitMean = mean(rho(idx), 'omitnan');
analysis.rhoFitMin = min(rho(idx));
analysis.rhoFitMax = max(rho(idx));

fields = opts.pressureFields;
for k = 1:numel(fields)
    fld = char(string(fields{k}));
    if ismember(fld, T.Properties.VariableNames)
        P = T.(fld);
    else
        P = nan(n,1);
    end
    fr = fit_one(rho, P, idx);
    safeName = matlab.lang.makeValidName(fld);
    analysis.([safeName 'Slope_dPdrho']) = fr.slope;
    analysis.(['Keff_' safeName]) = fr.Keff;
    analysis.(['R2_' safeName]) = fr.R2;
    analysis.(['PMean_' safeName]) = fr.PMean;
    analysis.(['PMin_' safeName]) = fr.PMin;
    analysis.(['PMax_' safeName]) = fr.PMax;
end

% Short aliases used by existing runners and console prints.
analysis.KeffPkin = getf(analysis, 'Keff_PkinMean', NaN);
analysis.R2Pkin = getf(analysis, 'R2_PkinMean', NaN);
analysis.KeffPkinTopLayer = getf(analysis, 'Keff_PkinTopLayerMean', NaN);
analysis.R2PkinTopLayer = getf(analysis, 'R2_PkinTopLayerMean', NaN);
analysis.KeffTopWallTotal = getf(analysis, 'Keff_pressureTopWallTotal', NaN);
analysis.R2TopWallTotal = getf(analysis, 'R2_pressureTopWallTotal', NaN);
analysis.KeffPvir = getf(analysis, 'Keff_PvirMean', NaN);
analysis.R2Pvir = getf(analysis, 'R2_PvirMean', NaN);
analysis.KeffPtot = getf(analysis, 'Keff_PtotMean', NaN);
analysis.R2Ptot = getf(analysis, 'R2_PtotMean', NaN);
analysis.KeffPdrive = getf(analysis, 'Keff_PdriveMean', NaN);
analysis.R2Pdrive = getf(analysis, 'R2_PdriveMean', NaN);
analysis.ok = true;
end

function fr = fit_one(rho, P, idx)
fr = struct('slope', NaN, 'intercept', NaN, 'Keff', NaN, 'R2', NaN, ...
            'PMean', NaN, 'PMin', NaN, 'PMax', NaN);
idx = idx & isfinite(rho) & isfinite(P);
if nnz(idx) < 3 || max(rho(idx)) <= min(rho(idx))
    return;
end
x = rho(idx);
y = P(idx);
coef = polyfit(x, y, 1);
yfit = polyval(coef, x);
ssRes = sum((y-yfit).^2, 'omitnan');
ssTot = sum((y-mean(y,'omitnan')).^2, 'omitnan');
fr.slope = coef(1);
fr.intercept = coef(2);
fr.Keff = mean(x, 'omitnan') * coef(1);
if ssTot > 0
    fr.R2 = 1 - ssRes/ssTot;
end
fr.PMean = mean(y, 'omitnan');
fr.PMin = min(y);
fr.PMax = max(y);
end

function analysis = fill_defaults(analysis)
analysis.KeffPkin = NaN;
analysis.R2Pkin = NaN;
analysis.KeffPkinTopLayer = NaN;
analysis.R2PkinTopLayer = NaN;
analysis.KeffTopWallTotal = NaN;
analysis.R2TopWallTotal = NaN;
analysis.KeffPvir = NaN;
analysis.R2Pvir = NaN;
analysis.KeffPtot = NaN;
analysis.R2Ptot = NaN;
analysis.KeffPdrive = NaN;
analysis.R2Pdrive = NaN;
end

function val = safe_get(T, fld, idx)
if ismember(fld, T.Properties.VariableNames) && idx >= 1 && idx <= height(T)
    val = T.(fld)(idx);
else
    val = NaN;
end
end

function opts = parse_inputs(varargin)
opts = struct();
if nargin == 1 && isstruct(varargin{1})
    opts = varargin{1};
    return;
end
if mod(nargin, 2) ~= 0
    error('Use name-value pairs or a single struct.');
end
for i = 1:2:nargin
    opts.(char(varargin{i})) = varargin{i+1};
end
end

function value = getf(s, name, defaultValue)
if nargin < 3
    defaultValue = NaN;
end
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
