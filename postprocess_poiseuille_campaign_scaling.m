function summaryTable = postprocess_poiseuille_campaign_scaling(campaignInput, varargin)
%POSTPROCESS_POISEUILLE_CAMPAIGN_SCALING Recompute raw/scaled Poiseuille viscosities.
%
%   summaryTable = postprocess_poiseuille_campaign_scaling('campaign_results.mat')
%   summaryTable = postprocess_poiseuille_campaign_scaling(results)
%
% Re-analyzes saved Poiseuille outputs with the current
% analyze_projection_poiseuille_viscosity.m and writes a compact CSV with
% nuEffRaw, nuEffScaledNyRef, nuEffCellY, dy and acceleration diagnostics.
%
% Useful for campaign_results.mat files generated before the scaled
% diagnostics were added.

p = inputParser;
addRequired(p, 'campaignInput');
addParameter(p, 'outputCsv', '');
addParameter(p, 'nuReferenceNy', 64);
addParameter(p, 'nuScalingEnable', true);
addParameter(p, 'fitWindowTime', []);
addParameter(p, 'excludeWallCells', []);
addParameter(p, 'accelerationWindowTime', 5);
parse(p, campaignInput, varargin{:});
opt = p.Results;

results = load_results(campaignInput);
outs = results.outputs;
fields = fieldnames(outs);
rows = repmat(empty_row(), numel(fields), 1);

for i = 1:numel(fields)
    label = fields{i};
    out = outs.(label);
    if isfield(out, 'out')
        out = out.out;
    elseif isfield(out, 'outLight')
        out = out.outLight;
    end
    if ~isstruct(out) || ~isfield(out, 'params')
        rows(i).label = {label};
        rows(i).ok = false;
        rows(i).errorMessage = {'No output struct found'};
        continue;
    end

    out.params.poiseuilleNuReferenceNy = opt.nuReferenceNy;
    out.params.poiseuilleNuScalingEnable = opt.nuScalingEnable;
    if ~isfield(out.params, 'poiseuillePhysicalLy') || isempty(out.params.poiseuillePhysicalLy)
        out.params.poiseuillePhysicalLy = get_field(out.params, 'Ly', 1.0);
    end

    fitArgs = {};
    if ~isempty(opt.excludeWallCells)
        fitArgs = [fitArgs, {'excludeWallCells', opt.excludeWallCells}]; %#ok<AGROW>
    elseif isfield(out.params, 'excludeWallCellsFit') && ~isempty(out.params.excludeWallCellsFit)
        fitArgs = [fitArgs, {'excludeWallCells', out.params.excludeWallCellsFit}]; %#ok<AGROW>
    end
    if ~isempty(opt.fitWindowTime)
        fitArgs = [fitArgs, {'fitWindowTime', opt.fitWindowTime}]; %#ok<AGROW>
    elseif isfield(out.params, 'fitWindowTime') && ~isempty(out.params.fitWindowTime)
        fitArgs = [fitArgs, {'fitWindowTime', out.params.fitWindowTime}]; %#ok<AGROW>
    elseif isfield(out.params, 'fitStartFraction') && ~isempty(out.params.fitStartFraction)
        fitArgs = [fitArgs, {'fitStartFraction', out.params.fitStartFraction}]; %#ok<AGROW>
    end
    fitArgs = [fitArgs, {'nuReferenceNy', opt.nuReferenceNy, 'nuScalingEnable', opt.nuScalingEnable}]; %#ok<AGROW>

    try
        visc = analyze_projection_poiseuille_viscosity(out, fitArgs{:});
        acc = poiseuille_acceleration_diagnostics(out, opt.accelerationWindowTime);

        row = empty_row();
        row.label = {label};
        row.method = {char(string(get_field(out.params, 'method', '')))};
        row.ok = true;
        row.Nx = get_field(out.params, 'Nx', NaN);
        row.Ny = get_field(out.params, 'Ny', NaN);
        row.gamma = get_field(out.params, 'gamma', NaN);
        row.nSteps = get_field(out.params, 'nSteps', NaN);
        row.actualStep = get_field(out, 'actualLastStep', NaN);
        row.nuEffRaw = visc.nuEffRaw;
        row.nuEffScaledNyRef = visc.nuEffScaledNyRef;
        row.nuEffCellY = visc.nuEffCellY;
        row.nuScaleFactorNyRef = visc.nuScaleFactorNyRef;
        row.nuReferenceNy = visc.nuReferenceNy;
        row.dy = visc.dy;
        row.R2 = visc.R2;
        row.SNR = visc.signalToNoise;
        row.fitT0 = visc.tMin;
        row.fitT1 = visc.tMax;
        row.centerWallFit = visc.centerMinusWall;
        row.ratioAccelRecent = acc.ratioAccel;
        row.slopeRecent = acc.slopeMeanUx;
        row.accelT0 = acc.t0;
        row.accelT1 = acc.t1;
        if isfield(out, 'diagTable') && istable(out.diagTable) && height(out.diagTable) > 0
            row.finalMeanUx = last_value(out.diagTable, 'meanUx');
            row.finalCenterWall = last_value(out.diagTable, 'centerMinusWall');
            row.finalLowK = last_value(out.diagTable, 'lowKDensityEnergy');
            row.finalKBT = last_value(out.diagTable, 'kBTCell');
        end
        rows(i) = row;
    catch ME
        rows(i).label = {label};
        rows(i).ok = false;
        rows(i).errorMessage = {ME.message};
    end
end

summaryTable = struct2table(rows);

if isempty(opt.outputCsv)
    if ischar(campaignInput) || isstring(campaignInput)
        [folder,~,~] = fileparts(char(string(campaignInput)));
        if isempty(folder), folder = pwd; end
        opt.outputCsv = fullfile(folder, 'summary_scaled.csv');
    else
        opt.outputCsv = fullfile(pwd, 'summary_scaled.csv');
    end
end
writetable(summaryTable, opt.outputCsv);
fprintf('Scaled Poiseuille summary written to: %s\n', opt.outputCsv);
disp(summaryTable);
end

function results = load_results(x)
if ischar(x) || isstring(x)
    S = load(char(string(x)));
    if isfield(S, 'results')
        results = S.results;
    else
        error('MAT file does not contain variable ''results''.');
    end
elseif isstruct(x)
    results = x;
else
    error('campaignInput must be a MAT filename or a results struct.');
end
if ~isfield(results, 'outputs')
    error('Input results struct has no outputs field.');
end
end

function row = empty_row()
row = struct('label', {{''}}, 'method', {{''}}, 'ok', false, ...
    'Nx', NaN, 'Ny', NaN, 'gamma', NaN, 'nSteps', NaN, 'actualStep', NaN, ...
    'nuEffRaw', NaN, 'nuEffScaledNyRef', NaN, 'nuEffCellY', NaN, ...
    'nuScaleFactorNyRef', NaN, 'nuReferenceNy', NaN, 'dy', NaN, ...
    'R2', NaN, 'SNR', NaN, 'fitT0', NaN, 'fitT1', NaN, 'centerWallFit', NaN, ...
    'slopeRecent', NaN, 'ratioAccelRecent', NaN, 'accelT0', NaN, 'accelT1', NaN, ...
    'finalMeanUx', NaN, 'finalCenterWall', NaN, 'finalLowK', NaN, 'finalKBT', NaN, ...
    'errorMessage', {{''}});
end

function val = last_value(T, name)
val = NaN;
if istable(T) && ismember(name, T.Properties.VariableNames) && height(T) >= 1
    x = T.(name);
    val = x(end);
end
end

function val = get_field(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
end
end
