function results = run_wallvp_equilibrium_pressure_validation(varargin)
%RUN_WALLVP_EQUILIBRIUM_PRESSURE_VALIDATION Stationary-wall pressure sanity test.
%
%   results = run_wallvp_equilibrium_pressure_validation(...)
%
% This is the first diagnostic to run before returning to the piston
% staircase.  It removes piston motion and virial forcing, compares wallVP off
% and wallVP on, and uses cumulative wall impulses rather than sparse
% instantaneous samples.  In equilibrium the total wall pressure
% impact+wallVP should be close to rho*kBT.  The wallVP contribution alone is
% expected to average near zero or remain small compared with the total for a
% stationary isothermal gas; it is mainly a no-slip/thermal wall collision
% model, not the complete ideal pressure estimator by itself.

opts = parse_inputs(varargin{:});
opts = set_default(opts, 'Nx', 32);
opts = set_default(opts, 'Ny', 32);
opts = set_default(opts, 'gamma', 20);
opts = set_default(opts, 'Lx', 1.0);
opts = set_default(opts, 'Ly', 1.0);
opts = set_default(opts, 'dt', 1.0e-3);
opts = set_default(opts, 'kBT', 0.01);
opts = set_default(opts, 'nSteps', 20000);
opts = set_default(opts, 'sampleEvery', 100);
opts = set_default(opts, 'progressEvery', 1000);
opts = set_default(opts, 'seed', 11);
opts = set_default(opts, 'densityFactor', 0.8);
opts = set_default(opts, 'wallVirtualParticlesGeometryMode', 'shifted_solid_fraction');
opts = set_default(opts, 'wallVirtualParticlesStochasticCount', true);
opts = set_default(opts, 'visualEnable', false);
opts = set_default(opts, 'saveOutput', false);
opts = set_default(opts, 'outputRoot', fullfile(pwd, ['wallvp_equilibrium_' datestr(now,'yyyymmdd_HHMMSS')]));

cases = struct([]);
cases(1).name = 'wallVP_off';
cases(1).enable = false;
cases(2).name = 'wallVP_on';
cases(2).enable = true;

rows = cell(numel(cases), 1);
outputs = struct();
for ic = 1:numel(cases)
    p = struct();
    p.method = 'classic';
    p.Nx = opts.Nx;
    p.Ny = opts.Ny;
    p.gamma = opts.gamma;
    p.Lx = opts.Lx;
    p.Ly = opts.Ly;
    p.Ly0 = opts.Ly;
    p.pistonY0 = opts.Ly;
    p.pistonYMin = opts.Ly;
    p.pistonVy = 0.0;
    p.pistonMotionMode = 'linear';
    p.compressionTarget = 0.0;
    p.pistonCompressionTarget = 0.0;
    p.dt = opts.dt;
    p.kBT = opts.kBT;
    p.nSteps = opts.nSteps;
    p.sampleEvery = opts.sampleEvery;
    p.progressEvery = opts.progressEvery;
    p.seed = opts.seed;
    p.initialPopulationMode = 'exact_per_cell';
    p.wallModeY = 'bounceback';
    p.pistonTangentialMode = 'bounceback';
    p.wallVirtualParticlesEnable = cases(ic).enable;
    p.wallVirtualParticlesGeometryMode = opts.wallVirtualParticlesGeometryMode;
    p.wallVirtualParticlesDensityFactor = opts.densityFactor;
    p.wallVirtualParticlesThermal = true;
    p.wallVirtualParticlesKBT = opts.kBT;
    p.wallVirtualParticlesForceRandomShiftY = true;
    p.wallVirtualParticlesStochasticCount = opts.wallVirtualParticlesStochasticCount;
    p.virialDiagnosticsEnable = false;
    p.virialKickEnable = false;
    p.visualEnable = opts.visualEnable;
    p.storeDensityMaps = false;
    p.computeFullDiagnosticsEveryStep = false;

    fprintf('\n=== wallVP equilibrium case: %s ===\n', cases(ic).name);
    out = run_projection_piston_wallvp_demo(p);
    outputs.(matlab.lang.makeValidName(cases(ic).name)) = out;
    rows{ic} = make_summary_row(cases(ic).name, out);
end

summary = vertcat(rows{:});
results = struct();
results.opts = opts;
results.summary = summary;
results.outputs = outputs;

fprintf('\n=== wallVP equilibrium pressure summary ===\n');
disp(summary);

if opts.saveOutput
    if ~exist(opts.outputRoot, 'dir')
        mkdir(opts.outputRoot);
    end
    writetable(summary, fullfile(opts.outputRoot, 'wallvp_equilibrium_summary.csv'));
    save(fullfile(opts.outputRoot, 'wallvp_equilibrium_results.mat'), 'results', '-v7.3');
end
end

function row = make_summary_row(caseName, out)
s = out.summary;
Pideal = s.finalRhoPhysicalMean * s.finalKBTCell;
PtopTotal = getf(s, 'finalTopWallPressureTotalCumulativeMean', NaN);
PbotTotal = getf(s, 'finalBottomWallPressureTotalCumulativeMean', NaN);
PtopImpact = getf(s, 'finalTopWallPressureImpactCumulativeMean', NaN);
PtopVP = getf(s, 'finalTopWallPressureVPCumulativeMean', NaN);
PbotImpact = getf(s, 'finalBottomWallPressureImpactCumulativeMean', NaN);
PbotVP = getf(s, 'finalBottomWallPressureVPCumulativeMean', NaN);
row = table(string(caseName), out.actualStep, s.finalRhoPhysicalMean, s.finalKBTCell, Pideal, ...
    PtopImpact, PtopVP, PtopTotal, PtopTotal / max(Pideal, eps), ...
    PbotImpact, PbotVP, PbotTotal, PbotTotal / max(Pideal, eps), ...
    abs(PtopTotal - PbotTotal) / max(Pideal, eps), ...
    'VariableNames', {'caseName','actualStep','rho','kBTCell','Pideal', ...
    'PtopImpactCum','PtopVPCum','PtopTotalCum','PtopTotalOverIdeal', ...
    'PbottomImpactCum','PbottomVPCum','PbottomTotalCum','PbottomTotalOverIdeal', ...
    'topBottomRelativeGap'});
end

function opts = parse_inputs(varargin)
if nargin == 1 && isstruct(varargin{1})
    opts = varargin{1};
    return;
end
opts = struct();
if mod(nargin, 2) ~= 0
    error('Use name-value pairs or a single struct.');
end
for k = 1:2:nargin
    opts.(char(varargin{k})) = varargin{k+1};
end
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
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
