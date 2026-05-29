function out = run_resamp_fluid_physics_validation_suite(varargin)
%RUN_RESAMP_FLUID_PHYSICS_VALIDATION_SUITE Launch TG and Poiseuille validations.
%
% This wrapper is intentionally small.  It first runs the fully periodic
% forced Taylor--Green validation, then the bounded-y Poiseuille validation.
%
% Example:
%   out = run_resamp_fluid_physics_validation_suite('quick',true);

opts = parse_options(varargin{:});
if opts.quick
    tgSteps = min(opts.tgSteps, 1000);
    poiSteps = min(opts.poiseuilleSteps, 1000);
    tgSummaryEvery = min(opts.tgSummaryEvery, 50);
    poiSummaryEvery = min(opts.poiseuilleSummaryEvery, 100);
else
    tgSteps = opts.tgSteps;
    poiSteps = opts.poiseuilleSteps;
    tgSummaryEvery = opts.tgSummaryEvery;
    poiSummaryEvery = opts.poiseuilleSummaryEvery;
end
root = opts.outputRoot;
if ~exist(root,'dir'), mkdir(root); end
outTG = run_resamp_tg_physics_validation( ...
    'outputRoot', fullfile(root,'tg'), ...
    'Nx', opts.NxTG, 'Ny', opts.NyTG, 'gamma', opts.gamma, ...
    'steps', tgSteps, 'summaryEvery', tgSummaryEvery, 'visualEvery', opts.visualEvery, ...
    'dt', opts.dtTG, 'kBT', opts.kBT, ...
    'NMin', opts.NMin, 'NTarget', opts.gamma, 'NMax', opts.NMax, ...
    'ThermostatStrength', opts.thermostatStrength, ...
    'TaylorGreenForcingAmplitude', opts.tgForceAmplitude, ...
    'rngSeed', opts.rngSeed);
outPoi = run_resamp_poiseuille_physics_validation( ...
    'outputRoot', fullfile(root,'poiseuille'), ...
    'Nx', opts.NxPoi, 'Ny', opts.NyPoi, 'gamma', opts.gamma, ...
    'steps', poiSteps, 'sampleEvery', opts.poiseuilleSampleEvery, 'summaryEvery', poiSummaryEvery, ...
    'dt', opts.dtPoi, 'kBT', opts.kBT, 'bodyForceX', opts.bodyForceX, ...
    'NMin', opts.NMin, 'NTarget', opts.gamma, 'NMax', opts.NMax, ...
    'ThermostatStrength', opts.thermostatStrength, ...
    'rngSeed', opts.rngSeed);
out = struct('options',opts,'tg',outTG,'poiseuille',outPoi,'outputRoot',root);
save(fullfile(root,'resamp_fluid_physics_validation_suite.mat'),'out','-v7.3');
end
function opts=parse_options(varargin)
opts=struct('outputRoot',fullfile('runs','resamp_fluid_physics_validation_suite'),'quick',false,'gamma',20,'NMin',14,'NMax',26,'kBT',0.01,'thermostatStrength',0.25,'rngSeed',12345,'visualEvery',0,'NxTG',64,'NyTG',64,'dtTG',1e-3,'tgSteps',3000,'tgSummaryEvery',50,'tgForceAmplitude',0.12,'NxPoi',64,'NyPoi',32,'dtPoi',0.005,'poiseuilleSteps',5000,'poiseuilleSampleEvery',50,'poiseuilleSummaryEvery',100,'bodyForceX',0.005);
if mod(numel(varargin),2)~=0, error('Options must be name/value pairs.'); end
for k=1:2:numel(varargin)
    key=lower(char(string(varargin{k}))); val=varargin{k+1};
    switch key
        case 'outputroot', opts.outputRoot=char(string(val));
        case 'quick', opts.quick=logical(val);
        case 'gamma', opts.gamma=val; opts.NMin=ceil(0.7*val); opts.NMax=ceil(1.3*val);
        case 'nmin', opts.NMin=val; case 'nmax', opts.NMax=val;
        case 'kbt', opts.kBT=val; case 'thermostatstrength', opts.thermostatStrength=val; case 'rngseed', opts.rngSeed=val; case 'visualevery', opts.visualEvery=val;
        case 'nxtg', opts.NxTG=val; case 'nytg', opts.NyTG=val; case 'dttg', opts.dtTG=val; case 'tgsteps', opts.tgSteps=val; case 'tgsummaryevery', opts.tgSummaryEvery=val; case {'tgforceamplitude','taylorgreenforcingamplitude'}, opts.tgForceAmplitude=val;
        case 'nxpoiseuille', opts.NxPoi=val; case 'nypoiseuille', opts.NyPoi=val; case 'dtpoiseuille', opts.dtPoi=val; case 'poiseuillesteps', opts.poiseuilleSteps=val; case 'poiseuillesampleevery', opts.poiseuilleSampleEvery=val; case 'poiseuillesummaryevery', opts.poiseuilleSummaryEvery=val; case 'bodyforcex', opts.bodyForceX=val;
        otherwise, error('Unknown option: %s', key);
    end
end
end
