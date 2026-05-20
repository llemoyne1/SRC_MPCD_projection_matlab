function results = run_q9_poiseuille_wall_virtual_particles_validation(varargin)
%RUN_Q9_POISEUILLE_WALL_VIRTUAL_PARTICLES_VALIDATION Fast wall-VP validation.
%
%   results = run_q9_poiseuille_wall_virtual_particles_validation(...)
%
% Default target: diagnose the wall no-slip problem in SRC classic first.
% Q6/Q9 can be included with 'methods', {'classic','q6','q9'}; they inherit
% the same classic collision substep.

p = inputParser;
addParameter(p, 'methods', {'classic'});
addParameter(p, 'Nx', 32);
addParameter(p, 'Ny', 48);
addParameter(p, 'gamma', 20);
addParameter(p, 'nSteps', 15000);
addParameter(p, 'sampleEvery', 100);
addParameter(p, 'progressEvery', 500);
addParameter(p, 'visualEnable', true);
addParameter(p, 'visualEvery', 500);
addParameter(p, 'bodyForceX', 0.005);
addParameter(p, 'wallModeY', 'bounceback');
addParameter(p, 'geometryMode', 'shifted_solid_fraction');
addParameter(p, 'densityFactor', 1.0);
addParameter(p, 'perCell', []);
addParameter(p, 'nLayers', 1);
addParameter(p, 'thermal', true);
addParameter(p, 'stochasticCount', true);
addParameter(p, 'forceRandomShiftY', true);
addParameter(p, 'initialPoiseuille', true);
addParameter(p, 'nuGuess', 0.021);
addParameter(p, 'fitWindowTime', 5);
addParameter(p, 'accelerationWindowTime', 5);
addParameter(p, 'nuReferenceNy', 64);
addParameter(p, 'nuScalingEnable', true);
addParameter(p, 'physicalLy', []);
addParameter(p, 'massFluxProjectionOperator', 'general_bc');
addParameter(p, 'massFluxTargetFilter', 'elliptic_lowpass');
addParameter(p, 'massFluxFinalVelocityProjectionCleanup', false);
parse(p, varargin{:});

methods = cellstr(p.Results.methods);
results = struct();
summaryRows = repmat(empty_row(), numel(methods), 1);

for im = 1:numel(methods)
    method = lower(strrep(methods{im}, '-', '_'));
    params = struct();
    params.method = method;
    params.Nx = p.Results.Nx;
    params.Ny = p.Results.Ny;
    params.gamma = p.Results.gamma;
    params.nSteps = p.Results.nSteps;
    params.sampleEvery = p.Results.sampleEvery;
    params.progressEvery = p.Results.progressEvery;
    params.visualEnable = p.Results.visualEnable;
    params.visualEvery = p.Results.visualEvery;
    params.bodyForceX = p.Results.bodyForceX;
    params.wallModeY = p.Results.wallModeY;
    params.fitWindowTime = p.Results.fitWindowTime;
    params.accelerationWindowTime = p.Results.accelerationWindowTime;
    params.visualAccelerationWindowTime = p.Results.accelerationWindowTime;
    params.poiseuilleNuReferenceNy = p.Results.nuReferenceNy;
    params.poiseuilleNuScalingEnable = p.Results.nuScalingEnable;
    if ~isempty(p.Results.physicalLy)
        params.poiseuillePhysicalLy = p.Results.physicalLy;
    end

    % Wall virtual particles v2.
    params.wallVirtualParticlesEnable = true;
    params.wallVirtualParticlesGeometryMode = p.Results.geometryMode;
    params.wallVirtualParticlesDensityFactor = p.Results.densityFactor;
    params.wallVirtualParticlesPerCell = p.Results.perCell;
    params.wallVirtualParticlesCells = p.Results.nLayers;
    params.wallVirtualParticlesThermal = p.Results.thermal;
    params.wallVirtualParticlesStochasticCount = p.Results.stochasticCount;
    params.wallVirtualParticlesForceRandomShiftY = p.Results.forceRandomShiftY;

    % Useful for Poiseuille calibration; not essential for wall diagnosis.
    params.initialPoiseuilleProfileEnable = p.Results.initialPoiseuille;
    params.initialPoiseuilleNuGuess = p.Results.nuGuess;
    params.initialPoiseuilleScale = 1.0;
    params.visualReferenceNu = p.Results.nuGuess;

    if strcmp(method, 'q9')
        params.massFluxProjectionOperator = char(string(p.Results.massFluxProjectionOperator));
        params.massFluxTargetFilter = char(string(p.Results.massFluxTargetFilter));
        params.massFluxFinalVelocityProjectionCleanup = logical(p.Results.massFluxFinalVelocityProjectionCleanup);
    end

    fprintf('\n=== Wall VP v2 validation: %s ===\n', upper(method));
    fprintf('grid=%dx%d, gamma=%g, force=%g, geom=%s, densityFactor=%g\n', ...
        params.Nx, params.Ny, params.gamma, params.bodyForceX, ...
        params.wallVirtualParticlesGeometryMode, params.wallVirtualParticlesDensityFactor);

    out = run_projection_poiseuille_medium64_demo(params);
    if isfield(out, 'accelerationRecent')
        acc = out.accelerationRecent;
    else
        acc = poiseuille_acceleration_diagnostics(out, p.Results.accelerationWindowTime);
    end

    row = empty_row();
    row.method = {char(method)};
    row.ok = true;
    row.actualStep = out.actualLastStep;
    row.elapsed = out.elapsedWallClock;
    row.slopeMeanUx = acc.slopeMeanUx;
    row.ratioAccel = acc.ratioAccel;
    row.finalMeanUx = out.diagTable.meanUx(end);
    row.finalCenterWall = out.diagTable.centerMinusWall(end);
    row.finalLowK = out.diagTable.lowKDensityEnergy(end);
    row.finalKBT = out.diagTable.kBTCell(end);
    row.nuEff = out.viscosity.nuEff;
    row.nuEffRaw = get_field(out.viscosity, 'nuEffRaw', out.viscosity.nuEff);
    row.nuEffScaledNyRef = get_field(out.viscosity, 'nuEffScaledNyRef', NaN);
    row.nuEffCellY = get_field(out.viscosity, 'nuEffCellY', NaN);
    row.nuScaleFactorNyRef = get_field(out.viscosity, 'nuScaleFactorNyRef', NaN);
    row.nuReferenceNy = get_field(out.viscosity, 'nuReferenceNy', NaN);
    row.dy = get_field(out.viscosity, 'dy', NaN);
    row.R2 = out.viscosity.R2;
    row.SNR = out.viscosity.signalToNoise;
    summaryRows(im) = row;

    results.(method) = out;
    fprintf('%s wallVP-v2: d<ux>/dt = %.6g, ratio = %.6g, final center-wall = %.6g, final lowK = %.3e, nuRaw = %.6g, nu@Ny%d = %.6g\n', ...
        upper(method), acc.slopeMeanUx, acc.ratioAccel, row.finalCenterWall, row.finalLowK, ...
        row.nuEffRaw, round(row.nuReferenceNy), row.nuEffScaledNyRef);
end

results.summary = struct2table(summaryRows);
disp(results.summary);
end

function row = empty_row()
row = struct('method', {{''}}, 'ok', false, 'actualStep', NaN, ...
    'elapsed', NaN, 'slopeMeanUx', NaN, 'ratioAccel', NaN, ...
    'finalMeanUx', NaN, 'finalCenterWall', NaN, 'finalLowK', NaN, ...
    'finalKBT', NaN, 'nuEff', NaN, 'nuEffRaw', NaN, 'nuEffScaledNyRef', NaN, ...
    'nuEffCellY', NaN, 'nuScaleFactorNyRef', NaN, 'nuReferenceNy', NaN, 'dy', NaN, ...
    'R2', NaN, 'SNR', NaN);
end

function val = get_field(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
end
end
