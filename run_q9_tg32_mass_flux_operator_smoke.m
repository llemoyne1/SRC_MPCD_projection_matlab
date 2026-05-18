function results = run_q9_tg32_mass_flux_operator_smoke(varargin)
%RUN_Q9_TG32_MASS_FLUX_OPERATOR_SMOKE Compare Q9 mass-flux projectors on one state.
%
% This fast non-dynamic test initializes one TG state and applies:
%   1. historical periodic FFT Q9 projector,
%   2. FV/general-BC projector with elliptic low-pass target filter,
%   3. FV/general-BC projector with the exact periodic FFT low-k target mask.
%
% Case (2) is the pure no-FFT elliptic-filter path. It is not expected to
% match the FFT projector exactly because its low-pass filter is not the same
% idempotent spectral mask.
%
% Case (3) is a diagnostic hybrid. It uses the same low-k target subspace as
% the FFT projector, but solves the correction with the general FV elliptic
% operator. It separates target-filter differences from elliptic-inverse
% differences.
%
% Example:
%   results = run_q9_tg32_mass_flux_operator_smoke();

opts = parse_options(varargin{:});
params = make_base_params(opts);
rng(params.seed);
[state, initInfo] = projection_initialize_particles_taylor_green_forced(params); %#ok<ASGLU>
G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
    'periodicX', true, 'periodicY', true, 'minCount', 1);

pFFT = params;
pFFT.massFluxProjectionOperator = 'periodic_fft';
pFFT.massFluxTargetFilter = 'lowpass_fft';

projFFT = projection_project_mass_flux_periodic_fft(G.N, G.Ux, G.Uy, pFFT, ...
    'mode', params.massFluxProjectionMode, 'relaxationBeta', params.massFluxDensityRelaxationBeta);

cases = make_general_cases(params, opts);
rows = cell(numel(cases), 1);
projs = struct();
for ic = 1:numel(cases)
    c = cases(ic);
    projGEN = projection_project_mass_flux_general_bc(G.N, G.Ux, G.Uy, c.params, ...
        'mode', params.massFluxProjectionMode, 'relaxationBeta', params.massFluxDensityRelaxationBeta);
    D = compare_projection_fields(projFFT, projGEN);
    rows{ic} = table({c.label}, opts.Nx, opts.Ny, opts.gamma, ...
        projFFT.rmsDivMassBefore, projFFT.rmsDivMassAfter, projFFT.rmsDivMassResidual, ...
        projGEN.rmsDivMassBefore, projGEN.rmsDivMassAfter, projGEN.rmsDivMassResidual, ...
        projGEN.rmsDivMassBeforeFull, projGEN.rmsDivMassAfterFull, projGEN.rmsDivMassResidualFull, ...
        D.rmsDUxDiff, D.rmsDUyDiff, D.rmsDUDiff, D.relRmsDUDiff, D.maxAbsDUDiff, D.rmsTargetDiff, D.rmsResidualDiff, ...
        'VariableNames', {'case','Nx','Ny','gamma', ...
        'fftRmsDivBefore','fftRmsDivAfter','fftRmsResidual', ...
        'genRmsDivBefore','genRmsDivAfter','genRmsResidual', ...
        'genRmsDivBeforeFull','genRmsDivAfterFull','genRmsResidualFull', ...
        'rmsDUxDiff','rmsDUyDiff','rmsDUDiff','relRmsDUDiff','maxAbsDUDiff','rmsTargetDiff','rmsResidualDiff'});
    projs.(c.fieldName) = projGEN;
end
summary = vertcat(rows{:});

fprintf('\n=== Q9 TG32 one-state mass-flux operator smoke v7 ===\n');
disp(summary);

results = struct();
results.params = params;
results.grid = G;
results.projFFT = projFFT;
results.projGeneralBC = projs;
results.summary = summary;
end

function cases = make_general_cases(params, opts)
base = params;
base.massFluxProjectionOperator = 'general_bc';
base.massFluxEllipticAlphaMode = 'constant';
base.massFluxEllipticSolver = opts.generalEllipticSolver;
base.massFluxEllipticUseSolveCache = logical(opts.generalEllipticUseSolveCache);
base.massFluxEllipticFactorization = opts.generalEllipticFactorization;
base.massFluxBCLeft = 'periodic';
base.massFluxBCRight = 'periodic';
base.massFluxBCBottom = 'periodic';
base.massFluxBCTop = 'periodic';

pEll = base;
pEll.massFluxTargetFilter = 'elliptic_lowpass';

pFFTMask = base;
pFFTMask.massFluxTargetFilter = 'periodic_fft_lowk';

cases = repmat(struct('label','','fieldName','','params',struct()), 2, 1);
cases(1).label = 'fft_vs_general_bc_elliptic_lowpass';
cases(1).fieldName = 'general_bc_elliptic_lowpass';
cases(1).params = pEll;
cases(2).label = 'fft_vs_general_bc_periodic_fft_lowk';
cases(2).fieldName = 'general_bc_periodic_fft_lowk';
cases(2).params = pFFTMask;
end

function opts = parse_options(varargin)
opts = struct();
opts.Nx = 32;
opts.Ny = 32;
opts.gamma = 20;
opts.seed = 11;
opts.Lx = 1.0;
opts.Ly = 1.0;
opts.dt = 1.0e-3;
opts.kBT = 0.01;
opts.alphaDeg = 90;
opts.taylorGreenForceAmplitude = 0.12;
opts.taylorGreenInitialAmplitude = 0.10;
opts.massFluxDensityRelaxationBeta = 5.0e-4;
opts.massFluxLowKMaxIndex = 2;
opts.generalEllipticSolver = 'backslash';
opts.generalEllipticUseSolveCache = true;
opts.generalEllipticFactorization = 'auto';
if mod(numel(varargin),2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = char(string(varargin{k}));
    if ~isfield(opts, name)
        error('Unknown option: %s', name);
    end
    opts.(name) = varargin{k+1};
end
end

function params = make_base_params(opts)
params = struct();
params.method = 'q9';
params.Lx = opts.Lx;
params.Ly = opts.Ly;
params.Nx = opts.Nx;
params.Ny = opts.Ny;
params.gamma = opts.gamma;
params.seed = opts.seed;
params.dt = opts.dt;
params.kBT = opts.kBT;
params.alphaDeg = opts.alphaDeg;
params.initialPopulationMode = 'exact_per_cell';
params.initialVelocityZeroGlobalMean = true;
params.useRandomGridShift = true;
params.bodyForceX = 0;
params.bodyForceY = 0;
params.taylorGreenInitialAmplitude = opts.taylorGreenInitialAmplitude;
params.taylorGreenAmplitude = opts.taylorGreenInitialAmplitude;
params.taylorGreenThermalNoise = true;
params.taylorGreenForceEnable = true;
params.taylorGreenForceAmplitude = opts.taylorGreenForceAmplitude;
params.taylorGreenForceZeroMeanKick = true;
params.taylorGreenModeX = 1;
params.taylorGreenModeY = 1;
params.massFluxProjectionMode = 'relax_to_uniform_lowk';
params.massFluxProjectionStrength = 1.0;
params.massFluxDensityRelaxationBeta = opts.massFluxDensityRelaxationBeta;
params.massFluxApplyAfterVelocityProjection = true;
params.massFluxLowKMaxIndex = opts.massFluxLowKMaxIndex;
params.lowKMaxIndex = opts.massFluxLowKMaxIndex;
params.massFluxMinCellCount = 1.0;
end

function D = compare_projection_fields(A, B)
dUx = A.dUx - B.dUx;
dUy = A.dUy - B.dUy;
D = struct();
D.rmsDUxDiff = sqrt(mean(dUx(:).^2, 'omitnan'));
D.rmsDUyDiff = sqrt(mean(dUy(:).^2, 'omitnan'));
D.rmsDUDiff = sqrt(mean(dUx(:).^2 + dUy(:).^2, 'omitnan'));
D.rmsDURef = sqrt(mean(A.dUx(:).^2 + A.dUy(:).^2, 'omitnan'));
D.relRmsDUDiff = D.rmsDUDiff / max(D.rmsDURef, eps);
D.maxAbsDUDiff = max(max(abs(dUx(:))), max(abs(dUy(:))));
D.rmsTargetDiff = rms_field_diff(A.targetDivMass, B.targetDivMass);
D.rmsResidualDiff = rms_field_diff(A.divMassResidual, B.divMassResidual);
end

function v = rms_field_diff(A, B)
A = double(A(:));
B = double(B(:));
v = sqrt(mean((A - B).^2, 'omitnan'));
end
