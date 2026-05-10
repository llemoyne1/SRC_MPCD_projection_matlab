function out = run_projection_periodic_particles_demo(params)
%RUN_PROJECTION_PERIODIC_PARTICLES_DEMO Demo of particle-grid pressure projection.
%
%   out = run_projection_periodic_particles_demo()
%   out = run_projection_periodic_particles_demo(params)
%
% This is a small standalone run, independent from the legacy redistribution
% simulator. It is intended to validate the coupling:
%   particles -> grid -> pressure projection -> particles.

if nargin < 1 || isempty(params)
    params = struct();
end

params = set_default(params, 'Lx', 1.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 32);
params = set_default(params, 'gamma', 20);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 1.0);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.0);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'useRandomGridShift', true);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'nSteps', 100);
params = set_default(params, 'seed', 1);
params = set_default(params, 'makeFigures', true);

rng(params.seed);
Np = round(params.gamma * params.Nx * params.Ny);

state = struct();
state.x = [params.Lx * rand(Np, 1), params.Ly * rand(Np, 1)];
state.v = sqrt(params.kBT) * randn(Np, 2);
state.v = state.v - mean(state.v, 1);

% Add a smooth compressible velocity perturbation so the projection has a
% clear signal in early tests.
state.v(:, 1) = state.v(:, 1) + 0.25 * sin(2*pi*state.x(:, 1)/params.Lx);
state.v(:, 2) = state.v(:, 2) + 0.25 * sin(2*pi*state.x(:, 2)/params.Ly);

hist = initialize_history(params.nSteps);
lastDiag = [];
for step = 1:params.nSteps
    [state, diag] = mpcd_step_projection_periodic(state, params);
    lastDiag = diag;
    hist.step(step) = step;
    hist.kBTBeforeProjection(step) = diag.kBTBeforeProjection;
    hist.kBTAfterProjection(step) = diag.kBTAfterProjection;
    hist.rmsDivBefore(step) = diag.rmsDivBefore;
    hist.rmsDivProjectedAfter(step) = diag.rmsDivProjectedAfter;
    hist.rmsDivParticleAfter(step) = diag.rmsDivParticleAfter;
    hist.divReductionProjected(step) = diag.divReductionProjected;
    hist.divReductionParticle(step) = diag.divReductionParticle;
    hist.meanVx(step) = diag.meanVxAfterProjection;
    hist.meanVy(step) = diag.meanVyAfterProjection;
    hist.nEmptyCells(step) = diag.nEmptyCellsAfter;
    hist.dvRms(step) = diag.dvRms;
end

out = struct();
out.params = params;
out.state = state;
out.history = hist;
out.lastDiag = lastDiag;

fprintf('\n=== run_projection_periodic_particles_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('projectionStrength         : %.6g\n', params.projectionStrength);
fprintf('last rms div before        : %.12e\n', hist.rmsDivBefore(end));
fprintf('last rms div projected     : %.12e\n', hist.rmsDivProjectedAfter(end));
fprintf('last rms div particle grid : %.12e\n', hist.rmsDivParticleAfter(end));
fprintf('last kBT after projection  : %.12e\n', hist.kBTAfterProjection(end));

if params.makeFigures
    make_demo_figures(out);
end
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function hist = initialize_history(nSteps)
nanv = nan(nSteps, 1);
hist = struct();
hist.step = nanv;
hist.kBTBeforeProjection = nanv;
hist.kBTAfterProjection = nanv;
hist.rmsDivBefore = nanv;
hist.rmsDivProjectedAfter = nanv;
hist.rmsDivParticleAfter = nanv;
hist.divReductionProjected = nanv;
hist.divReductionParticle = nanv;
hist.meanVx = nanv;
hist.meanVy = nanv;
hist.nEmptyCells = nanv;
hist.dvRms = nanv;
end

function make_demo_figures(out)
h = out.history;

figure('Name', 'Projection divergence history');
semilogy(h.step, h.rmsDivBefore, '-', 'DisplayName', 'before projection');
hold on;
semilogy(h.step, h.rmsDivProjectedAfter, '-', 'DisplayName', 'grid projected');
semilogy(h.step, h.rmsDivParticleAfter, '-', 'DisplayName', 'after particle correction');
xlabel('step');
ylabel('RMS divergence');
title('Periodic particle-grid pressure projection');
grid on;
legend('Location', 'best');

figure('Name', 'Projection temperature history');
plot(h.step, h.kBTBeforeProjection, '-', 'DisplayName', 'before projection');
hold on;
plot(h.step, h.kBTAfterProjection, '-', 'DisplayName', 'after projection');
xlabel('step');
ylabel('estimated kBT');
title('Temperature diagnostic');
grid on;
legend('Location', 'best');
end
