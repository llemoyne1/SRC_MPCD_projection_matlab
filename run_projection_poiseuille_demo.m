function out = run_projection_poiseuille_demo(params)
%RUN_PROJECTION_POISEUILLE_DEMO Poiseuille prototype with pressure projection.
%
%   out = run_projection_poiseuille_demo(params)
%
% Standalone Q3 demo: x-periodic channel, y walls, body force, classic SRD,
% algebraic pressure projection, and population diagnostics before/after the
% projection substep.

if nargin < 1 || isempty(params)
    params = struct();
end
params = set_default_params(params);

rng(params.seed);
Np = round(params.gamma * params.Nx * params.Ny);
state = struct();
state.x = [params.Lx * rand(Np, 1), params.Ly * rand(Np, 1)];
state.v = sqrt(params.kBT) * randn(Np, 2);
state.v(:, 1) = state.v(:, 1) - mean(state.v(:, 1));
state.v(:, 2) = state.v(:, 2) - mean(state.v(:, 2));

nSamples = floor(params.nSteps / params.sampleEvery) + 1;
sampleTimes = nan(nSamples, 1);
sampleSteps = nan(nSamples, 1);
UxProfiles = nan(params.Ny, nSamples);
UyProfiles = nan(params.Ny, nSamples);
NProfiles = nan(params.Ny, nSamples);
diagHistory = nan(nSamples, 15);
% columns: step,t,rmsDivBefore,rmsDivProjected,rmsDivParticle,kBTBefore,kBTAfter,
%          popStdBefore,popStdClassic,popStdProj,emptyBefore,emptyClassic,emptyProj,
%          popDeltaProjMax,dvRms

isamp = 0;
for it = 0:params.nSteps
    if mod(it, params.sampleEvery) == 0
        isamp = isamp + 1;
        G = projection_deposit_particles_to_grid(state.x, state.v, params, ...
            'periodicX', true, 'periodicY', false, 'minCount', 1);
        UxProfiles(:, isamp) = mean(G.Ux, 1).';
        UyProfiles(:, isamp) = mean(G.Uy, 1).';
        NProfiles(:, isamp) = mean(G.N, 1).';
        sampleTimes(isamp) = it * params.dt;
        sampleSteps(isamp) = it;
    end

    if it == params.nSteps
        break;
    end

    [state, stepDiag] = mpcd_step_projection_poiseuille(state, params);

    if mod(it + 1, params.sampleEvery) == 0
        % Fill the row that has just been created on next loop pass. The
        % diagnostic arrays are aligned to sample step after this step.
        row = floor((it + 1) / params.sampleEvery) + 1;
        if row <= nSamples
            diagHistory(row, :) = [it+1, (it+1)*params.dt, ...
                stepDiag.rmsDivBefore, stepDiag.rmsDivProjectedAfter, stepDiag.rmsDivParticleAfter, ...
                stepDiag.kBTBeforeProjection, stepDiag.kBTAfterProjection, ...
                stepDiag.populationBeforeStep.stdN, stepDiag.populationAfterClassic.stdN, stepDiag.populationAfterProjection.stdN, ...
                stepDiag.populationBeforeStep.nEmptyCells, stepDiag.populationAfterClassic.nEmptyCells, stepDiag.populationAfterProjection.nEmptyCells, ...
                stepDiag.populationProjectionDeltaMaxAbs, stepDiag.dvRms];
        end
    end
end

% Trim, in case sample schedule changed.
valid = isfinite(sampleTimes);
sampleTimes = sampleTimes(valid);
sampleSteps = sampleSteps(valid);
UxProfiles = UxProfiles(:, valid);
UyProfiles = UyProfiles(:, valid);
NProfiles = NProfiles(:, valid);
diagHistory = diagHistory(isfinite(diagHistory(:, 1)), :);

yCenters = ((0:params.Ny-1).' + 0.5) * params.Ly / params.Ny;

out = struct();
out.params = params;
out.state = state;
out.sampleTimes = sampleTimes;
out.sampleSteps = sampleSteps;
out.yCenters = yCenters;
out.UxProfiles = UxProfiles;
out.UyProfiles = UyProfiles;
out.NProfiles = NProfiles;
out.diagHistory = diagHistory;
out.diagColumns = {'step','t','rmsDivBefore','rmsDivProjected','rmsDivParticle', ...
    'kBTBefore','kBTAfter','popStdBefore','popStdClassic','popStdProjection', ...
    'emptyBefore','emptyClassic','emptyProjection','popDeltaProjectionMaxAbs','dvRms'};
out.summary = summarize_run(out);
out.viscosity = analyze_projection_poiseuille_viscosity(out, 'excludeWallCells', params.excludeWallCellsFit);

fprintf('\n=== run_projection_poiseuille_demo ===\n');
fprintf('Np                         : %d\n', Np);
fprintf('grid                       : %d x %d\n', params.Nx, params.Ny);
fprintf('wallModeY                  : %s\n', params.wallModeY);
fprintf('projectionStrength         : %.6g\n', params.projectionStrength);
fprintf('last rms div before        : %.12e\n', out.summary.lastRmsDivBefore);
fprintf('last rms div projected     : %.12e\n', out.summary.lastRmsDivProjected);
fprintf('last rms div particle grid : %.12e\n', out.summary.lastRmsDivParticle);
fprintf('last kBT after projection  : %.12e\n', out.summary.lastKBTAfter);
fprintf('last pop std classic/proj  : %.6g / %.6g\n', out.summary.lastPopStdClassic, out.summary.lastPopStdProjection);
fprintf('max pop delta projection   : %.6g\n', out.summary.maxPopDeltaProjection);
fprintf('nu_eff fit                 : %.12e\n', out.viscosity.nuEff);
fprintf('fit R2                     : %.6f\n', out.viscosity.R2);

if params.makeFigures
    make_figures(out);
end
end

function params = set_default_params(params)
params = set_default(params, 'Lx', 2.0);
params = set_default(params, 'Ly', 1.0);
params = set_default(params, 'Nx', 32);
params = set_default(params, 'Ny', 16);
params = set_default(params, 'gamma', 10);
params = set_default(params, 'dt', 1.0e-3);
params = set_default(params, 'kBT', 1.0);
params = set_default(params, 'alphaDeg', 90);
params = set_default(params, 'bodyForceX', 0.02);
params = set_default(params, 'bodyForceY', 0.0);
params = set_default(params, 'nSteps', 500);
params = set_default(params, 'sampleEvery', 10);
params = set_default(params, 'projectionEnable', true);
params = set_default(params, 'projectionStrength', 1.0);
params = set_default(params, 'projectionInterpolationMethod', 'nearest');
params = set_default(params, 'wallModeY', 'bounceback');
params = set_default(params, 'Ubottom', 0.0);
params = set_default(params, 'Utop', 0.0);
params = set_default(params, 'useRandomGridShiftX', true);
params = set_default(params, 'useRandomGridShiftY', false);
params = set_default(params, 'seed', 2);
params = set_default(params, 'makeFigures', true);
params = set_default(params, 'excludeWallCellsFit', 2);
end

function params = set_default(params, name, value)
if ~isfield(params, name) || isempty(params.(name))
    params.(name) = value;
end
end

function summary = summarize_run(out)
H = out.diagHistory;
summary = struct();
if isempty(H)
    summary.lastRmsDivBefore = NaN;
    summary.lastRmsDivProjected = NaN;
    summary.lastRmsDivParticle = NaN;
    summary.lastKBTAfter = NaN;
    summary.lastPopStdClassic = NaN;
    summary.lastPopStdProjection = NaN;
    summary.maxPopDeltaProjection = NaN;
    return;
end
summary.lastRmsDivBefore = H(end, 3);
summary.lastRmsDivProjected = H(end, 4);
summary.lastRmsDivParticle = H(end, 5);
summary.lastKBTAfter = H(end, 7);
summary.lastPopStdClassic = H(end, 9);
summary.lastPopStdProjection = H(end, 10);
summary.maxPopDeltaProjection = max(H(:, 14));
summary.meanPopStdClassic = mean(H(:, 9), 'omitnan');
summary.meanPopStdProjection = mean(H(:, 10), 'omitnan');
summary.meanDivReductionParticle = mean(H(:, 5) ./ max(H(:, 3), eps), 'omitnan');
end

function make_figures(out)
figure('Name', 'Projection Poiseuille profile');
hold on;
plot(out.yCenters, out.viscosity.UxMean, 'o-', 'DisplayName', 'time-avg Ux');
plot(out.yCenters, out.viscosity.UxFit, '-', 'LineWidth', 1.5, 'DisplayName', sprintf('fit nu=%.4g', out.viscosity.nuEff));
xlabel('y'); ylabel('Ux'); grid on; legend('Location', 'best');
title('Poiseuille profile with pressure projection');

if ~isempty(out.diagHistory)
    H = out.diagHistory;
    figure('Name', 'Projection Poiseuille diagnostics');
    tiledlayout(3, 1);
    nexttile;
    semilogy(H(:, 2), H(:, 3), '-', H(:, 2), H(:, 5), '-');
    ylabel('rms div'); grid on; legend('before', 'particle after');
    nexttile;
    plot(H(:, 2), H(:, 7), '-'); ylabel('kBT after'); grid on;
    nexttile;
    plot(H(:, 2), H(:, 9), '-', H(:, 2), H(:, 10), '--');
    xlabel('t'); ylabel('population std'); grid on; legend('after classic', 'after projection');
end
end
