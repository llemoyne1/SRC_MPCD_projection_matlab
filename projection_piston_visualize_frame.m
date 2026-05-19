function projection_piston_visualize_frame(state, params, diagTable, step, sampleDiag, fig)
%PROJECTION_PISTON_VISUALIZE_FRAME Live dashboard for moving-piston validation.
%
% The visualization is intentionally lightweight but always shows the signals
% that detect piston/moving-wall failures early: yTop/compression, density
% homogeneity, pressure/kBT, and mean normal velocity.

if nargin < 6 || isempty(fig) || ~isgraphics(fig)
    fig = figure(get_param(params, 'visualFigureId', 620));
end
figure(fig);
if ~isappdata(fig, 'stopRequested')
    setappdata(fig, 'stopRequested', false);
end
clf(fig);
set(fig, 'Name', sprintf('Piston wallVP-v2 step %d', step));

G = sampleDiag.G;
Np = size(state.x, 1);
maxParticles = min(Np, round(get_param(params, 'visualMaxParticles', 5000)));
if maxParticles < Np
    ids = randperm(Np, maxParticles);
else
    ids = 1:Np;
end

tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% Particles and moving wall.
nexttile;
scatter(state.x(ids,1), state.x(ids,2), 4, state.v(ids,2), 'filled');
hold on;
yTop = sampleDiag.row.yTop;
plot([0 params.Lx], [0 0], 'k-', 'LineWidth', 1.0);
plot([0 params.Lx], [yTop yTop], 'k-', 'LineWidth', 1.5);
hold off;
axis([0 params.Lx 0 max(params.Ly0, params.Ly)]);
axis equal tight;
colorbar;
title(sprintf('particles v_y, step %d', step));
xlabel('x'); ylabel('y');

% Density map.
nexttile;
imagesc((G.N.' ./ max(params.gamma, eps)) - 1.0);
axis image;
colorbar;
title('N/\gamma - 1');
xlabel('x cell'); ylabel('y cell');

% Velocity / density profiles.
nexttile;
y = linspace(0, yTop, params.Ny).';
plot(sampleDiag.UyProfile, y, '-o', 'DisplayName','<U_y>'); hold on;
plot(sampleDiag.UxProfile, y, '-s', 'DisplayName','<U_x>');
plot((sampleDiag.NProfile/params.gamma - 1), y, '-.', 'DisplayName','N/\gamma-1');
hold off; grid on;
xlabel('profile'); ylabel('y');
legend('Location','best');
title('vertical profiles');

% Kinematics.
nexttile;
if ~isempty(diagTable)
    plot(diagTable.t, diagTable.yTop, '-', 'DisplayName','yTop'); hold on;
    plot(diagTable.t, diagTable.compression, '--', 'DisplayName','compression');
    hold off; grid on; legend('Location','best');
end
xlabel('t'); title('piston kinematics');

% Density / pressure / temperature / wall mechanics.
nexttile;
if ~isempty(diagTable)
    yyaxis left;
    plot(diagTable.t, diagTable.rhoPhysicalMean, '-', 'DisplayName','rho'); hold on;
    plot(diagTable.t, diagTable.PkinMean, '--', 'DisplayName','Pkin');
    if ismember('PkinIdealMean', diagTable.Properties.VariableNames)
        plot(diagTable.t, diagTable.PkinIdealMean, ':', 'DisplayName','rho*kBT');
    end
    if ismember('pressureTopWallTotal', diagTable.Properties.VariableNames) && any(isfinite(diagTable.pressureTopWallTotal))
        plot(diagTable.t, diagTable.pressureTopWallTotal, '-.', 'DisplayName','Ptop wall');
    end
    ylabel('rho, P');
    yyaxis right;
    plot(diagTable.t, diagTable.kBTCell, '-', 'DisplayName','kBT');
    ylabel('kBT');
    grid on; hold off;
end
xlabel('t'); title('rho / pressure / wall / kBT');

% Homogeneity / Q9 residual.
nexttile;
if ~isempty(diagTable)
    semilogy(diagTable.t, max(diagTable.lowKDensity, eps), '-', 'DisplayName','lowK density'); hold on;
    semilogy(diagTable.t, max(diagTable.relStdN, eps), '--', 'DisplayName','std(N)/mean(N)');
    if any(isfinite(diagTable.massFluxResidual))
        semilogy(diagTable.t, max(abs(diagTable.massFluxResidual), eps), '-.', 'DisplayName','mass flux residual');
    end
    if ismember('massFluxAfter', diagTable.Properties.VariableNames) && any(isfinite(diagTable.massFluxAfter))
        semilogy(diagTable.t, max(abs(diagTable.massFluxAfter), eps), ':', 'DisplayName','mass flux after');
    end
    hold off; grid on; legend('Location','best');
end
xlabel('t'); title('density homogeneity');

mf = sampleDiag.row.massFluxResidual;
if ~isfinite(mf), mf = NaN; end
Ptop = NaN;
if ismember('pressureTopWallTotal', sampleDiag.row.Properties.VariableNames)
    Ptop = sampleDiag.row.pressureTopWallTotal;
end
sgtitle(sprintf('%s | yTop=%.4g | comp=%.3g | rhoErr=%.1e | lowK=%.2e | Ptop=%.3g | mf=%.2e | kBT=%.4g', ...
    upper(char(string(get_param(params, 'method', 'classic')))), ...
    sampleDiag.row.yTop, sampleDiag.row.compression, sampleDiag.row.rhoPhysicalRelError, ...
    sampleDiag.row.lowKDensity, Ptop, mf, sampleDiag.row.kBTCell));

% Minimal stop button, recreated after clf.
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop run', ...
    'Units', 'normalized', 'Position', [0.88 0.94 0.10 0.045], ...
    'Callback', @(src,evt) setappdata(fig, 'stopRequested', true));
end

function value = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = defaultValue;
end
end
