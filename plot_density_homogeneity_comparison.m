function figs = plot_density_homogeneity_comparison(cmp, varargin)
%PLOT_DENSITY_HOMOGENEITY_COMPARISON Figures for density homogeneity comparison.

saveFigures = false;
outputDir = 'density_projection_output';
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "savefigures"
            saveFigures = logical(val);
        case "outputdir"
            outputDir = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end
if saveFigures && ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

mC = cmp.metricsClassic;
mP = cmp.metricsProjected;
t = mC.sampleTimes;
figs = gobjects(0);

fig = figure('Name', 'Density homogeneity metrics');
figs(end+1) = fig;
tiledlayout(4, 1);
nexttile;
plot(t, mC.stdN, '-', t, mP.stdN, '--', 'LineWidth', 1.2);
ylabel('std(N)'); grid on; legend('classic', 'projected', 'Location', 'best');
nexttile;
plot(t, mC.rmsRel, '-', t, mP.rmsRel, '--', 'LineWidth', 1.2);
ylabel('RMS(N/gamma-1)'); grid on; legend('classic', 'projected', 'Location', 'best');
nexttile;
semilogy(t, mC.lowKEnergy, '-', t, mP.lowKEnergy, '--', 'LineWidth', 1.2);
ylabel('low-k density energy'); grid on; legend('classic', 'projected', 'Location', 'best');
nexttile;
plot(t, mC.outBandFraction, '-', t, mP.outBandFraction, '--', 'LineWidth', 1.2);
ylabel('out-band fraction'); xlabel('t'); grid on; legend('classic', 'projected', 'Location', 'best');
sgtitle('Cell-population homogeneity');
save_if_requested(fig, outputDir, 'fig_density_homogeneity_timeseries.png', saveFigures);

fig = figure('Name', 'Density maps classic vs projected');
figs(end+1) = fig;
clim = max(abs([mC.timeAvgRel(:); mP.timeAvgRel(:); mC.finalRel(:); mP.finalRel(:)]));
if ~isfinite(clim) || clim <= 0
    clim = 1;
end
tiledlayout(2, 3);
nexttile; imagesc(mC.timeAvgRel.'); axis image; colorbar; caxis([-clim clim]); title('classic mean rel'); xlabel('x'); ylabel('y');
nexttile; imagesc(mP.timeAvgRel.'); axis image; colorbar; caxis([-clim clim]); title('projected mean rel'); xlabel('x'); ylabel('y');
nexttile; imagesc((mP.timeAvgRel - mC.timeAvgRel).'); axis image; colorbar; caxis([-clim clim]); title('projected - classic'); xlabel('x'); ylabel('y');
nexttile; imagesc(mC.finalRel.'); axis image; colorbar; caxis([-clim clim]); title('classic final rel'); xlabel('x'); ylabel('y');
nexttile; imagesc(mP.finalRel.'); axis image; colorbar; caxis([-clim clim]); title('projected final rel'); xlabel('x'); ylabel('y');
nexttile; imagesc((mP.finalRel - mC.finalRel).'); axis image; colorbar; caxis([-clim clim]); title('final diff'); xlabel('x'); ylabel('y');
sgtitle('Relative density maps N/gamma - 1');
save_if_requested(fig, outputDir, 'fig_density_maps_classic_projected.png', saveFigures);

fig = figure('Name', 'Density histograms');
figs(end+1) = fig;
edges = (0:max([mC.finalN(:); mP.finalN(:)])+1) - 0.5;
histogram(mC.finalN(:), edges, 'Normalization', 'probability');
hold on;
histogram(mP.finalN(:), edges, 'Normalization', 'probability');
xlabel('N per cell'); ylabel('probability'); grid on; legend('classic final', 'projected final');
title('Final cell-population histogram');
save_if_requested(fig, outputDir, 'fig_density_histograms.png', saveFigures);

fig = figure('Name', 'Density transport comparison');
figs(end+1) = fig;
HC = cmp.classic.diagHistory;
HP = cmp.projected.diagHistory;
plot(HC(:, 2), HC(:, 44), '-', HP(:, 2), HP(:, 44), '--', 'LineWidth', 1.2);
xlabel('t'); ylabel('continuous rho transport RMS'); grid on;
legend('classic run', 'projected run', 'Location', 'best');
title('Continuous density-transport proxy after each run correction');
save_if_requested(fig, outputDir, 'fig_density_transport_comparison.png', saveFigures);
end

function save_if_requested(fig, outputDir, name, saveFigures)
if ~saveFigures
    return;
end
try
    exportgraphics(fig, fullfile(outputDir, name), 'Resolution', 160);
catch
    saveas(fig, fullfile(outputDir, name));
end
end
