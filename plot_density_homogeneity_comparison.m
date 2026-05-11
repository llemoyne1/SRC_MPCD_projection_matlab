function figs = plot_density_homogeneity_comparison(cmp, varargin)
%PLOT_DENSITY_HOMOGENEITY_COMPARISON Figures for density homogeneity comparison.
% Supports Q5 two-way cmp.classic/cmp.projected, Q7 cmp.virial, and
% Q7c cmp.repair three-way comparison structs.

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
hasVirial = isfield(cmp, 'metricsVirial') && ~isempty(cmp.metricsVirial);
hasRepair = isfield(cmp, 'metricsRepair') && ~isempty(cmp.metricsRepair);
hasThird = hasVirial || hasRepair;
if hasVirial
    mV = cmp.metricsVirial;
    thirdLabel = 'projected+virial';
    thirdField = 'virial';
elseif hasRepair
    mV = cmp.metricsRepair;
    thirdLabel = 'projected+repair';
    thirdField = 'repair';
end
t = mC.sampleTimes;
figs = gobjects(0);

fig = figure('Name', 'Density homogeneity metrics');
figs(end+1) = fig;
tiledlayout(4, 1);
nexttile;
if hasThird
    plot(t, mC.stdN, '-', t, mP.stdN, '--', t, mV.stdN, ':', 'LineWidth', 1.2);
    legend('classic', 'projected', thirdLabel, 'Location', 'best');
else
    plot(t, mC.stdN, '-', t, mP.stdN, '--', 'LineWidth', 1.2);
    legend('classic', 'projected', 'Location', 'best');
end
ylabel('std(N)'); grid on;
nexttile;
if hasThird
    plot(t, mC.rmsRel, '-', t, mP.rmsRel, '--', t, mV.rmsRel, ':', 'LineWidth', 1.2);
    legend('classic', 'projected', thirdLabel, 'Location', 'best');
else
    plot(t, mC.rmsRel, '-', t, mP.rmsRel, '--', 'LineWidth', 1.2);
    legend('classic', 'projected', 'Location', 'best');
end
ylabel('RMS(N/gamma-1)'); grid on;
nexttile;
if hasThird
    semilogy(t, mC.lowKEnergy, '-', t, mP.lowKEnergy, '--', t, mV.lowKEnergy, ':', 'LineWidth', 1.2);
    legend('classic', 'projected', thirdLabel, 'Location', 'best');
else
    semilogy(t, mC.lowKEnergy, '-', t, mP.lowKEnergy, '--', 'LineWidth', 1.2);
    legend('classic', 'projected', 'Location', 'best');
end
ylabel('low-k density energy'); grid on;
nexttile;
if hasThird
    plot(t, mC.outBandFraction, '-', t, mP.outBandFraction, '--', t, mV.outBandFraction, ':', 'LineWidth', 1.2);
    legend('classic', 'projected', thirdLabel, 'Location', 'best');
else
    plot(t, mC.outBandFraction, '-', t, mP.outBandFraction, '--', 'LineWidth', 1.2);
    legend('classic', 'projected', 'Location', 'best');
end
ylabel('out-band fraction'); xlabel('t'); grid on;
sgtitle('Cell-population homogeneity');
save_if_requested(fig, outputDir, 'fig_density_homogeneity_timeseries.png', saveFigures);

fig = figure('Name', 'Density maps');
figs(end+1) = fig;
if hasThird
    clim = max(abs([mC.timeAvgRel(:); mP.timeAvgRel(:); mV.timeAvgRel(:); ...
        mC.finalRel(:); mP.finalRel(:); mV.finalRel(:)]));
    if ~isfinite(clim) || clim <= 0
        clim = 1;
    end
    tiledlayout(2, 4);
    nexttile; imagesc(mC.timeAvgRel.'); axis image; colorbar; caxis([-clim clim]); title('classic mean rel'); xlabel('x'); ylabel('y');
    nexttile; imagesc(mP.timeAvgRel.'); axis image; colorbar; caxis([-clim clim]); title('projected mean rel'); xlabel('x'); ylabel('y');
    nexttile; imagesc(mV.timeAvgRel.'); axis image; colorbar; caxis([-clim clim]); title([thirdLabel ' mean rel']); xlabel('x'); ylabel('y');
    nexttile; imagesc((mV.timeAvgRel - mP.timeAvgRel).'); axis image; colorbar; caxis([-clim clim]); title([thirdLabel ' - projected']); xlabel('x'); ylabel('y');
    nexttile; imagesc(mC.finalRel.'); axis image; colorbar; caxis([-clim clim]); title('classic final rel'); xlabel('x'); ylabel('y');
    nexttile; imagesc(mP.finalRel.'); axis image; colorbar; caxis([-clim clim]); title('projected final rel'); xlabel('x'); ylabel('y');
    nexttile; imagesc(mV.finalRel.'); axis image; colorbar; caxis([-clim clim]); title([thirdLabel ' final rel']); xlabel('x'); ylabel('y');
    nexttile; imagesc((mV.finalRel - mP.finalRel).'); axis image; colorbar; caxis([-clim clim]); title('final third - P'); xlabel('x'); ylabel('y');
    sgtitle('Relative density maps N/gamma - 1');
else
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
end
save_if_requested(fig, outputDir, 'fig_density_maps_classic_projected.png', saveFigures);

fig = figure('Name', 'Density histograms');
figs(end+1) = fig;
if hasThird
    edges = (0:max([mC.finalN(:); mP.finalN(:); mV.finalN(:)])+1) - 0.5;
    histogram(mC.finalN(:), edges, 'Normalization', 'probability');
    hold on;
    histogram(mP.finalN(:), edges, 'Normalization', 'probability');
    histogram(mV.finalN(:), edges, 'Normalization', 'probability');
    legend('classic final', 'projected final', [thirdLabel ' final']);
else
    edges = (0:max([mC.finalN(:); mP.finalN(:)])+1) - 0.5;
    histogram(mC.finalN(:), edges, 'Normalization', 'probability');
    hold on;
    histogram(mP.finalN(:), edges, 'Normalization', 'probability');
    legend('classic final', 'projected final');
end
xlabel('N per cell'); ylabel('probability'); grid on;
title('Final cell-population histogram');
save_if_requested(fig, outputDir, 'fig_density_histograms.png', saveFigures);

fig = figure('Name', 'Density transport comparison');
figs(end+1) = fig;
HC = cmp.classic.diagHistory;
HP = cmp.projected.diagHistory;
if hasThird
    HV = cmp.(thirdField).diagHistory;
    plot(HC(:, 2), HC(:, 44), '-', HP(:, 2), HP(:, 44), '--', HV(:, 2), HV(:, 44), ':', 'LineWidth', 1.2);
    legend('classic run', 'projected run', [thirdLabel ' run'], 'Location', 'best');
else
    plot(HC(:, 2), HC(:, 44), '-', HP(:, 2), HP(:, 44), '--', 'LineWidth', 1.2);
    legend('classic run', 'projected run', 'Location', 'best');
end
xlabel('t'); ylabel('continuous rho transport RMS'); grid on;
title('Continuous density-transport proxy after each run correction');
save_if_requested(fig, outputDir, 'fig_density_transport_comparison.png', saveFigures);

if hasThird
    HV = cmp.(thirdField).diagHistory;
    if hasRepair
        fig = figure('Name', 'Density repair diagnostics');
        figs(end+1) = fig;
        tiledlayout(3, 1);
        nexttile;
        plot(HV(:, 2), HV(:, 65), '-', HV(:, 2), HV(:, 66), '--', 'LineWidth', 1.2);
        ylabel('repair |dx|'); grid on; legend('rms', 'max');
        nexttile;
        plot(HV(:, 2), HV(:, 67), '-', HV(:, 2), HV(:, 68), '--', 'LineWidth', 1.2);
        ylabel('std(N) repair'); grid on; legend('before', 'after');
        nexttile;
        semilogy(HV(:, 2), HV(:, 59), '-', HV(:, 2), HV(:, 74), '--', HV(:, 2), HV(:, 5), ':', 'LineWidth', 1.2);
        xlabel('t'); ylabel('rms div'); grid on; legend('after pressure', 'after repair', 'final');
        sgtitle('Local virial density-position repair diagnostics');
        save_if_requested(fig, outputDir, 'fig_density_repair_diagnostics.png', saveFigures);
    else
        fig = figure('Name', 'Virial kick diagnostics');
        figs(end+1) = fig;
        tiledlayout(3, 1);
        nexttile;
        plot(HV(:, 2), HV(:, 55), '-', HV(:, 2), HV(:, 56), '--', 'LineWidth', 1.2);
        ylabel('virial |dv|'); grid on; legend('rms', 'max');
        nexttile;
        plot(HV(:, 2), HV(:, 57), '-', HV(:, 2), HV(:, 58), '--', 'LineWidth', 1.2);
        ylabel('Pvir'); grid on; legend('rms', 'max abs');
        nexttile;
        semilogy(HV(:, 2), HV(:, 59), '-', HV(:, 2), HV(:, 60), '--', HV(:, 2), HV(:, 5), ':', 'LineWidth', 1.2);
        xlabel('t'); ylabel('rms div'); grid on; legend('after pressure', 'after virial', 'final');
        sgtitle('Local virial density-kick diagnostics');
        save_if_requested(fig, outputDir, 'fig_virial_kick_diagnostics.png', saveFigures);
    end
end
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
