%RUN_Q9_FORCED_TAYLOR_GREEN_VALIDATION
% Forced periodic Taylor-Green validation for classic/Q6/Q9.
%
% The goal is to impose a clean divergence-free vortical structure in a
% fully periodic bulk domain, then test whether Q6/Q9 suppress compressible
% density modes without destroying the forced vortex mode.

clearvars
%clc

%% User controls ----------------------------------------------------------
runMode = 'compare_classic_q6_q9';  % 'classic_only','q6_only','q9_only','compare_q6_q9','compare_classic_q9','compare_classic_q6_q9'
preset  = 'medium64';              % 'quick48','medium64','long96'

visualEnable = true;
visualSaveFrames = false;
meanFieldEnable = true;
meanFieldVisualEnable = true;
meanFieldStartStep = [];           % [] uses preset default; e.g. 3000
meanFieldSaveFrames = false;

% Forced TG parameters.  Start with a visible vortex and maintain it with a
% divergence-free body force of the same shape.
tgInitialAmplitude = 0.10;
tgForceAmplitude   = 0.08*1.5;
tgModeX = 1;
tgModeY = 1;

% Q9 reference parameters.
q9Beta = 5.0e-4;
q9LowKMaxIndex = 2;
q9CleanupStrength = 0.5;

% Leave [] to use preset values.
overrideNSteps = [];
overrideSampleEvery = [];
overrideVisualEvery = [];

%% Output setup ----------------------------------------------------------
stamp = datestr(now, 'yyyymmdd_HHMMSS');
outputDir = sprintf('q9_forced_taylor_green_%s_%s_%s', preset, runMode, stamp);
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
logFile = fullfile(outputDir, 'console_log.txt');
try, diary(logFile); catch ME, warning('Could not start diary: %s', ME.message); end
cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('=== Forced Taylor-Green validation ===\n');
fprintf('runMode : %s\n', runMode);
fprintf('preset  : %s\n', preset);
fprintf('output  : %s\n\n', outputDir);

%% Parameters ------------------------------------------------------------
base = make_base_tg_params(preset);
if ~isempty(overrideNSteps), base.nSteps = overrideNSteps; end
if ~isempty(overrideSampleEvery), base.sampleEvery = overrideSampleEvery; end
if ~isempty(overrideVisualEvery), base.visualEvery = overrideVisualEvery; end
if ~isempty(meanFieldStartStep), base.meanFieldStartStep = meanFieldStartStep; end

base.visualEnable = visualEnable;
base.visualSaveFrames = visualSaveFrames;
base.visualFigureId = 410;
base.visualFigureName = 'Forced Taylor-Green live';
base.visualFrameDir = fullfile(outputDir, 'frames');
base.visualFramePrefix = 'tg';
base.visualPause = 0.001;
base.visualMaxParticles = 8000;
base.visualDensityCLim = 0.6;
base.visualOmegaCLim = NaN;
base.visualQuiverStrideX = max(2, round(base.Nx/18));
base.visualQuiverStrideY = max(2, round(base.Ny/12));
base.visualQuiverScale = 1.2;

base.meanFieldEnable = meanFieldEnable;
base.meanFieldVisualEnable = meanFieldVisualEnable;
base.meanFieldFigureId = base.visualFigureId + 20;
base.meanFieldFigureName = 'Forced Taylor-Green running mean';
base.meanFieldSaveFrames = meanFieldSaveFrames;
base.meanFieldFrameDir = fullfile(outputDir, 'mean_frames');
base.meanFieldFramePrefix = 'tg_mean';

base.taylorGreenInitialAmplitude = tgInitialAmplitude;
base.taylorGreenAmplitude = tgInitialAmplitude;
base.taylorGreenForceAmplitude = tgForceAmplitude;
base.taylorGreenModeX = tgModeX;
base.taylorGreenModeY = tgModeY;

base.massFluxDensityRelaxationBeta = q9Beta;
base.massFluxLowKMaxIndex = q9LowKMaxIndex;
base.lowKMaxIndex = q9LowKMaxIndex;
base.massFluxFinalVelocityProjectionStrength = q9CleanupStrength;
base.projectionMomentumCorrectionEnable = true;
base.projectionMomentumCorrectionMode = 'particle_global_exact';

%% Methods ---------------------------------------------------------------
switch lower(runMode)
    case 'classic_only'
        methods = {'classic'};
    case 'q6_only'
        methods = {'q6'};
    case 'q9_only'
        methods = {'q9'};
    case 'compare_q6_q9'
        methods = {'q6','q9'};
    case 'compare_classic_q9'
        methods = {'classic','q9'};
    case 'compare_classic_q6_q9'
        methods = {'classic','q6','q9'};
    otherwise
        error('Unknown runMode: %s', runMode);
end

outs = struct();
for im = 1:numel(methods)
    method = methods{im};
    p = configure_method(base, method);
    p.visualFigureId = base.visualFigureId + im - 1;
    p.visualFigureName = sprintf('Forced Taylor-Green %s', upper(method));
    p.visualTitleSuffix = upper(method);
    p.visualFrameDir = fullfile(outputDir, 'frames', method);
    p.visualFramePrefix = sprintf('tg_%s', method);
    p.meanFieldFigureId = base.meanFieldFigureId + im - 1;
    p.meanFieldFigureName = sprintf('Forced Taylor-Green running mean %s', upper(method));
    p.meanFieldFrameDir = fullfile(outputDir, 'mean_frames', method);
    p.meanFieldFramePrefix = sprintf('tg_mean_%s', method);

    fprintf('\n--- Running forced TG %s ---\n', upper(method));
    out = run_projection_forced_taylor_green_demo(p);
    outs.(method) = out;
    tsFile = fullfile(outputDir, sprintf('forced_tg_%s_timeseries.csv', method));
    write_tg_timeseries_csv(tsFile, out);
    fprintf('Timeseries written to: %s\n', tsFile);
end

summaryTable = build_tg_summary_table(outs, methods);
disp(summaryTable);
csvFile = fullfile(outputDir, 'forced_tg_summary.csv');
txtFile = fullfile(outputDir, 'forced_tg_summary.txt');
matFile = fullfile(outputDir, 'forced_tg_validation.mat');
figFile = fullfile(outputDir, 'forced_tg_timeseries.png');
writetable(summaryTable, csvFile);
make_tg_summary_figure(outs, methods, figFile);
fitFigFile = fullfile(outputDir, 'forced_tg_viscosity_fit.png');
make_tg_viscosity_fit_figure(outs, methods, fitFigFile);
write_tg_summary_txt(txtFile, base, methods, summaryTable, matFile, csvFile, figFile, fitFigFile);
save(matFile, 'base', 'methods', 'outs', 'summaryTable', '-v7.3');

fprintf('\nSummary written to:\n%s\n', txtFile);
fprintf('CSV written to:\n%s\n', csvFile);
fprintf('MAT written to:\n%s\n', matFile);
fprintf('Figure written to:\n%s\n', figFile);
fprintf('Viscosity fit figure written to:\n%s\n', fitFigFile);

%% Local functions -------------------------------------------------------
function p = make_base_tg_params(preset)
p = struct();
p.Lx = 1.0; p.Ly = 1.0; p.gamma = 20; p.seed = 11;
p.dt = 1.0e-3; p.kBT = 0.01; p.alphaDeg = 90;
p.initialPopulationMode = 'exact_per_cell';
p.initialVelocityZeroGlobalMean = true;
p.taylorGreenThermalNoise = true;
p.taylorGreenForceEnable = true;
p.taylorGreenForceZeroMeanKick = true;
p.bodyForceX = 0; p.bodyForceY = 0;
p.useRandomGridShift = true;
p.projectionInterpolationMethod = 'nearest';
% Use the same post-step thermostat for classic, Q6 and Q9.
% In Q6/Q9 it is applied after projection/momentum correction; in classic
% it is applied after the SRD collision.
p.thermostatAfterStep = true;
p.thermostatAfterProjection = true;
p.computeFullDiagnosticsEveryStep = false;
p.abortOnUnstable = true;
p.maxStableKBTCell = 0.2;
p.maxStableHighKFraction = 0.9;
p.maxStableDensityRelRMS = 1.2;
p.maxWallClockSeconds = Inf;
p.taylorGreenViscosityFitFraction = 0.60;
switch lower(preset)
    case 'quick48'
        p.Nx = 48; p.Ny = 48; p.nSteps = 10000; p.sampleEvery = 50; p.progressEvery = 500; p.visualEvery = 250;
    case 'medium64'
        p.Nx = 64; p.Ny = 64; p.nSteps = 5000; p.sampleEvery = 100; p.progressEvery = 1000; p.visualEvery = 250;
    case 'long96'
        p.Nx = 96; p.Ny = 96; p.nSteps = 50000; p.sampleEvery = 200; p.progressEvery = 1000; p.visualEvery = 1000;
    otherwise
        error('Unknown preset: %s', preset);
end
p.meanFieldStartStep = round(0.2*p.nSteps);
p.massFluxProjectionMode = 'relax_to_uniform_lowk';
p.massFluxProjectionStrength = 1.0;
p.massFluxDensityRelaxationBeta = 5e-4;
p.massFluxApplyAfterVelocityProjection = true;
p.massFluxLowKMaxIndex = 2;
p.massFluxFinalVelocityProjectionCleanup = true;
p.massFluxFinalVelocityProjectionStrength = 0.5;
p.projectionStrength = 1.0;
p.projectionEnable = true;
end

function p = configure_method(base, method)
p = base; p.method = lower(method);
switch lower(method)
    case 'classic'
        p.projectionEnable = false; p.projectionStrength = 0; p.massFluxProjectionMode = 'off'; p.massFluxProjectionStrength = 0; p.massFluxDensityRelaxationBeta = 0; p.massFluxApplyAfterVelocityProjection = false; p.massFluxFinalVelocityProjectionCleanup = false; p.massFluxFinalVelocityProjectionStrength = 0; p.projectionMomentumCorrectionEnable = false;
    case 'q6'
        p.projectionEnable = true; p.projectionStrength = 1; p.massFluxProjectionMode = 'off'; p.massFluxProjectionStrength = 0; p.massFluxDensityRelaxationBeta = 0; p.massFluxApplyAfterVelocityProjection = false; p.massFluxFinalVelocityProjectionCleanup = false; p.massFluxFinalVelocityProjectionStrength = 0;
    case 'q9'
        p.projectionEnable = true; p.projectionStrength = 1; p.massFluxProjectionMode = 'relax_to_uniform_lowk'; p.massFluxProjectionStrength = 1; p.massFluxApplyAfterVelocityProjection = true; p.massFluxFinalVelocityProjectionCleanup = true;
    otherwise
        error('Unknown method: %s', method);
end
end

function write_tg_timeseries_csv(filename, out)
S = out.series;
T = table(S.step, S.t, S.tgAmplitude, S.tgCoherence, S.tgModeEnergy, S.tgResidualEnergy, S.enstrophy, S.omegaRms, S.highKVelocityFraction, S.lowKVelocityFraction, S.lowKDensity, S.densityRelRMS, S.rmsDivGrid, S.divParticleAfter, S.massFluxResidual, S.kBTCell, S.NStd, S.NOutBand, S.dAmplitudeDt, S.nuEffForced, S.momentumRawDV, S.momentumResidualDV, ...
    'VariableNames', {'step','t','tgAmplitude','tgCoherence','tgModeEnergy','tgResidualEnergy','enstrophy','omegaRms','highKVelocityFraction','lowKVelocityFraction','lowKDensity','densityRelRMS','rmsDivGrid','divParticleAfter','massFluxResidual','kBTCell','NStd','NOutBand','dAmplitudeDt','nuEffForced','momentumRawDV','momentumResidualDV'});
writetable(T, filename);
end

function T = build_tg_summary_table(outs, methods)
rows = cell(numel(methods),1); vals = nan(numel(methods),29); pref = cell(numel(methods),1);
for i=1:numel(methods)
    rows{i}=upper(methods{i}); s=outs.(methods{i}).summary;
    pref{i}=char(string(getfield_with_default(s,'viscosityFitMethod',''))); %#ok<GFLD>
    vals(i,:) = [ ...
        getfield_with_default(s,'actualLastStep',NaN), ...
        double(getfield_with_default(s,'stoppedEarly',NaN)), ...
        getfield_with_default(s,'meanAmplitude',NaN), ...
        getfield_with_default(s,'finalAmplitude',NaN), ...
        getfield_with_default(s,'meanCoherence',NaN), ...
        getfield_with_default(s,'finalCoherence',NaN), ...
        getfield_with_default(s,'meanModeEnergy',NaN), ...
        getfield_with_default(s,'meanEnstrophy',NaN), ...
        getfield_with_default(s,'finalEnstrophy',NaN), ...
        getfield_with_default(s,'meanHighKFraction',NaN), ...
        getfield_with_default(s,'finalHighKFraction',NaN), ...
        getfield_with_default(s,'meanLowKDensity',NaN), ...
        getfield_with_default(s,'finalLowKDensity',NaN), ...
        getfield_with_default(s,'meanDensityRelRMS',NaN), ...
        getfield_with_default(s,'meanKBTCell',NaN), ...
        getfield_with_default(s,'finalKBTCell',NaN), ...
        getfield_with_default(s,'meanNuEffForced',NaN), ...
        getfield_with_default(s,'finalNuEffForced',NaN), ...
        getfield_with_default(s,'nuEffFitPreferred',NaN), ...
        getfield_with_default(s,'nuEffFitDerivativeKnownF',NaN), ...
        getfield_with_default(s,'nuEffFitDerivativeFreeF',NaN), ...
        getfield_with_default(s,'tgForceFit',NaN), ...
        getfield_with_default(s,'nuEffFitExpKnownF',NaN), ...
        getfield_with_default(s,'nuEffFitPlateau',NaN), ...
        getfield_with_default(s,'viscosityFitR2',NaN), ...
        getfield_with_default(s,'viscosityFitT0',NaN), ...
        getfield_with_default(s,'viscosityFitT1',NaN), ...
        getfield_with_default(s,'meanMomentumRawDV',NaN), ...
        getfield_with_default(s,'meanMomentumResidualDV',NaN)];
end
T = array2table(vals, 'VariableNames', {'actualLastStep','stoppedEarly','meanAmplitude','finalAmplitude','meanCoherence','finalCoherence','meanModeEnergy','meanEnstrophy','finalEnstrophy','meanHighKFraction','finalHighKFraction','meanLowKDensity','finalLowKDensity','meanDensityRelRMS','meanKBTCell','finalKBTCell','meanNuEffForced','finalNuEffForced','nuEffFitPreferred','nuEffFitDerivativeKnownF','nuEffFitDerivativeFreeF','tgForceFit','nuEffFitExpKnownF','nuEffFitPlateau','viscosityFitR2','viscosityFitT0','viscosityFitT1','meanMomentumRawDV','meanMomentumResidualDV'});
T.viscosityFitMethod = pref;
T.method = rows; T = movevars(T, 'method', 'Before', 1);
end

function make_tg_summary_figure(outs, methods, figFile)
figure(1); clf; set(gcf,'Name','Forced Taylor-Green validation diagnostics','Color','w');
plots = {'tgAmplitude','tgCoherence','enstrophy','highKVelocityFraction','lowKDensity','kBTCell'};
titles = {'TG amplitude','mode coherence','enstrophy','high-k velocity fraction','low-k density','kBT cell'};
for ip=1:numel(plots)
    subplot(2,3,ip); hold on;
    for im=1:numel(methods)
        S=outs.(methods{im}).series; plot(S.t, S.(plots{ip}), 'DisplayName', upper(methods{im}));
    end
    hold off; grid on; title(titles{ip}); xlabel('t'); legend('Location','best');
end
exportgraphics(gcf, figFile);
end

function make_tg_viscosity_fit_figure(outs, methods, figFile)
figure(2); clf; set(gcf,'Name','Forced Taylor-Green viscosity fit','Color','w');
subplot(1,2,1); hold on;
for im=1:numel(methods)
    method=methods{im}; S=outs.(method).series;
    plot(S.t, S.tgAmplitude, 'DisplayName', upper(method));
    fit=outs.(method).summary.viscosityFit;
    if isstruct(fit) && isfield(fit,'lambdaExpKnownF') && isfinite(fit.lambdaExpKnownF)
        F=fit.forceKnown; lam=fit.lambdaExpKnownF; C=fit.CExpKnownF;
        tau=S.t-fit.tStart; mask=S.t>=fit.tStart & S.t<=fit.tEnd;
        Apred=F/max(lam,eps)+C*exp(-lam*(S.t-fit.tStart));
        plot(S.t(mask), Apred(mask), '--', 'HandleVisibility','off');
    end
end
hold off; grid on; xlabel('t'); ylabel('A(t)'); title('TG amplitude and fitted exponential'); legend('Location','best');
subplot(1,2,2);
nu = nan(numel(methods),3); labels=cell(numel(methods),1);
for im=1:numel(methods)
    method=methods{im}; s=outs.(method).summary; labels{im}=upper(method);
    nu(im,:)=[s.nuEffFitPreferred, s.nuEffFitDerivativeKnownF, s.nuEffFitPlateau];
end
bar(nu); grid on; set(gca,'XTickLabel',labels); ylabel('\nu_{eff}'); title('Effective viscosity estimates');
legend({'preferred','dA/dt fit','plateau'}, 'Location','best');
exportgraphics(gcf, figFile);
end

function write_tg_summary_txt(filename, base, methods, tableSummary, matFile, csvFile, figFile, fitFigFile)
fid=fopen(filename,'w'); if fid<0, warning('Could not write %s', filename); return; end
fprintf(fid,'=== Forced Taylor-Green validation ===\n');
fprintf(fid,'methods                         : %s\n', strjoin(cellfun(@upper, methods, 'UniformOutput', false), ', '));
fprintf(fid,'steps/sampleEvery/visualEvery   : %d / %d / %d\n', base.nSteps, base.sampleEvery, base.visualEvery);
fprintf(fid,'grid, gamma                     : %d x %d, %g\n', base.Nx, base.Ny, base.gamma);
fprintf(fid,'dt, kBT                         : %.12g, %.12g\n', base.dt, base.kBT);
fprintf(fid,'TG initial/force amplitude       : %.12g / %.12g\n', base.taylorGreenInitialAmplitude, base.taylorGreenForceAmplitude);
fprintf(fid,'TG mode                         : (%d,%d)\n', base.taylorGreenModeX, base.taylorGreenModeY);
fprintf(fid,'Q9 beta/lowK/cleanup            : %.12g / %d / %.12g\n\n', base.massFluxDensityRelaxationBeta, base.massFluxLowKMaxIndex, base.massFluxFinalVelocityProjectionStrength);
fprintf(fid,'--- Summary table ---\n');
for i=1:height(tableSummary)
    fprintf(fid,'\nmethod: %s\n', tableSummary.method{i});
    fprintf(fid,'actual last step / stopped early : %.0f / %.0f\n', tableSummary.actualLastStep(i), tableSummary.stoppedEarly(i));
    fprintf(fid,'mean/final amplitude             : %.12g / %.12g\n', tableSummary.meanAmplitude(i), tableSummary.finalAmplitude(i));
    fprintf(fid,'mean/final coherence             : %.12g / %.12g\n', tableSummary.meanCoherence(i), tableSummary.finalCoherence(i));
    fprintf(fid,'mean/final enstrophy             : %.12g / %.12g\n', tableSummary.meanEnstrophy(i), tableSummary.finalEnstrophy(i));
    fprintf(fid,'mean/final high-k fraction        : %.12g / %.12g\n', tableSummary.meanHighKFraction(i), tableSummary.finalHighKFraction(i));
    fprintf(fid,'mean/final low-k density          : %.12g / %.12g\n', tableSummary.meanLowKDensity(i), tableSummary.finalLowKDensity(i));
    fprintf(fid,'mean density rel RMS              : %.12g\n', tableSummary.meanDensityRelRMS(i));
    fprintf(fid,'mean/final kBT cell               : %.12g / %.12g\n', tableSummary.meanKBTCell(i), tableSummary.finalKBTCell(i));
    fprintf(fid,'mean nuEff forced estimate        : %.12g\n', tableSummary.meanNuEffForced(i));
    fprintf(fid,'nuEff fit preferred/method        : %.12g / %s\n', tableSummary.nuEffFitPreferred(i), tableSummary.viscosityFitMethod{i});
    fprintf(fid,'nuEff derivative known/free F     : %.12g / %.12g\n', tableSummary.nuEffFitDerivativeKnownF(i), tableSummary.nuEffFitDerivativeFreeF(i));
    fprintf(fid,'F fitted / F configured           : %.12g / %.12g\n', tableSummary.tgForceFit(i), base.taylorGreenForceAmplitude);
    fprintf(fid,'nuEff exp/plateau/R2              : %.12g / %.12g / %.12g\n', tableSummary.nuEffFitExpKnownF(i), tableSummary.nuEffFitPlateau(i), tableSummary.viscosityFitR2(i));
    fprintf(fid,'fit window                        : %.12g -- %.12g\n', tableSummary.viscosityFitT0(i), tableSummary.viscosityFitT1(i));
    fprintf(fid,'mean raw momentum kick per part.  : %.12g\n', table_value_default(tableSummary,'meanMomentumRawDV',i,NaN));
    fprintf(fid,'mean residual momentum kick/part. : %.12g\n', table_value_default(tableSummary,'meanMomentumResidualDV',i,NaN));
end
fprintf(fid,'\n--- Files ---\nMAT: %s\nCSV: %s\nFigure: %s\nViscosity fit figure: %s\n', matFile, csvFile, figFile, fitFigFile);
fclose(fid);
end



function v = table_value_default(T, name, row, defaultValue)
if istable(T) && any(strcmp(T.Properties.VariableNames, name)) && row <= height(T)
    v = T.(name)(row);
else
    v = defaultValue;
end
end

function v = getfield_with_default(s, name, defaultValue)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
