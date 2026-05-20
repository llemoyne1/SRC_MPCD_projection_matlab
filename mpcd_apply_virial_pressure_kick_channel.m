function [stateOut, info] = mpcd_apply_virial_pressure_kick_channel(stateIn, params)
%MPCD_APPLY_VIRIAL_PRESSURE_KICK_CHANNEL Virial EOS diagnostics and optional kick.
%
%   [stateOut, info] = mpcd_apply_virial_pressure_kick_channel(stateIn, params)
%
% This is a clean/Q6-Q9 transplant of the historical liquid-closure virial
% idea, without redistribution or velocity repair.  It is intended to be
% called after the validated Q6/Q9 block in the moving-piston channel.
%
% The diagnostic EOS pressure is
%      PtotEOS = Pkin + PvirEOS,
%      PvirEOS = Kvirial * (rho - rhoEOSRef),
% where rhoEOSRef defaults to the initial physical density gamma/(dx0*dy0).
%
% The active kick is optional and uses
%      dv = - betaVirial * dt/rhoKick * grad(Pdrive),
% with Pdrive = Pkin + PvirDrive.  Since the reference density is spatially
% uniform in the current piston geometry, choosing initial or current target
% only changes the pressure offset, not the kick gradient.  The default kick
% target is therefore the current uniform density, as in the historical
% closure logic.

if nargin < 2
    error('Usage: [stateOut, info] = mpcd_apply_virial_pressure_kick_channel(stateIn, params)');
end
validate_state(stateIn);

stateOut = stateIn;
enableDiag = logical(getp(params, 'virialDiagnosticsEnable', false));
enableKick = logical(getp(params, 'virialKickEnable', false));
Kvirial = getp(params, 'Kvirial', getp(params, 'virialK', 0.0));
betaVirial = getp(params, 'virialBeta', getp(params, 'betaEOS', 0.0));

info = empty_info(params, Kvirial, betaVirial, enableDiag, enableKick);
if ~enableDiag && ~enableKick
    return;
end

periodicX = logical(getp(params, 'virialPeriodicX', true));
periodicY = logical(getp(params, 'virialPeriodicY', false));
interpMethod = char(string(getp(params, 'virialInterpolationMethod', getp(params, 'projectionInterpolationMethod', 'nearest'))));
minCount = getp(params, 'virialMinCellCount', getp(params, 'massFluxMinCellCount', 1));

G = projection_deposit_particles_to_grid(stateIn.x, stateIn.v, params, ...
    'periodicX', periodicX, 'periodicY', periodicY, 'minCount', minCount);
Tmap = virial_temperature_map(stateIn.x, stateIn.v, params, G, periodicX, periodicY);
Pkin = G.rho .* max(Tmap, 0.0);
Pkin(~isfinite(Pkin)) = 0.0;

rhoMean = size(stateIn.x, 1) / max(params.Lx * params.Ly, eps);
rhoEOSRef = resolve_rho_eos_ref(params, stateIn, rhoMean);
rhoUniformNow = resolve_rho_uniform_now(params, stateIn, rhoMean);

PvirEOS = Kvirial .* (G.rho - rhoEOSRef);
PtotEOS = Pkin + PvirEOS;
PvirEOS(~isfinite(PvirEOS)) = 0.0;
PtotEOS(~isfinite(PtotEOS)) = 0.0;

driveTargetMode = lower(strrep(char(string(getp(params, 'virialDriveTargetMode', 'current_uniform'))), '-', '_'));
switch driveTargetMode
    case {'current_uniform','uniform_now','rho_uniform_now','historical'}
        rhoDriveRef = rhoUniformNow;
    case {'eos_ref','initial','initial_density','rho_eos_ref'}
        rhoDriveRef = rhoEOSRef;
    case {'zero','none'}
        rhoDriveRef = 0.0;
    otherwise
        error('Unknown virialDriveTargetMode: %s', driveTargetMode);
end
PvirDrive = Kvirial .* (G.rho - rhoDriveRef);
Pdrive = Pkin + PvirDrive;
PvirDrive(~isfinite(PvirDrive)) = 0.0;
Pdrive(~isfinite(Pdrive)) = 0.0;

[dPdx, dPdy] = virial_gradient_boundary_aware(Pdrive, params, periodicX, periodicY);
rhoKick = resolve_rho_kick(params, G.rho, rhoUniformNow);
maskKick = isfinite(rhoKick) & rhoKick > 0 & isfinite(dPdx) & isfinite(dPdy);

duXRaw = zeros(size(Pdrive));
duYRaw = zeros(size(Pdrive));
if any(maskKick(:))
    duXRaw(maskKick) = -(params.dt ./ rhoKick(maskKick)) .* dPdx(maskKick);
    duYRaw(maskKick) = -(params.dt ./ rhoKick(maskKick)) .* dPdy(maskKick);
end

duXApplied = betaVirial .* duXRaw;
duYApplied = betaVirial .* duYRaw;
[duXApplied, duYApplied, limiterInfo] = apply_du_limiter(duXApplied, duYApplied, params);

info.enabled = enableDiag || enableKick;
info.diagnosticsEnabled = enableDiag;
info.kickEnabled = enableKick;
info.Kvirial = Kvirial;
info.betaVirial = betaVirial;
info.rhoEOSRef = rhoEOSRef;
info.rhoUniformNow = rhoUniformNow;
info.rhoDriveRef = rhoDriveRef;
info.rhoMean = rhoMean;
info.rhoDefectRms = sqrt(mean((G.rho(:) - rhoUniformNow).^2, 'omitnan'));
info.rhoDefectRelRms = info.rhoDefectRms / max(abs(rhoUniformNow), eps);
info.rhoMin = min(G.rho(:));
info.rhoMax = max(G.rho(:));
info.PkinMean = mean(Pkin(:), 'omitnan');
info.PvirMean = mean(PvirEOS(:), 'omitnan');
info.PtotMean = mean(PtotEOS(:), 'omitnan');
info.PdriveMean = mean(Pdrive(:), 'omitnan');
info.PkinMin = min(Pkin(:));
info.PkinMax = max(Pkin(:));
info.PvirMin = min(PvirEOS(:));
info.PvirMax = max(PvirEOS(:));
info.PtotMin = min(PtotEOS(:));
info.PtotMax = max(PtotEOS(:));
info.gradPdriveRms = sqrt(mean(dPdx(:).^2 + dPdy(:).^2, 'omitnan'));
info.gradPdriveMaxAbs = max(sqrt(dPdx(:).^2 + dPdy(:).^2));
info.duVirialRawRms = sqrt(mean(duXRaw(:).^2 + duYRaw(:).^2, 'omitnan'));
info.duVirialAppliedRms = sqrt(mean(duXApplied(:).^2 + duYApplied(:).^2, 'omitnan'));
info.duVirialAppliedMaxAbs = max(sqrt(duXApplied(:).^2 + duYApplied(:).^2));
info.duVirialOverThermalRms = info.duVirialAppliedRms / max(sqrt(getp(params, 'kBT', 1.0)), eps);
info.limiter = limiterInfo;
info.limitedCellFraction = limiterInfo.limitedCellFraction;
info.limitedCellCount = getf(limiterInfo, 'nLimitedCells', NaN);
info.duLimiterMax = getf(limiterInfo, 'duMax', NaN);
info.duRmsBeforeLimiter = getf(limiterInfo, 'duRmsBeforeLimiter', NaN);
info.duRmsAfterLimiter = getf(limiterInfo, 'duRmsAfterLimiter', NaN);
info.duMaxFractionThermal = getp(params, 'virialMaxDuFractionThermal', NaN);
info.Pkin = Pkin;
info.PvirEOS = PvirEOS;
info.PtotEOS = PtotEOS;
info.PvirDrive = PvirDrive;
info.Pdrive = Pdrive;
info.dPdx = dPdx;
info.dPdy = dPdy;
info.duVirialX = duXApplied;
info.duVirialY = duYApplied;

if enableKick && betaVirial ~= 0 && Kvirial ~= 0
    dvRaw = projection_interpolate_grid_delta_to_particles(stateIn.x, duXApplied, duYApplied, params, ...
        'periodicX', periodicX, 'periodicY', periodicY, 'method', interpMethod);
    [stateOut.v, momentumInfo] = projection_apply_global_momentum_correction(stateIn.v, dvRaw, params, 'stage', 'virial');
    info.momentumCorrection = momentumInfo;
    info.rawMomentumKick = momentumInfo.rawDeltaP;
    info.rawMomentumKickNorm = momentumInfo.rawDeltaPNorm;
    info.residualMomentumKick = momentumInfo.residualDeltaP;
    info.residualMomentumKickNorm = momentumInfo.residualDeltaPNorm;
    info.dvParticleRawRms = momentumInfo.dvRawRms;
    info.dvParticleCorrectedRms = momentumInfo.dvCorrectedRms;
else
    info.momentumCorrection = empty_momentum_info(params, size(stateIn.v,1));
end

if logical(getp(params, 'virialStoreMaps', false))
    info.grid = G;
else
    info = rmfield_safe(info, {'Pkin','PvirEOS','PtotEOS','PvirDrive','Pdrive','dPdx','dPdy','duVirialX','duVirialY'});
end
end

function validate_state(state)
if ~isstruct(state) || ~isfield(state,'x') || ~isfield(state,'v') || size(state.x,2) ~= 2 || size(state.v,2) ~= 2 || size(state.x,1) ~= size(state.v,1)
    error('state must contain matching Np-by-2 x and v arrays.');
end
end

function T = virial_temperature_map(x, v, params, G, periodicX, periodicY)
Nx = params.Nx;
Ny = params.Ny;
dx = params.Lx / Nx;
dy = params.Ly / Ny;
xp = x(:,1);
yp = x(:,2);
if periodicX
    xp = mod(xp, params.Lx);
else
    xp = min(max(xp, 0), params.Lx - eps(params.Lx));
end
if periodicY
    yp = mod(yp, params.Ly);
else
    yp = min(max(yp, 0), params.Ly - eps(params.Ly));
end
ix = floor(xp / dx) + 1;
iy = floor(yp / dy) + 1;
ix = min(max(ix,1),Nx);
iy = min(max(iy,1),Ny);
cellId = iy + Ny*(ix-1);
UxVec = reshape(G.Ux.', [], 1);
UyVec = reshape(G.Uy.', [], 1);
NVec = reshape(G.N.', [], 1);
ux = UxVec(cellId);
uy = UyVec(cellId);
rel2 = (v(:,1)-ux).^2 + (v(:,2)-uy).^2;
sumRel2 = accumarray(cellId, rel2, [Nx*Ny, 1], @sum, 0);
dof = 2 * max(double(NVec)-1, 0);
Tvec = nan(Nx*Ny, 1);
mask = dof > 0;
Tvec(mask) = sumRel2(mask) ./ dof(mask);
T = reshape(Tvec, [Ny, Nx]).';
end

function rho = resolve_rho_eos_ref(params, state, rhoMean)
mode = lower(strrep(char(string(getp(params, 'virialRhoEOSRefMode', 'initial_physical_density'))), '-', '_'));
switch mode
    case {'initial_physical_density','initial','rho0'}
        Ly0 = getp(params, 'Ly0', getp(params, 'LyReference', params.Ly));
        rho = params.Nx * params.Ny * getp(params, 'gamma', size(state.x,1)/max(params.Nx*params.Ny,1)) / max(params.Lx * Ly0, eps);
    case {'current_uniform','uniform_now','rho_mean'}
        rho = rhoMean;
    case {'explicit','user'}
        rho = getp(params, 'virialRhoEOSRef', rhoMean);
    otherwise
        error('Unknown virialRhoEOSRefMode: %s', mode);
end
end

function rho = resolve_rho_uniform_now(params, state, rhoMean)
mode = lower(strrep(char(string(getp(params, 'virialRhoUniformMode', 'reference_gamma_current_volume'))), '-', '_'));
switch mode
    case {'reference_gamma_current_volume','current_volume','gamma_current_volume'}
        rho = params.Nx * params.Ny * getp(params, 'gamma', size(state.x,1)/max(params.Nx*params.Ny,1)) / max(params.Lx * params.Ly, eps);
    case {'particle_mean','rho_mean','actual'}
        rho = rhoMean;
    case {'explicit','user'}
        rho = getp(params, 'virialRhoUniformNow', rhoMean);
    otherwise
        error('Unknown virialRhoUniformMode: %s', mode);
end
end

function rhoKick = resolve_rho_kick(params, rhoLocal, rhoUniformNow)
mode = lower(strrep(char(string(getp(params, 'virialRhoKickMode', 'uniform_now'))), '-', '_'));
switch mode
    case {'uniform_now','current_uniform','constant'}
        rhoKick = rhoUniformNow .* ones(size(rhoLocal));
    case {'local','cell','rho_local'}
        rhoMin = getp(params, 'virialRhoKickMinFraction', 0.1) * max(rhoUniformNow, eps);
        rhoKick = max(rhoLocal, rhoMin);
    otherwise
        error('Unknown virialRhoKickMode: %s', mode);
end
end

function [dFdx, dFdy] = virial_gradient_boundary_aware(F, params, periodicX, periodicY)
dx = params.Lx / params.Nx;
dy = params.Ly / params.Ny;
[Nx, Ny] = size(F);
dFdx = zeros(Nx, Ny);
dFdy = zeros(Nx, Ny);
if periodicX
    dFdx = (circshift(F, [-1, 0]) - circshift(F, [1, 0])) / (2 * dx);
elseif Nx > 1
    dFdx(1,:) = (F(2,:) - F(1,:)) / dx;
    dFdx(Nx,:) = (F(Nx,:) - F(Nx-1,:)) / dx;
    if Nx > 2
        dFdx(2:Nx-1,:) = (F(3:Nx,:) - F(1:Nx-2,:)) / (2 * dx);
    end
end
if periodicY
    dFdy = (circshift(F, [0, -1]) - circshift(F, [0, 1])) / (2 * dy);
elseif Ny > 1
    dFdy(:,1) = (F(:,2) - F(:,1)) / dy;
    dFdy(:,Ny) = (F(:,Ny) - F(:,Ny-1)) / dy;
    if Ny > 2
        dFdy(:,2:Ny-1) = (F(:,3:Ny) - F(:,1:Ny-2)) / (2 * dy);
    end
end
end

function [duX, duY, info] = apply_du_limiter(duX, duY, params)
info = struct('enabled', false, 'duMax', Inf, 'nLimitedCells', 0, 'limitedCellFraction', 0.0, ...
    'duRmsBeforeLimiter', sqrt(mean(duX(:).^2 + duY(:).^2, 'omitnan')), ...
    'duRmsAfterLimiter', NaN);
if ~logical(getp(params, 'virialLimiterEnable', true))
    info.duRmsAfterLimiter = info.duRmsBeforeLimiter;
    return;
end
frac = getp(params, 'virialMaxDuFractionThermal', 0.1);
if ~isfinite(frac) || frac <= 0
    info.duRmsAfterLimiter = info.duRmsBeforeLimiter;
    return;
end
duMax = frac * sqrt(max(getp(params, 'kBT', 1.0), eps));
mag = sqrt(duX.^2 + duY.^2);
mask = mag > duMax & isfinite(mag);
if any(mask(:))
    scale = duMax ./ max(mag(mask), eps);
    duX(mask) = duX(mask) .* scale;
    duY(mask) = duY(mask) .* scale;
end
info.enabled = true;
info.duMax = duMax;
info.nLimitedCells = nnz(mask);
info.limitedCellFraction = nnz(mask) / max(numel(mask), 1);
info.duRmsAfterLimiter = sqrt(mean(duX(:).^2 + duY(:).^2, 'omitnan'));
end

function info = empty_info(params, Kvirial, betaVirial, enableDiag, enableKick)
info = struct();
info.enabled = false;
info.diagnosticsEnabled = enableDiag;
info.kickEnabled = enableKick;
info.Kvirial = Kvirial;
info.betaVirial = betaVirial;
info.rhoEOSRef = NaN;
info.rhoUniformNow = NaN;
info.rhoDriveRef = NaN;
info.rhoMean = NaN;
info.rhoDefectRms = NaN;
info.rhoDefectRelRms = NaN;
info.rhoMin = NaN;
info.rhoMax = NaN;
info.PkinMean = NaN;
info.PvirMean = NaN;
info.PtotMean = NaN;
info.PdriveMean = NaN;
info.PkinMin = NaN;
info.PkinMax = NaN;
info.PvirMin = NaN;
info.PvirMax = NaN;
info.PtotMin = NaN;
info.PtotMax = NaN;
info.gradPdriveRms = NaN;
info.gradPdriveMaxAbs = NaN;
info.duVirialRawRms = NaN;
info.duVirialAppliedRms = NaN;
info.duVirialAppliedMaxAbs = NaN;
info.duVirialOverThermalRms = NaN;
info.limitedCellFraction = NaN;
info.limitedCellCount = NaN;
info.duLimiterMax = NaN;
info.duRmsBeforeLimiter = NaN;
info.duRmsAfterLimiter = NaN;
info.duMaxFractionThermal = NaN;
info.limiter = struct();
info.momentumCorrection = empty_momentum_info(params, 0);
info.rawMomentumKick = [NaN NaN];
info.rawMomentumKickNorm = NaN;
info.residualMomentumKick = [NaN NaN];
info.residualMomentumKickNorm = NaN;
info.dvParticleRawRms = NaN;
info.dvParticleCorrectedRms = NaN;
end

function info = empty_momentum_info(params, Np)
info = struct();
info.stage = 'virial';
info.enabled = logical(getp(params, 'projectionMomentumCorrectionEnable', true));
info.applied = false;
info.mode = char(string(getp(params, 'projectionMomentumCorrectionMode', 'particle_global_exact')));
info.nParticles = Np;
info.nActiveParticles = Np;
info.rawDeltaP = [0 0];
info.rawDeltaPNorm = 0;
info.residualDeltaP = [0 0];
info.residualDeltaPNorm = 0;
info.dvRawRms = 0;
info.dvCorrectedRms = 0;
end

function s = rmfield_safe(s, names)
for i = 1:numel(names)
    if isfield(s, names{i})
        s = rmfield(s, names{i});
    end
end
end


function v = getf(s, name, defaultValue)
%GETF Local safe field getter. Kept here to make this standalone function
% independent of helper subfunctions defined in other piston scripts.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end

function v = getp(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
