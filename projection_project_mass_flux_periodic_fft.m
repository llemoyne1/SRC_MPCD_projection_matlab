function proj = projection_project_mass_flux_periodic_fft(N, Ux, Uy, params, varargin)
%PROJECTION_PROJECT_MASS_FLUX_PERIODIC_FFT Project mass flux M=N*u in a periodic box.
%
% Supports Q9 mode 'relax_to_uniform_lowk', where the low-k divergence of
% mass flux is relaxed toward beta/dt*(N-gamma).  It also supports a pure
% conservative div(M)=0 mode.

if nargin < 4 || isempty(params)
    params = struct();
end
if ~isequal(size(N), size(Ux), size(Uy))
    error('N, Ux and Uy must have identical size.');
end

mode = char(string(get_param_default(params, 'massFluxProjectionMode', 'relax_to_uniform_lowk')));
beta = get_param_default(params, 'massFluxDensityRelaxationBeta', 0.0);
lowKMaxIndex = get_param_default(params, 'massFluxLowKMaxIndex', get_param_default(params, 'lowKMaxIndex', 2));
targetFilter = char(string(get_param_default(params, 'massFluxTargetFilter', 'lowpass_fft')));
minCellCount = get_param_default(params, 'massFluxMinCellCount', 1.0);
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "mode"
            mode = char(string(val));
        case {"relaxationbeta", "beta"}
            beta = val;
        case {"lowkmaxindex", "lowk"}
            lowKMaxIndex = val;
        case "targetfilter"
            targetFilter = char(string(val));
        otherwise
            error('Unknown option: %s', string(key));
    end
end
modeLower = lower(strrep(mode, '-', '_'));

Lx = get_param_default(params, 'Lx', size(Ux, 1));
Ly = get_param_default(params, 'Ly', size(Ux, 2));
dt = get_param_default(params, 'dt', 1.0);
gamma = get_param_default(params, 'gamma', mean(double(N(:)), 'omitnan'));
if dt <= 0 || Lx <= 0 || Ly <= 0
    error('dt, Lx and Ly must be strictly positive.');
end

[Nx, Ny] = size(Ux);
kx1 = periodic_spectral_wavenumbers(Nx, Lx);
ky1 = periodic_spectral_wavenumbers(Ny, Ly);
[KX, KY] = ndgrid(kx1, ky1);
K2 = KX.^2 + KY.^2;
nonzeroK = K2 > 0;
lowKMask = make_low_k_mask(Nx, Ny, lowKMaxIndex);
lowKMask(1,1) = false;

N = double(N);
Ux = double(Ux);
Uy = double(Uy);
Mx = N .* Ux;
My = N .* Uy;
MxHat = fft2(Mx);
MyHat = fft2(My);
divBeforeHat = 1i*KX.*MxHat + 1i*KY.*MyHat;
divBeforeFull = real(ifft2(divBeforeHat));

lowKCorrectionOnly = false;
switch modeLower
    case {"conservative", "divrho", "divrho_zero", "mass_conservative"}
        targetRequested = zeros(Nx, Ny);
        target = targetRequested;
        canonicalMode = 'conservative';
    case {"relax_to_uniform", "relax", "density_relax", "rho_relax"}
        targetRaw = (beta / dt) * (N - gamma);
        targetRaw = targetRaw - mean(targetRaw(:), 'omitnan');
        targetRequested = apply_mass_flux_target_filter(targetRaw, targetFilter, lowKMaxIndex);
        targetRequested = targetRequested - mean(targetRequested(:), 'omitnan');
        target = targetRequested;
        canonicalMode = 'relax_to_uniform';
    case {"relax_to_uniform_lowk", "relax_lowk", "density_relax_lowk", "rho_relax_lowk"}
        targetRaw = (beta / dt) * (N - gamma);
        targetRaw = targetRaw - mean(targetRaw(:), 'omitnan');
        targetRequested = apply_mass_flux_target_filter(targetRaw, 'lowpass_fft', lowKMaxIndex);
        targetRequested = targetRequested - mean(targetRequested(:), 'omitnan');
        target = targetRequested;
        canonicalMode = 'relax_to_uniform_lowk';
        targetFilter = 'lowpass_fft';
        lowKCorrectionOnly = true;
    otherwise
        error('Unknown massFluxProjectionMode: %s', mode);
end

targetHat = fft2(target);
rhsHatFull = divBeforeHat - targetHat;
rhsHatUsed = rhsHatFull;
if lowKCorrectionOnly
    rhsHatUsed(~lowKMask) = 0;
    rhsHatUsed(1,1) = 0;
end

lambdaHat = zeros(Nx, Ny);
lambdaHat(nonzeroK) = -rhsHatUsed(nonzeroK) ./ K2(nonzeroK);
lambdaHat(~nonzeroK) = 0;

dMxHat = -1i*KX .* lambdaHat;
dMyHat = -1i*KY .* lambdaHat;
MxProj = real(ifft2(MxHat + dMxHat));
MyProj = real(ifft2(MyHat + dMyHat));

divAfterHatFull = 1i*KX.*fft2(MxProj) + 1i*KY.*fft2(MyProj);
divAfterFull = real(ifft2(divAfterHatFull));
residualFull = divAfterFull - target;

if lowKCorrectionOnly
    divBeforeDiag = real(ifft2(mask_fft(divBeforeHat, lowKMask)));
    divAfterDiag = real(ifft2(mask_fft(divAfterHatFull, lowKMask)));
    residualDiag = divAfterDiag - target;
else
    divBeforeDiag = divBeforeFull;
    divAfterDiag = divAfterFull;
    residualDiag = residualFull;
end

Nsafe = max(N, minCellCount);
UxProj = MxProj ./ Nsafe;
UyProj = MyProj ./ Nsafe;
UxProj(N <= 0) = 0;
UyProj(N <= 0) = 0;
dUx = UxProj - Ux;
dUy = UyProj - Uy;

proj = struct();
proj.mode = canonicalMode;
proj.beta = beta;
proj.N = N;
proj.gamma = gamma;
proj.Mx = Mx;
proj.My = My;
proj.MxProjected = MxProj;
proj.MyProjected = MyProj;
proj.Ux = UxProj;
proj.Uy = UyProj;
proj.dUx = dUx;
proj.dUy = dUy;
proj.lambda = real(ifft2(lambdaHat));
proj.pi = proj.lambda / dt;
proj.targetDivMassRequested = targetRequested;
proj.targetDivMass = target;
proj.divMassBefore = divBeforeDiag;
proj.divMassAfter = divAfterDiag;
proj.divMassResidual = residualDiag;
proj.divMassBeforeFull = divBeforeFull;
proj.divMassAfterFull = divAfterFull;
proj.divMassResidualFull = residualFull;
proj.rmsDivMassBefore = sqrt(mean(divBeforeDiag(:).^2, 'omitnan'));
proj.rmsDivMassAfter = sqrt(mean(divAfterDiag(:).^2, 'omitnan'));
proj.rmsTargetDivMassRequested = sqrt(mean(targetRequested(:).^2, 'omitnan'));
proj.rmsTargetDivMass = sqrt(mean(target(:).^2, 'omitnan'));
proj.rmsDivMassResidual = sqrt(mean(residualDiag(:).^2, 'omitnan'));
proj.maxAbsDivMassBefore = max(abs(divBeforeDiag(:)));
proj.maxAbsDivMassAfter = max(abs(divAfterDiag(:)));
proj.maxAbsTargetDivMassRequested = max(abs(targetRequested(:)));
proj.maxAbsTargetDivMass = max(abs(target(:)));
proj.maxAbsDivMassResidual = max(abs(residualDiag(:)));
proj.rmsDivMassBeforeFull = sqrt(mean(divBeforeFull(:).^2, 'omitnan'));
proj.rmsDivMassAfterFull = sqrt(mean(divAfterFull(:).^2, 'omitnan'));
proj.rmsDivMassResidualFull = sqrt(mean(residualFull(:).^2, 'omitnan'));
proj.lowKCorrectionOnly = lowKCorrectionOnly;
proj.divMassReduction = proj.rmsDivMassResidual / max(proj.rmsDivMassBefore, eps);
proj.massDeltaBefore = sum(divBeforeFull(:));
proj.massDeltaAfter = sum(divAfterFull(:));
proj.massDeltaTargetRequested = sum(targetRequested(:));
proj.massDeltaTarget = sum(target(:));
proj.dt = dt;
proj.Lx = Lx;
proj.Ly = Ly;
proj.Nx = Nx;
proj.Ny = Ny;
proj.dx = Lx / Nx;
proj.dy = Ly / Ny;
proj.minCellCount = minCellCount;
proj.targetFilter = targetFilter;
proj.lowKMaxIndex = lowKMaxIndex;
proj.method = 'mass_flux_periodic_fft';
end

function B = mask_fft(A, mask)
B = A;
B(~mask) = 0;
B(1,1) = 0;
end

function target = apply_mass_flux_target_filter(targetRaw, targetFilter, lowKMaxIndex)
filterLower = lower(strrep(char(string(targetFilter)), '-', '_'));
[Nx, Ny] = size(targetRaw);
switch filterLower
    case {"none", "off", "identity", "raw"}
        target = double(targetRaw);
    case {"lowpass_fft", "lowk", "fft_lowk"}
        mask = make_low_k_mask(Nx, Ny, lowKMaxIndex);
        targetHat = fft2(double(targetRaw));
        targetHat(~mask) = 0;
        target = real(ifft2(targetHat));
    otherwise
        error('Unknown massFluxTargetFilter: %s', targetFilter);
end
end

function mask = make_low_k_mask(Nx, Ny, kmax)
ix = [0:floor(Nx/2), -ceil(Nx/2)+1:-1];
iy = [0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix) > Nx, ix = ix(1:Nx); end
if numel(iy) > Ny, iy = iy(1:Ny); end
[KX, KY] = ndgrid(ix, iy);
mask = sqrt(double(KX).^2 + double(KY).^2) <= kmax;
end

function k = periodic_spectral_wavenumbers(N, L)
if mod(N, 2) == 0
    k = (2*pi/L) * [0:(N/2-1), 0, (-N/2+1):-1];
else
    k = (2*pi/L) * [0:((N-1)/2), -((N-1)/2):-1];
end
end

function val = get_param_default(params, name, defaultValue)
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultValue;
end
end
