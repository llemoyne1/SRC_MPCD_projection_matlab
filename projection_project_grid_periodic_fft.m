function proj = projection_project_grid_periodic_fft(Ux, Uy, params, varargin)
%PROJECTION_PROJECT_GRID_PERIODIC_FFT Periodic incompressible projection on a 2D grid.

if nargin < 3 || isempty(params)
    params = struct();
end
if ~isnumeric(Ux) || ~isnumeric(Uy) || ~isequal(size(Ux), size(Uy))
    error('Ux and Uy must be numeric arrays with identical size.');
end

rho = get_param_default(params, 'rho0', 1.0);
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "rho"
            rho = val;
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Lx = get_param_default(params, 'Lx', size(Ux, 1));
Ly = get_param_default(params, 'Ly', size(Ux, 2));
dt = get_param_default(params, 'dt', 1.0);
if rho <= 0 || dt <= 0 || Lx <= 0 || Ly <= 0
    error('rho, dt, Lx and Ly must be strictly positive.');
end

[Nx, Ny] = size(Ux);
kx1 = periodic_spectral_wavenumbers(Nx, Lx);
ky1 = periodic_spectral_wavenumbers(Ny, Ly);
[KX, KY] = ndgrid(kx1, ky1);
K2 = KX.^2 + KY.^2;

UxHat = fft2(double(Ux));
UyHat = fft2(double(Uy));
divHatBefore = 1i*KX.*UxHat + 1i*KY.*UyHat;
rhsHat = (rho/dt) * divHatBefore;

pHat = zeros(size(rhsHat));
mask = K2 > 0;
pHat(mask) = -rhsHat(mask) ./ K2(mask);
pHat(~mask) = 0;

gpxHat = 1i*KX .* pHat;
gpyHat = 1i*KY .* pHat;
dUx = real(ifft2(-(dt/rho) * gpxHat));
dUy = real(ifft2(-(dt/rho) * gpyHat));
UxProj = double(Ux) + dUx;
UyProj = double(Uy) + dUy;

divHatAfter = 1i*KX.*fft2(UxProj) + 1i*KY.*fft2(UyProj);
divBefore = real(ifft2(divHatBefore));
divAfter = real(ifft2(divHatAfter));

proj = struct();
proj.Ux = UxProj;
proj.Uy = UyProj;
proj.dUx = dUx;
proj.dUy = dUy;
proj.p = real(ifft2(pHat));
proj.divBefore = divBefore;
proj.divAfter = divAfter;
proj.rmsDivBefore = sqrt(mean(divBefore(:).^2, 'omitnan'));
proj.rmsDivAfter = sqrt(mean(divAfter(:).^2, 'omitnan'));
proj.maxAbsDivBefore = max(abs(divBefore(:)));
proj.maxAbsDivAfter = max(abs(divAfter(:)));
proj.meanDivBefore = mean(divBefore(:), 'omitnan');
proj.meanDivAfter = mean(divAfter(:), 'omitnan');
proj.rho = rho;
proj.dt = dt;
proj.Lx = Lx;
proj.Ly = Ly;
proj.Nx = Nx;
proj.Ny = Ny;
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
