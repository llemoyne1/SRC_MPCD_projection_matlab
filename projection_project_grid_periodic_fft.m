function proj = projection_project_grid_periodic_fft(Ux, Uy, params, varargin)
%PROJECTION_PROJECT_GRID_PERIODIC_FFT Periodic incompressible projection on a 2D grid.
%
%   proj = projection_project_grid_periodic_fft(Ux, Uy, params)
%
% This is a deliberately small, grid-only prototype for the future
% SRC/MPCD pressure-projection branch. It projects a cell-centered velocity
% field (Ux, Uy) onto a divergence-free field on a fully periodic domain.
%
% The pressure is defined by:
%       laplacian(p) = rho/dt * div(u*)
% and the corrected velocity is:
%       u = u* - dt/rho * grad(p).
%
% The FFT implementation is useful as a reference test because the periodic
% projection should reduce spectral divergence to round-off error.
%
% Required params fields:
%   Lx, Ly, dt
%
% Optional params/arguments:
%   params.rho0 or name-value 'rho'   default 1
%
% Output fields:
%   Ux, Uy              projected velocity
%   dUx, dUy            velocity correction added to input field
%   p                   pressure field, zero mean
%   divBefore, divAfter spectral divergence before/after projection
%   rmsDivBefore, rmsDivAfter
%   rho, dt, Lx, Ly

if nargin < 3 || isempty(params)
    params = struct();
end
if ~isnumeric(Ux) || ~isnumeric(Uy) || ~isequal(size(Ux), size(Uy))
    error('Ux and Uy must be numeric arrays with identical size.');
end

rho = get_param_default(params, 'rho0', 1.0);
for k = 1:2:numel(varargin)
    key = varargin{k};
    val = varargin{k+1};
    switch lower(string(key))
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

kx1 = (2*pi/Lx) * [0:floor((Nx-1)/2), -ceil((Nx-1)/2):-1];
ky1 = (2*pi/Ly) * [0:floor((Ny-1)/2), -ceil((Ny-1)/2):-1];
[KX, KY] = ndgrid(kx1, ky1);
K2 = KX.^2 + KY.^2;

UxHat = fft2(Ux);
UyHat = fft2(Uy);

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
UxProj = Ux + dUx;
UyProj = Uy + dUy;

UxProjHat = fft2(UxProj);
UyProjHat = fft2(UyProj);
divHatAfter = 1i*KX.*UxProjHat + 1i*KY.*UyProjHat;

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
proj.rmsDivBefore = sqrt(mean(divBefore(:).^2));
proj.rmsDivAfter = sqrt(mean(divAfter(:).^2));
proj.maxAbsDivBefore = max(abs(divBefore(:)));
proj.maxAbsDivAfter = max(abs(divAfter(:)));
proj.rho = rho;
proj.dt = dt;
proj.Lx = Lx;
proj.Ly = Ly;
proj.Nx = Nx;
proj.Ny = Ny;
end

function val = get_param_default(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultValue;
end
end
