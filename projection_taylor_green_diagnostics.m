function tg = projection_taylor_green_diagnostics(G, params)
%PROJECTION_TAYLOR_GREEN_DIAGNOSTICS Structure metrics for periodic Taylor-Green flow.
Ux = double(G.Ux); Uy = double(G.Uy); N = double(G.N);
[Nx, Ny] = size(Ux);
Lx = getf(params, 'Lx', G.Lx); Ly = getf(params, 'Ly', G.Ly);
dx = Lx/Nx; dy = Ly/Ny;
mx = getf(params, 'taylorGreenModeX', 1); my = getf(params, 'taylorGreenModeY', 1);
kx = 2*pi*mx/Lx; ky = 2*pi*my/Ly; k2mode = kx^2 + ky^2;
A0 = getf(params, 'taylorGreenAmplitude', getf(params, 'taylorGreenInitialAmplitude', NaN));

xc = ((0:Nx-1)+0.5)*dx; yc = ((0:Ny-1)+0.5)*dy;
[X,Y] = ndgrid(xc,yc);
phiX = sin(kx*X).*cos(ky*Y);
phiY = -cos(kx*X).*sin(ky*Y);
ampYFactor = kx/max(ky, eps);
phiYScaled = ampYFactor*phiY;

num = sum(Ux(:).*phiX(:), 'omitnan') + sum(Uy(:).*phiYScaled(:), 'omitnan');
den = sum(phiX(:).^2, 'omitnan') + sum(phiYScaled(:).^2, 'omitnan');
modeAmplitude = num/max(den, eps);
modeUxAmplitude = sum(Ux(:).*phiX(:), 'omitnan')/max(sum(phiX(:).^2, 'omitnan'), eps);
modeUyAmplitude = sum(Uy(:).*phiYScaled(:), 'omitnan')/max(sum(phiYScaled(:).^2, 'omitnan'), eps);
UmodeX = modeAmplitude*phiX;
UmodeY = modeAmplitude*phiYScaled;
modeEnergy = 0.5*mean(UmodeX(:).^2 + UmodeY(:).^2, 'omitnan');
totalHydroEnergy = 0.5*mean(Ux(:).^2 + Uy(:).^2, 'omitnan');
residualEnergy = 0.5*mean((Ux(:)-UmodeX(:)).^2 + (Uy(:)-UmodeY(:)).^2, 'omitnan');
modeCoherence = modeEnergy/max(totalHydroEnergy, eps);

kx1 = periodic_spectral_wavenumbers(Nx, Lx); ky1 = periodic_spectral_wavenumbers(Ny, Ly);
[KX,KY] = ndgrid(kx1,ky1);
omega = real(ifft2(1i*KX.*fft2(Uy) - 1i*KY.*fft2(Ux)));
enstrophy = 0.5*mean(omega(:).^2, 'omitnan');
omegaRms = sqrt(mean(omega(:).^2, 'omitnan'));

UxHat = fft2(Ux); UyHat = fft2(Uy);
energySpec = 0.5*(abs(UxHat).^2 + abs(UyHat).^2)/(Nx*Ny)^2;
modeMask = false(Nx,Ny);
modeMask = mark_mode(modeMask, mx, my);
modeMask = mark_mode(modeMask, -mx, my);
modeMask = mark_mode(modeMask, mx, -my);
modeMask = mark_mode(modeMask, -mx, -my);
lowKMaxIndex = getf(params, 'lowKMaxIndex', getf(params, 'massFluxLowKMaxIndex', 2));
lowKMask = make_low_k_mask(Nx, Ny, lowKMaxIndex); lowKMask(1,1)=false;
highKCut = getf(params, 'taylorGreenHighKCut', max(4, min(Nx,Ny)/4));
highKMask = make_high_k_mask(Nx, Ny, highKCut);
nonzeroMask = true(Nx,Ny); nonzeroMask(1,1)=false;
spectralEnergyTotal = sum(energySpec(nonzeroMask), 'omitnan');
spectralEnergyLowK = sum(energySpec(lowKMask), 'omitnan');
spectralEnergyHighK = sum(energySpec(highKMask), 'omitnan');
spectralEnergyTGMode = sum(energySpec(modeMask), 'omitnan');

projDiv = projection_project_grid_periodic_fft(Ux, Uy, params);

gamma = getf(params, 'gamma', mean(N(:), 'omitnan'));
rhoRel = N/max(gamma, eps) - 1;
rhoRel = rhoRel - mean(rhoRel(:), 'omitnan');
rhoHat = fft2(rhoRel);
rhoSpec = abs(rhoHat).^2/(Nx*Ny)^2;
densityLowKEnergy = sum(rhoSpec(lowKMask), 'omitnan');
densityRelRMS = sqrt(mean(rhoRel(:).^2, 'omitnan'));

tg = struct();
tg.modeAmplitude = modeAmplitude;
tg.modeUxAmplitude = modeUxAmplitude;
tg.modeUyAmplitude = modeUyAmplitude;
tg.targetAmplitude = A0;
tg.amplitudeRatioToTarget = modeAmplitude/max(abs(A0), eps);
tg.modeEnergy = modeEnergy;
tg.totalHydroEnergy = totalHydroEnergy;
tg.residualEnergy = residualEnergy;
tg.modeCoherence = modeCoherence;
tg.omega = omega;
tg.enstrophy = enstrophy;
tg.omegaRms = omegaRms;
tg.spectralEnergyTotal = spectralEnergyTotal;
tg.spectralEnergyLowK = spectralEnergyLowK;
tg.spectralEnergyHighK = spectralEnergyHighK;
tg.spectralEnergyTGMode = spectralEnergyTGMode;
tg.lowKEnergyFractionVelocity = spectralEnergyLowK/max(spectralEnergyTotal, eps);
tg.highKEnergyFractionVelocity = spectralEnergyHighK/max(spectralEnergyTotal, eps);
tg.tgModeEnergyFractionVelocity = spectralEnergyTGMode/max(spectralEnergyTotal, eps);
tg.rmsDiv = projDiv.rmsDivBefore;
tg.maxAbsDiv = projDiv.maxAbsDivBefore;
tg.modeX = mx; tg.modeY = my;
tg.kx = kx; tg.ky = ky; tg.k2mode = k2mode;
tg.dx = dx; tg.dy = dy;
tg.densityLowKEnergy = densityLowKEnergy;
tg.densityRelRMS = densityRelRMS;
tg.densityRel = rhoRel;
end

function mask = mark_mode(mask,kx,ky)
[Nx,Ny]=size(mask); mask(mode_to_index(kx,Nx), mode_to_index(ky,Ny)) = true;
end
function id = mode_to_index(k,N), k=mod(k,N); id=k+1; end
function mask = make_low_k_mask(Nx,Ny,kmax)
ix=[0:floor(Nx/2), -ceil(Nx/2)+1:-1]; iy=[0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix)>Nx, ix=ix(1:Nx); end; if numel(iy)>Ny, iy=iy(1:Ny); end
[KX,KY]=ndgrid(ix,iy); mask=sqrt(double(KX).^2+double(KY).^2)<=kmax;
end
function mask = make_high_k_mask(Nx,Ny,kcut)
ix=[0:floor(Nx/2), -ceil(Nx/2)+1:-1]; iy=[0:floor(Ny/2), -ceil(Ny/2)+1:-1];
if numel(ix)>Nx, ix=ix(1:Nx); end; if numel(iy)>Ny, iy=iy(1:Ny); end
[KX,KY]=ndgrid(ix,iy); mask=sqrt(double(KX).^2+double(KY).^2)>=kcut;
end
function k = periodic_spectral_wavenumbers(N,L)
if mod(N,2)==0, k=(2*pi/L)*[0:(N/2-1), 0, (-N/2+1):-1]; else, k=(2*pi/L)*[0:((N-1)/2), -((N-1)/2):-1]; end
end
function v=getf(s,name,defaultValue)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=defaultValue; end
end
