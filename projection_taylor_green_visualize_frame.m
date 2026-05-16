function projection_taylor_green_visualize_frame(stateOrG, params, step, varargin)
%PROJECTION_TAYLOR_GREEN_VISUALIZE_FRAME Live visualization for forced TG.
figId = getf(params, 'visualFigureId', 410);
figName = getf(params, 'visualFigureName', 'Forced Taylor-Green');
titleSuffix = getf(params, 'visualTitleSuffix', '');
isMean = false;
saveFrame = false;
frameDir = getf(params, 'visualFrameDir', 'frames_tg');
framePrefix = getf(params, 'visualFramePrefix', 'tg');
for k=1:2:numel(varargin)
    key=lower(string(varargin{k})); val=varargin{k+1};
    switch key
        case "ismean", isMean=logical(val);
        case "saveframe", saveFrame=logical(val);
        case "framedir", frameDir=char(string(val));
        case "frameprefix", framePrefix=char(string(val));
    end
end

if isstruct(stateOrG) && isfield(stateOrG,'Ux')
    G = stateOrG; x = []; v = [];
else
    x = stateOrG.x; v = stateOrG.v;
    G = projection_deposit_particles_to_grid(x, v, params, 'periodicX', true, 'periodicY', true, 'minCount', 1);
end
tg = projection_taylor_green_diagnostics(G, params);
[Nx,Ny] = size(G.Ux); dx=params.Lx/Nx; dy=params.Ly/Ny;
xc=((0:Nx-1)+0.5)*dx; yc=((0:Ny-1)+0.5)*dy; [X,Y]=ndgrid(xc,yc);

figure(figId); clf;
set(gcf, 'Name', figName, 'Color', 'w');
kind = 'instantaneous'; if isMean, kind = 'running mean'; end
sgtitle(sprintf('%s TG %s step %d, t=%.3g %s | A=%.4g coh=%.3g enst=%.3g', ...
    kind, upper(getf(params,'method','')), step, step*params.dt, titleSuffix, tg.modeAmplitude, tg.modeCoherence, tg.enstrophy));

subplot(2,3,1);
if isempty(x)
    imagesc(xc,yc,G.N'); axis xy equal tight; colorbar; title('N'); xlabel('x'); ylabel('y');
else
    speedP = sqrt(sum(v.^2,2));
    maxP = getf(params,'visualMaxParticles',8000);
    if size(x,1)>maxP
        idx = randperm(size(x,1), maxP);
    else
        idx = 1:size(x,1);
    end
    scatter(x(idx,1), x(idx,2), 4, speedP(idx), 'filled'); axis equal tight; xlim([0 params.Lx]); ylim([0 params.Ly]); colorbar; title('particles |v|'); xlabel('x'); ylabel('y');
end

subplot(2,3,2);
rhoRel = G.N / max(params.gamma,eps) - 1;
imagesc(xc,yc,rhoRel'); axis xy equal tight; colorbar; title('density N/\gamma - 1'); xlabel('x'); ylabel('y');
cl = getf(params,'visualDensityCLim',NaN); if isfinite(cl), caxis([-cl cl]); end

subplot(2,3,3);
speed = sqrt(G.Ux.^2+G.Uy.^2);
imagesc(xc,yc,speed'); axis xy equal tight; hold on; colorbar; title('|u| + quiver'); xlabel('x'); ylabel('y');
strideX=max(1,getf(params,'visualQuiverStrideX',max(1,round(Nx/18)))); strideY=max(1,getf(params,'visualQuiverStrideY',max(1,round(Ny/12))));
ix=1:strideX:Nx; iy=1:strideY:Ny;
quiver(X(ix,iy),Y(ix,iy),G.Ux(ix,iy),G.Uy(ix,iy),getf(params,'visualQuiverScale',1.2),'k'); hold off;

subplot(2,3,4);
omega=tg.omega; imagesc(xc,yc,omega'); axis xy equal tight; colorbar; title('\omega = \partial_x u_y - \partial_y u_x'); xlabel('x'); ylabel('y');
omegaClim = getf(params,'visualOmegaCLim',NaN); if isfinite(omegaClim), caxis([-omegaClim omegaClim]); else, c=max(abs(prctile(omega(:),[2 98]))); if isfinite(c)&&c>0, caxis([-c c]); end; end

subplot(2,3,5);
resSpeed = sqrt((G.Ux - tg.modeAmplitude*mode_phi_x(X,Y,params)).^2 + (G.Uy - tg.modeAmplitude*mode_phi_y_scaled(X,Y,params)).^2);
imagesc(xc,yc,resSpeed'); axis xy equal tight; colorbar; title('|u - u_{TG}|'); xlabel('x'); ylabel('y');

subplot(2,3,6);
bar([tg.modeEnergy, tg.residualEnergy, tg.spectralEnergyHighK]);
set(gca,'XTickLabel',{'mode','residual','high-k'}); title(sprintf('coh=%.3g, high-k=%.3g', tg.modeCoherence, tg.highKEnergyFractionVelocity)); ylabel('energy'); grid on;
drawnow;

if saveFrame
    if ~exist(frameDir,'dir'), mkdir(frameDir); end
    exportgraphics(gcf, fullfile(frameDir, sprintf('%s_%07d.png', framePrefix, step)));
end
end

function phi = mode_phi_x(X,Y,params)
kx=2*pi*getf(params,'taylorGreenModeX',1)/params.Lx; ky=2*pi*getf(params,'taylorGreenModeY',1)/params.Ly;
phi=sin(kx*X).*cos(ky*Y);
end
function phi = mode_phi_y_scaled(X,Y,params)
kx=2*pi*getf(params,'taylorGreenModeX',1)/params.Lx; ky=2*pi*getf(params,'taylorGreenModeY',1)/params.Ly;
phi=-(kx/max(ky,eps))*cos(kx*X).*sin(ky*Y);
end
function v=getf(s,name,defaultValue)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=defaultValue; end
end
