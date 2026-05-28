function cellId = resamp_cell_ids_periodic(x, params, varargin)
%RESAMP_CELL_IDS_PERIODIC Return nearest MPCD cell ids for particle positions.
%
% The function name says periodic because the first resampling prototype is
% periodic, but periodicX/periodicY can still be toggled for reuse.

if size(x, 2) < 2
    error('x must be an Np-by-2 array.');
end
periodicX = true;
periodicY = true;
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "periodicx"
            periodicX = logical(val);
        case "periodicy"
            periodicY = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;

xp = x(:,1);
yp = x(:,2);
if periodicX
    xp = mod(xp, Lx);
else
    xp = min(max(xp, 0), Lx - eps(Lx));
end
if periodicY
    yp = mod(yp, Ly);
else
    yp = min(max(yp, 0), Ly - eps(Ly));
end

ix = floor(xp / dx) + 1;
iy = floor(yp / dy) + 1;
ix = min(max(ix, 1), Nx);
iy = min(max(iy, 1), Ny);
cellId = iy + Ny * (ix - 1);
end
