function cellId = resamp_cell_ids_periodic(x, params, varargin)
%RESAMP_CELL_IDS_PERIODIC Return nearest MPCD cell ids for particle positions.
%
%   CELLID = RESAMP_CELL_IDS_PERIODIC(X, PARAMS) maps the N-by-2 particle
%   position array X to linear cell ids using the convention
%
%       cellId = iy + Ny*(ix-1).
%
%   Name/value options:
%       periodicX      true/false, default true
%       periodicY      true/false, default true
%       activeMask     optional N-by-1 logical mask.  If provided, only
%                      active entries are mapped and inactive entries get
%                      inactiveValue.
%       inactiveValue  default 0
%
%   The activeMask option is accepted so callers operating on pooled states
%   can pass full state arrays without first compacting them.  Older callers
%   that pass already-compacted arrays are unchanged.

if size(x, 2) < 2
    error('x must be an Np-by-2 array.');
end

periodicX = true;
periodicY = true;
activeMask = [];
inactiveValue = 0;

if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    val = varargin{k+1};
    switch key
        case 'periodicx'
            periodicX = logical(val);
        case 'periodicy'
            periodicY = logical(val);
        case {'activemask','active'}
            activeMask = logical(val(:));
        case {'inactivevalue','inactivecellid'}
            inactiveValue = val;
        otherwise
            error('Unknown option: %s', key);
    end
end

Np = size(x,1);
if isempty(activeMask)
    activeMask = true(Np,1);
elseif numel(activeMask) ~= Np
    error('activeMask must have one entry per particle position.');
end

cellId = inactiveValue .* ones(Np,1);
if ~any(activeMask)
    return;
end

Lx = params.Lx;
Ly = params.Ly;
Nx = params.Nx;
Ny = params.Ny;
dx = Lx / Nx;
dy = Ly / Ny;

xp = x(activeMask,1);
yp = x(activeMask,2);

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
cellId(activeMask) = iy + Ny * (ix - 1);
end
