function [state, info] = projection_initialize_particles(params, varargin)
%PROJECTION_INITIALIZE_PARTICLES Initialize particles for projection prototypes.
%
%   [state, info] = projection_initialize_particles(params)
%
% Supported population modes:
%   'random_uniform'  : Np=round(gamma*Nx*Ny) particles uniformly in domain.
%   'exact_per_cell'  : exactly gamma particles uniformly placed inside each
%                       cell. gamma must be an integer in this simple mode.
%
% The velocity initialization is intentionally simple: Maxwellian velocities at
% params.kBT followed by optional removal of the global mean velocity. This
% keeps the first density-homogeneity tests focused on the population field.

if nargin < 1 || isempty(params)
    error('params is required.');
end

mode = get_param(params, 'initialPopulationMode', 'random_uniform');
zeroGlobalMean = get_param(params, 'initialVelocityZeroGlobalMean', true);

for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "initialpopulationmode"
            mode = char(val);
        case "initialvelocityzeroglobalmean"
            zeroGlobalMean = logical(val);
        otherwise
            error('Unknown option: %s', string(key));
    end
end

mode = lower(strrep(char(mode), '-', '_'));
Nx = params.Nx;
Ny = params.Ny;
Lx = params.Lx;
Ly = params.Ly;
gamma = params.gamma;
kBT = params.kBT;

switch mode
    case {'random_uniform','random','uniform'}
        Np = round(gamma * Nx * Ny);
        x = [Lx * rand(Np, 1), Ly * rand(Np, 1)];
        canonicalMode = 'random_uniform';

    case {'exact_per_cell','uniform_by_cell','homogeneous_per_cell','cell_exact'}
        gammaInt = round(gamma);
        if abs(gamma - gammaInt) > 1e-12
            error(['initialPopulationMode=exact_per_cell requires integer gamma in this simple prototype. ', ...
                   'Received gamma=%.16g.'], gamma);
        end
        if gammaInt < 0
            error('gamma must be non-negative.');
        end
        Np = gammaInt * Nx * Ny;
        x = zeros(Np, 2);
        dx = Lx / Nx;
        dy = Ly / Ny;
        cursor = 0;
        for ix = 0:Nx-1
            for iy = 0:Ny-1
                ids = cursor + (1:gammaInt);
                if gammaInt > 0
                    x(ids, 1) = (ix + rand(gammaInt, 1)) * dx;
                    x(ids, 2) = (iy + rand(gammaInt, 1)) * dy;
                end
                cursor = cursor + gammaInt;
            end
        end
        canonicalMode = 'exact_per_cell';

    otherwise
        error('Unknown initialPopulationMode: %s', mode);
end

v = sqrt(kBT) * randn(Np, 2);
if zeroGlobalMean && Np > 0
    v(:, 1) = v(:, 1) - mean(v(:, 1));
    v(:, 2) = v(:, 2) - mean(v(:, 2));
end

state = struct();
state.x = x;
state.v = v;

pop = projection_population_diagnostics(state.x, params, 'periodicX', true, 'periodicY', false);
info = struct();
info.mode = canonicalMode;
info.Np = Np;
info.gamma = gamma;
info.zeroGlobalMean = zeroGlobalMean;
info.population = pop;
info.initialStdN = pop.stdN;
info.initialMeanN = pop.meanN;
info.initialMinN = pop.minN;
info.initialMaxN = pop.maxN;
info.initialEmptyCells = pop.nEmptyCells;
info.initialOutBand20 = out_band_fraction(pop.N, gamma, 0.20);
info.initialMaxAbsDeltaFromGamma = max(abs(double(pop.N(:)) - gamma));
info.initialRmsRel = sqrt(mean((double(pop.N(:))/gamma - 1).^2, 'omitnan'));
end

function val = get_param(params, name, defaultVal)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultVal;
end
end

function f = out_band_fraction(N, gamma, bandFraction)
lo = gamma * (1 - bandFraction);
hi = gamma * (1 + bandFraction);
n = double(N(:));
f = nnz(n < lo | n > hi) / numel(n);
end
