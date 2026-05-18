function proj = projection_project_mass_flux_channel_operator(N, Ux, Uy, params, varargin)
%PROJECTION_PROJECT_MASS_FLUX_CHANNEL_OPERATOR Dispatch Q9 mass-flux projection.
%
% Default behavior is intentionally unchanged: unless
% params.massFluxProjectionOperator requests the new general finite-volume
% operator, this wrapper calls the legacy periodic-x/bounded-y operator.
%
% New operator aliases:
%   'general_bc'
%   'finite_volume_general_bc'
%   'elliptic_general_bc'
%   'fv_general_bc'
%
% The new operator is a validation path for a compact, BC-aware, matrix-free
% future C++/OpenMP/GPU implementation.  The MATLAB version may still use a
% sparse matrix internally for robustness and diagnostics.

operator = char(string(get_param_default(params, 'massFluxProjectionOperator', 'legacy_periodic_x_neumann_y')));
operatorLower = lower(strrep(operator, '-', '_'));

switch operatorLower
    case {'legacy', 'legacy_periodic_x_neumann_y', 'periodic_x_neumann_y', ...
          'ddt_periodic_x_bounded_y', 'collocated_legacy'}
        proj = projection_project_mass_flux_periodic_x_neumann_y(N, Ux, Uy, params, varargin{:});
    case {'general_bc', 'finite_volume_general_bc', 'elliptic_general_bc', ...
          'fv_general_bc', 'general', 'finite_volume'}
        proj = projection_project_mass_flux_general_bc(N, Ux, Uy, params, varargin{:});
    otherwise
        error('Unknown massFluxProjectionOperator: %s', operator);
end
end

function val = get_param_default(params, name, defaultValue)
if nargin < 1 || isempty(params) || ~isstruct(params)
    val = defaultValue;
elseif isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = defaultValue;
end
end
