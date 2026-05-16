function [Ux, Uy] = projection_taylor_green_mode_at_points(x, params, amplitude)
%PROJECTION_TAYLOR_GREEN_MODE_AT_POINTS Evaluate divergence-free TG mode at points.
if nargin < 3 || isempty(amplitude)
    amplitude = get_param(params, 'taylorGreenAmplitude', get_param(params, 'taylorGreenInitialAmplitude', 0.1));
end
mx = get_param(params, 'taylorGreenModeX', 1);
my = get_param(params, 'taylorGreenModeY', 1);
kx = 2*pi*mx / params.Lx;
ky = 2*pi*my / params.Ly;
Ay = amplitude * kx / max(ky, eps);
Ux = amplitude * sin(kx*x(:,1)) .* cos(ky*x(:,2));
Uy = -Ay * cos(kx*x(:,1)) .* sin(ky*x(:,2));
end

function v = get_param(params, name, defaultValue)
if isfield(params, name) && ~isempty(params.(name))
    v = params.(name);
else
    v = defaultValue;
end
end
