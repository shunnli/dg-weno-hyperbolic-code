clc;
clear;
close all;

cd(fileparts(mfilename('fullpath')));
addpath(genpath('../../base'))
addpath('../common')
set_groot_plot_style

xleft = -2;
xright = 2;
dim = 3;
gamma = 1.4;

pk = 2; % Pk
gk = 5; % Gauss Points

limiter = tvd();

init_func = @(x) eulereqs_riemann_init(x, 1, 0, 1, 0.125, 0, 0.1);
exact_func = @(x, t) euler_riemann_exact(1, 0, 1, 0.125, 0, 0.1, gamma, x, 0, t);

% init_func = @(x) eulereqs_riemann_init(x, 0.445, 0.698, 3.528, 0.5, 0, 0.571);
% exact_func = @(x,t) euler_riemann_exact(0.445, 0.698, 3.528, 0.5, 0, 0.571, gamma, x, 0, t);

basis = MatLegendre(pk + 1);
basis_dx = MatLegendreDx(pk + 1);

%% plot test

t = 0.6;
nxlist = [80, 160];

for w = 1:2
    [x, dx] = mesh_init_1d(xleft, xright, nxlist(w));

    uh0 = dg_projection_eqs(init_func, x, dx, pk, gk, basis, dim);
    uh = dg_rk3_scheme_eqs(uh0, dx, t, @eulereqs_f, @eulereqs_fhat_LF, @eulereqs_get_alpha, pk, gk, basis, basis_dx, dim, true, limiter);

    vc = basis.eval(0, pk + 1); % column vector

    v1 = vc' * uh(1:(pk + 1), :);
    v2 = vc' * uh(pk + 2:2 * (pk + 1), :);
    v3 = vc' * uh(2 * (pk + 1) + 1:3 * (pk + 1), :);

    [rho_values, u_values, p_values] = eulereqs_trans2raw(v1, v2, v3);
    e_values = 1 / (gamma - 1) * p_values ./ rho_values;

    [rho_values_ref, u_values_ref, p_values_ref, ~] = exact_func(x, t);
    e_values_ref = 1 / (gamma - 1) * p_values_ref ./ rho_values_ref;

    primitive = {rho_values, u_values, p_values, e_values};
    primitive_ref = {rho_values_ref, u_values_ref, p_values_ref, e_values_ref};
    names = {'Density', 'Velocity', 'Pressure', 'Internal Energy'};

    figure();

    for s = 1:4
        subplot(2, 2, s);
        hold('on');

        q_ref = primitive_ref{s};
        plot(x, q_ref, 'b', DisplayName = 'u-ref');

        q = primitive{s};
        plot(x, q, DisplayName = sprintf('uh (n=%d)', nxlist(w)));

        qmax = max(q);
        qmin = min(q);
        qdiff = qmax - qmin;
        ylim([qmin - 0.1 * qdiff, qmax + 0.1 * qdiff]);

        hold('off');
        box('on');
        title(names{s});
        % legend(Location = 'best');
    end

    exportgraphics(gcf, sprintf('result_riemann_%d.png', w), Resolution = 600);
end

function [v1, v2, v3] = eulereqs_riemann_init(x, rho_l, u_l, p_l, rho_r, u_r, p_r)
    % v = [v1, v2, v3] = [rho, rho u, E]

    condition = x < 0;

    rho = zeros(size(x));
    u = zeros(size(x));
    p = zeros(size(x));

    rho(condition) = rho_l;
    u(condition) = u_l;
    p(condition) = p_l;

    rho(~condition) = rho_r;
    u(~condition) = u_r;
    p(~condition) = p_r;

    [v1, v2, v3] = eulereqs_trans2cv(rho, u, p);
end
