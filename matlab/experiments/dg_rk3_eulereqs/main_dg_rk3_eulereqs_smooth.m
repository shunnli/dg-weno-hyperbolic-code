clc;
clear;
close all;

cd(fileparts(mfilename('fullpath')));
addpath(genpath('../../base'))
addpath('../common')
set_groot_plot_style

alpha = 1;
beta = 0.2;
omega = 1;
phi = 0;
u0_ic = 1;
p0_ic = 1;

xleft = -pi;
xright = pi;

pk = 2; % Pk
gk = 5; % Gauss Points

dim = 3;

basis = MatLegendre(pk + 1);
basis_dx = MatLegendreDx(pk + 1);

%% order test

t1 = 0.5;
nxlist = [10, 20, 40, 80, 160, 320, 640];
n = numel(nxlist);
errors = zeros(3, n);

for w = 1:n
    [x, dx] = mesh_init_1d(xleft, xright, nxlist(w));

    init_func = @(s) eulereqs_sin_exact(s, 0, alpha, beta, omega, phi, u0_ic, p0_ic);
    uh0 = dg_projection_eqs(init_func, x, dx, pk, gk, basis, dim);
    uh = dg_rk3_scheme_eqs(uh0, dx, t1, @eulereqs_f, @eulereqs_fhat_LF, @eulereqs_get_alpha, pk, gk, basis, basis_dx, dim, false, false);

    u_exact_func = @(s) eulereqs_sin_exact(s, t1, alpha, beta, omega, phi, u0_ic, p0_ic);

    errors(:, w) = dg_errors_eqs(uh, u_exact_func, x, dx, pk, gk, basis, dim, @eulereqs_trans2raw);
end

show_results(nxlist, errors(1, :), errors(2, :), errors(3, :));

%% plot test

t2 = 0.6;
nxlist2 = [20, 80];

nx_ref = 320;
[x_ref, dx_ref] = mesh_init_1d(xleft, xright, nx_ref);

exact_ref_func = @(s) eulereqs_sin_exact(s, t2, alpha, beta, omega, phi, u0_ic, p0_ic);
% u_exact_ref = dg_projection_eqs(exact_ref_func, x_ref, dx_ref, pk, gk, basis, dim);
[rho_values_ref, ~, ~] = exact_ref_func(x_ref);

figure();
hold('on');
plot(x_ref, rho_values_ref, DisplayName = 'rho-ref')

for w = 1:2
    [x, dx] = mesh_init_1d(xleft, xright, nxlist2(w));

    init_func = @(s) eulereqs_sin_exact(s, 0, alpha, beta, omega, phi, u0_ic, p0_ic);
    uh0 = dg_projection_eqs(init_func, x, dx, pk, gk, basis, dim);
    uh = dg_rk3_scheme_eqs(uh0, dx, t2, @eulereqs_f, @eulereqs_fhat_LF, @eulereqs_get_alpha, pk, gk, basis, basis_dx, dim, false, false);

    v = basis.eval(0, pk + 1); % column vector

    plot(x, v' * uh(1:(pk + 1), :), DisplayName = sprintf('rhoh (n=%d)', nxlist2(w)))
end

hold('off');
box('on');
legend(Location = 'best');
exportgraphics(gcf, 'result_smooth.png', Resolution = 600);

function [v1, v2, v3] = eulereqs_sin_exact(x, t, alpha, beta, omega, phi, u0, p0)
    % v = [v1, v2, v3] = [rho, rho u, E]

    rho = alpha + beta * sin(omega * (x - u0 * t) + phi);
    u = u0 * ones(size(x));
    p = p0 * ones(size(x));

    [v1, v2, v3] = eulereqs_trans2cv(rho, u, p);
end
