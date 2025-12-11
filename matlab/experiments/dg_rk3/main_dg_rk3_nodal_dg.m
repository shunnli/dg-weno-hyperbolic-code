clc;
clear;
close all;

cd(fileparts(mfilename('fullpath')));
addpath(genpath('../../base'))
addpath('../common')
set_groot_plot_style

alpha = 0.5;
beta = 1;
xleft = -pi;
xright = pi;

pk = 2; % Pk
gk = 7; % Gauss Points

[points, ~] = gauss_lobatto(pk + 1);
basis = MatLagrange(points);
basis_dx = MatLagrangeDx(points);

% limiter = false; % disable limiter
% limiter = tvd();
% limiter = tvb(0.5);
limiter = tvb(1);

%% order test

t1 = 0.5;
nxlist = [10, 20, 40, 80, 160, 320];
n = numel(nxlist);
errors = zeros(3, n);

for w = 1:n
    [x, dx] = mesh_init_1d(xleft, xright, nxlist(w));

    init_func = @(s) alpha + beta * sin(s);
    uh0 = dg_projection(init_func, x, dx, pk, gk, basis);
    uh = dg_rk3_scheme(uh0, dx, t1, @burgers_f, @burgers_fhat_LF, @burgers_df, pk, gk, basis, basis_dx, limiter);

    u_exact_func = @(s) burgers_sin_exact(s, t1, alpha, beta);

    errors(:, w) = dg_errors(uh, u_exact_func, x, dx, pk, gk, basis);
end

show_results(nxlist, errors(1, :), errors(2, :), errors(3, :));

%% plot test

t2 = 1.5;
nxlist2 = [20, 80];

nx_ref = 320;
[x_ref, dx_ref] = mesh_init_1d(xleft, xright, nx_ref);

exact_ref_func = @(s) burgers_sin_exact(s, t2, alpha, beta);
% u_exact_ref = dg_projection(exact_ref_func, x_ref, dx_ref, pk, gk, basis);

figure();
hold('on');
plot(x_ref, exact_ref_func(x_ref), DisplayName = 'u-ref')

for w = 1:2
    [x, dx] = mesh_init_1d(xleft, xright, nxlist2(w));

    init_func = @(s) alpha + beta * sin(s);
    uh0 = dg_projection(init_func, x, dx, pk, gk, basis);
    uh = dg_rk3_scheme(uh0, dx, t2, @burgers_f, @burgers_fhat_LF, @burgers_df, pk, gk, basis, basis_dx, limiter);

    v = basis.eval(0, pk + 1); % column vector

    plot(x, v' * uh, DisplayName = sprintf('uh (n=%d)', nxlist2(w)))
end

hold('off');
box('on');
legend(Location = 'best');
exportgraphics(gcf, 'result_nodal.png', Resolution = 600);

function result = burgers_df(u)
    result = u;
end

function result = burgers_f(u)
    result = 0.5 * u .^ 2;
end

function result = burgers_get_alpha(ul, ur)
    result = max(abs(burgers_df(ul)), abs(burgers_df(ur)));
end

function result = burgers_fhat_LF(ul, ur)
    alpha = burgers_get_alpha(ul, ur);
    result = 0.25 * (ul .^ 2 + ur .^ 2) - 0.5 * alpha .* (ur - ul);
end
