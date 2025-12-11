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

%% order test

t1 = 0.5;
nxlist = [10, 20, 40, 80, 160, 320, 640];
n = numel(nxlist);
errors = zeros(3, n);

for w = 1:n
    [x, dx] = mesh_init_1d(xleft, xright, nxlist(w));

    init_func = @(s) alpha + beta * sin(s);
    uh0 = Quad().integrate(init_func, x - dx / 2, x + dx / 2) / dx;
    uh = godunov_scheme(uh0, dx, t1, @burgers_fhat_godunov, @burgers_df);

    exact_func = @(s) burgers_sin_exact(s, t1, alpha, beta);
    u_exact = Quad().integrate(exact_func, x - dx / 2, x + dx / 2) / dx;

    errors(1, w) = sum(abs(uh - u_exact)) * dx;
    errors(2, w) = sqrt(sum(abs(uh - u_exact) .^ 2) * dx);
    errors(3, w) = max(abs(uh - u_exact));
end

show_results(nxlist, errors(1, :), errors(2, :), errors(3, :));

%% plot test

t2 = 1.5;
nxlist2 = [20, 80];
nx_ref = 320;
[x_ref, dx_ref] = mesh_init_1d(xleft, xright, nx_ref);

exact_func = @(s) burgers_sin_exact(s, t2, alpha, beta);
u_exact_ref = Quad().integrate(exact_func, x_ref - dx_ref / 2, x_ref + dx_ref / 2) / dx_ref;

figure();
hold('on');
plot(x_ref, u_exact_ref, DisplayName = 'u-ref')

for w = 1:2
    [x, dx] = mesh_init_1d(xleft, xright, nxlist2(w));

    init_func = @(s) alpha + beta * sin(s);
    uh0 = Quad().integrate(init_func, x - dx / 2, x + dx / 2) / dx;
    uh = godunov_scheme(uh0, dx, t2, @burgers_fhat_godunov, @burgers_df);

    plot(x, uh, DisplayName = sprintf('uh (n=%d)', nxlist2(w)))
end

hold('off');
box('on');
legend(Location = 'best');
exportgraphics(gcf, 'result.png', Resolution = 600);

function result = burgers_fhat_godunov(ul, ur, ~)
    result = zeros(size(ul));

    % ul <= ur
    % result = 0 if ul < 0 < ur, otherwise min(ul^2/2, ur^2/2)
    condition = (ul <= ur);
    result(condition) = min(ul(condition) .^ 2/2, ur(condition) .^ 2/2) .* (ul(condition) .* ur(condition) > 0);

    % ul >= ur
    % max(ul^2/2, ur^2/2)
    result(~condition) = max(ul(~condition) .^ 2/2, ur(~condition) .^ 2/2);
end

function result = burgers_df(u)
    result = u;
end
