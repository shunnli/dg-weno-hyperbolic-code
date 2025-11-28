clc;
clear;
close all;

cd(fileparts(mfilename('fullpath')));
addpath(genpath('../../base'))
addpath('../common')

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

    uh0 = alpha + beta * sin(x);
    uh = fd_rk3_weno5_scheme(uh0, dx, t1, @burgers_fl, @burgers_fr, @burgers_df);

    u_exact = burgers_sin_exact(x, t1, alpha, beta);

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

u_exact_ref = burgers_sin_exact(x_ref, t2, alpha, beta);

figure;
hold on
plot(x_ref, u_exact_ref, DisplayName = 'u-ref')

for w = 1:2
    [x, dx] = mesh_init_1d(xleft, xright, nxlist2(w));

    uh0 = alpha + beta * sin(x);
    uh = fd_rk3_weno5_scheme(uh0, dx, t2, @burgers_fl, @burgers_fr, @burgers_df);

    plot(x, uh, DisplayName = sprintf('uh (n=%d)', nxlist2(w)))
end

hold off
legend('Location', 'best');

function result = burgers_fl(u, c)
    result = u .^ 2/4 + c / 2 .* u;
end

function result = burgers_fr(u, c)
    result = u .^ 2/4 - c / 2 .* u;
end

function result = burgers_df(u)
    result = u;
end
