function [res_ul, res_ur] = weno5(u)
    % WENO-JS 5

    % Domain cells reference:
    %
    %                |           |   u(i)    |           |
    %                |  u(i-1)   |___________|           |
    %                |___________|           |   u(i+1)  |
    %                |           |           |___________|
    %             ...|-----0-----|-----0-----|-----0-----|...
    %                |    i-1    |     i     |    i+1    |
    %                |-         +|-         +|-         +|
    %              i-3/2       i-1/2       i+1/2       i+3/2
    %
    % WENO stencils reference:
    %
    %
    %                               |___________S2__________|
    %                               |                       |
    %                       |___________S1__________|       |
    %                       |                       |       |
    %               |___________S0__________|       |       |
    %             ..|---o---|---o---|---o---|---o---|---o---|...
    %               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
    %                                      -|
    %                                     i+1/2
    %
    %
    %               |___________S0__________|
    %               |                       |
    %               |       |___________S1__________|
    %               |       |                       |
    %               |       |       |___________S2__________|
    %             ..|---o---|---o---|---o---|---o---|---o---|...
    %               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
    %                               |+
    %                             i-1/2
    %
    % WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]

    assert(size(u, 1) == 1);

    % linear weights
    d_l = [3.0/10; 3.0/5; 1.0/10];
    d_r = [1.0/10; 3.0/5; 3.0/10];

    weno_ep = 1e-6;

    % Note: by using circshift over our domain, we are implicitly creating
    % favorable code that includes periodical boundary conditions.
    ul_mat = [circshift(u, 2); circshift(u, 1); u];
    u_mat = [circshift(u, 1); u; circshift(u, -1)];
    ur_mat = [u; circshift(u, -1); circshift(u, -2)];

    % smooth indicator
    b = zeros(3, numel(u));
    b(1, :) = 13.0/12 * ([1, -2, 1] * ul_mat) .^ 2 + 1.0/4 * ([1, -4, 3] * ul_mat) .^ 2;
    b(2, :) = 13.0/12 * ([1, -2, 1] * u_mat) .^ 2 + 1.0/4 * ([1, 0, -1] * u_mat) .^ 2;
    b(3, :) = 13.0/12 * ([1, -2, 1] * ur_mat) .^ 2 + 1.0/4 * ([3, -4, 1] * ur_mat) .^ 2;

    % nonlinear weights
    b_ep = (weno_ep + b) .^ 2;
    b_ep_inv = 1 ./ b_ep;

    a_l = d_l .* b_ep_inv;
    a_r = d_r .* b_ep_inv;

    % normalize
    norm_factor_l = sum(a_l, 1);
    norm_factor_r = sum(a_r, 1);

    w_l = a_l ./ norm_factor_l;
    w_r = a_r ./ norm_factor_r;

    coeff_l = [-1/6, 5/6, 1/3;
               1/3, 5/6, -1/6;
               11/6, -7/6, 1/3];

    coeff_r = [1/3, -7/6, 11/6;
               -1/6, 5/6, 1/3;
               1/3, 5/6, -1/6];

    % reconstruct
    u_l = [coeff_l(1, :) * ul_mat;
           coeff_l(2, :) * u_mat;
           coeff_l(3, :) * ur_mat];
    u_r = [coeff_r(1, :) * ul_mat;
           coeff_r(2, :) * u_mat;
           coeff_r(3, :) * ur_mat];

    res_ul = sum(w_l .* u_l, 1);
    res_ur = sum(w_r .* u_r, 1);
end
