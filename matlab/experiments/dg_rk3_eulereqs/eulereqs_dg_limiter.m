function lim = eulereqs_dg_limiter(raw_lim, pk, gk, h, basis, is_riemann)

    if pk == 0
        lim = @(u) u;
        return;
    end

    vl = basis.eval(-1, pk + 1); % column vector, (k+1,1)
    vr = basis.eval(1, pk + 1);

    [points, weights] = gauss_legendre(gk);

    u_values = basis.eval(points, pk + 1); % (gk, k+1)
    vm = 0.5 * (weights' * u_values); % (1, k+1) !

    if pk == 1
        M = [vm(1), vm(2); vl(1), vl(2)];
    else
        M = [vm(1), vm(2), vm(3); vl(1), vl(2), vl(3); vr(1), vr(2), vr(3)];
    end

    lim = @(u) dg_limiter_inner(u, raw_lim, vl, vm, vr, inv(M), pk, h, is_riemann);
end

function result = dg_limiter_inner(u, raw_lim, vl, vm, vr, M_inv, pk, h, is_riemann)

    u1l = vl' * u(1:(pk + 1), :); % (1, n)
    u2l = vl' * u(pk + 2:2 * (pk + 1), :);
    u3l = vl' * u(2 * (pk + 1) + 1:3 * (pk + 1), :);

    u1r = vr' * u(1:(pk + 1), :);
    u2r = vr' * u(pk + 2:2 * (pk + 1), :);
    u3r = vr' * u(2 * (pk + 1) + 1:3 * (pk + 1), :);

    u1m = vm * u(1:(pk + 1), :);
    u2m = vm * u(pk + 2:2 * (pk + 1), :);
    u3m = vm * u(2 * (pk + 1) + 1:3 * (pk + 1), :);

    ul = [u1l; u2l; u3l]; % (3,n)
    ur = [u1r; u2r; u3r];
    um = [u1m; u2m; u3m];

    n = size(u, 2);
    ul_new = zeros(3, n);
    ur_new = zeros(3, n);
    um_new = zeros(3, n);

    for i = 1:n
        idxl = mod(i - 2, n) + 1;
        idxr = mod(i, n) + 1;

        % A(U) = R*D*L
        % U -> W = L U
        % W_t + D W_x = 0
        [L, ~, R] = eulereqs_DF(um(1, i), um(2, i), um(3, i));

        wm = L * um(:, i); % (3,1)
        wlm = L * um(:, idxl);
        wrm = L * um(:, idxr);

        wl = L * ul(:, i);
        wr = L * ur(:, i);

        wl_new = wm - raw_lim(wm - wl, wm - wlm, wrm - wm, h); % (3,1)
        wr_new = wm + raw_lim(wr - wm, wm - wlm, wrm - wm, h);

        ul_new(:, i) = R * wl_new; % (3,1)
        ur_new(:, i) = R * wr_new;
        um_new(:, i) = R * wm;
    end

    result = zeros(size(u)); % (3(pk+1), n)

    if pk == 1
        B1 = [um_new(1, :); ul_new(1, :)]; % (2,n)
        result(1:2, :) = M_inv * B1;

        B2 = [um_new(2, :); ul_new(2, :)];
        result(3:4, :) = M_inv * B2;

        B3 = [um_new(3, :); ul_new(3, :)];
        result(5:6, :) = M_inv * B3;
    else % pk >= 2
        B1 = [um_new(1, :); ul_new(1, :); ur_new(1, :)]; % (3,n)
        result(1:3, :) = M_inv * B1;

        B2 = [um_new(2, :); ul_new(2, :); ur_new(2, :)];
        result((pk + 1) + 1:(pk + 1) + 3, :) = M_inv * B2;

        B3 = [um_new(3, :); ul_new(3, :); ur_new(3, :)];
        result(2 * (pk + 1) + 1:2 * (pk + 1) + 3, :) = M_inv * B3;
    end

    if is_riemann
        result(:, 1) = u(:, 1);
        result(end, 1) = u(end, 1);
    end

end
