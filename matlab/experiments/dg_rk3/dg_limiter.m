function lim = dg_limiter(raw_lim, pk, gk, h, basis)

    if pk == 0
        lim = @(u) u;
        return;
    end

    vl = basis.eval(-1, pk + 1); % column vector
    vr = basis.eval(1, pk + 1);

    [points, weights] = gauss_legendre(gk);

    u_values = basis.eval(points, pk + 1); % (gk, k+1)
    vmean = 0.5 * (weights' * u_values); % (1, k+1)

    if pk == 1
        M = [vmean(1), vmean(2); vl(1), vl(2)];
    else
        M = [vmean(1), vmean(2), vmean(3); vl(1), vl(2), vl(3); vr(1), vr(2), vr(3)];
    end

    lim = @(u) dg_limiter_inner(u, raw_lim, vl, vmean, vr, inv(M), pk, h);
end

function result = dg_limiter_inner(u, raw_lim, vl, vmean, vr, M_inv, pk, h)

    ul = vl' * u; % (1,n)
    um = vmean * u; % (1,n)
    ur = vr' * u; % (1,n)

    delte_r = circshift(um, -1) - um;
    delte_l = um - circshift(um, 1);

    ul_new = um - raw_lim(um - ul, delte_l, delte_r, h);
    ur_new = um + raw_lim(ur - um, delte_l, delte_r, h);

    result = zeros(size(u));

    if pk == 1
        B = [um; ul_new]; % (2,n)
        result = M_inv * B;
    else % pk >= 2
        B = [um; ul_new; ur_new]; % (3,n)
        result(1:3, :) = M_inv * B;
    end

end
