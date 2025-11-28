function f = eulereqs_f(v)
    assert(mod(size(v, 1), 3) == 0)

    gk = size(v, 1) / 3;

    v1 = v(1:gk, :); % size = [gk, n]
    v2 = v(gk + 1:2 * gk, :);
    v3 = v(2 * gk + 1:3 * gk, :);

    [rho, u, p] = eulereqs_trans2raw(v1, v2, v3); % size(rho) = [gk, n]

    f1 = rho .* u;
    f2 = rho .* u .^ 2 + p;
    f3 = u .* (v3 + p); % E = v3;

    f = [f1; f2; f3]; % size = [3*gk, n]
end
