function u = dg_rk3_scheme(u, dx, tend, f, fhat, df, pk, gk, basis, basis_dx, limiter)
    % dg_rk3_scheme
    %
    % INPUT:
    %   u         - Initial solution, must be a numeric vector or matrix.
    %   dx        - Spatial step size, must be a non-negative scalar.
    %   tend      - Final time, must be a non-negative scalar.
    %   f         - Flux function handle, defined as f(u).
    %   fhat      - Numerical flux function handle, defined as fhat(u_L, u_R, c).
    %   df        - Derivative of the flux function, defined as df(u).
    %   pk        - Polynomial order, must be a positive integer.
    %   gk        - Number of Gauss quadrature points, must be a positive integer.
    %   basis     - Basis function object, must be an instance of MatBase or its subclass.
    %   basis_dx  - Derivative of the basis function object, must be an instance of MatBase or its subclass.
    %   limiter   - Limiter function handle or false.
    %
    % OUTPUT:
    %   u         - Numerical solution at final time.

    assert(isnumeric(u) && (isvector(u) || ismatrix(u)), 'u must be a numeric vector or matrix.');
    assert(isscalar(dx) && dx >= 0, 'dx must be a non-negative scalar.');
    assert(isscalar(tend) && tend >= 0, 'tend must be a non-negative scalar.');
    assert(isa(f, 'function_handle'), 'f must be a function handle.');
    assert(isa(fhat, 'function_handle'), 'fhat must be a function handle.');
    assert(isa(df, 'function_handle'), 'df must be a function handle.');
    assert(isnumeric(pk) && isscalar(pk) && pk >= 0 && mod(pk, 1) == 0, 'pk must be a non-negative integer.');
    assert(isnumeric(gk) && isscalar(gk) && gk > 0 && mod(gk, 1) == 0, 'gk must be a positive integer.');
    assert(isa(basis, 'MatBase'), 'basis must be an object of class MatBase or its subclass.');
    assert(isa(basis_dx, 'MatBase'), 'basis_dx must be an object of class MatBase or its subclass.');
    assert(isa(limiter, 'function_handle') || isequal(limiter, false), 'limiter must be a function handle or false.');

    assert(2 * gk >= 2 * basis.funcs_num, 'insufficient precision of numerical quadrature formula');
    [points, weights] = gauss_legendre(gk);

    M = basis.eval(points, pk + 1);
    W = diag(weights);
    InvMass = inv(dx / 2 * M' * W * M);
    DM = 2 / dx * basis_dx.eval(points, pk + 1);

    vl = basis.eval(-1, pk + 1); % column vector
    vc = basis.eval(0, pk + 1);
    vr = basis.eval(1, pk + 1);

    params = struct();
    params.M = M;
    params.W = W;
    params.InvMass = InvMass;
    params.DM = DM;
    params.vl = vl;
    params.vc = vc;
    params.vr = vr;

    if isequal(limiter, false) || pk == 0
        post_precessor = @(u) u;
    else
        post_precessor = dg_limiter(limiter, pk, gk, dx, basis);
    end

    tnow = 0;

    while tnow < tend
        uh_mid = vc' * u;
        dt = 1 / ((2 * pk + 1) * max(abs(df(uh_mid)))) * dx ^ (max((pk + 1) / 3, 1));
        dt = min([dt, tend - tnow]);

        u1_pre = u + dt * L_op(u, dx, f, fhat, params);
        u1 = post_precessor(u1_pre);
        u2_pre = (3/4) * u + (1/4) * (u1 + dt * L_op(u1, dx, f, fhat, params));
        u2 = post_precessor(u2_pre);
        u3_pre = (1/3) * u + (2/3) * (u2 + dt * L_op(u2, dx, f, fhat, params));
        u3 = post_precessor(u3_pre);

        u = u3;

        tnow = tnow + dt;
    end

end

function result = L_op(u, dx, f, fhat, params)

    ul = params.vl' * u; % row vector
    ur = params.vr' * u;

    fhat_l = fhat(circshift(ur, 1), ul);
    fhat_r = fhat(ur, circshift(ul, -1));

    main = dx / 2 * (params.DM') * params.W * f(params.M * u);
    left_boundary = params.vl * fhat_l;
    right_boundary = params.vr * fhat_r;

    result = params.InvMass * (main - right_boundary + left_boundary);
end
