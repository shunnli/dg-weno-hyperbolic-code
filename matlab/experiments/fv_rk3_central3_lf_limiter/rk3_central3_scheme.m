function u = rk3_central3_scheme(u, dx, tend, fhat, df, limiter)
    % rk3_central3_scheme
    %
    % INPUT:
    %   u         - Initial solution, must be a numeric vector or matrix.
    %   dx        - Spatial step size, must be a positive scalar.
    %   tend      - Final time, must be a positive scalar.
    %   fhat      - Numerical flux function handle, defined as fhat(u_L, u_R, c).
    %   df        - Derivative of the flux function, defined as df(u).
    %   limiter   - Limiter function handle or false.
    %
    % OUTPUT:
    %   u         - Numerical solution at final time.

    assert(isnumeric(u) && (isvector(u) || ismatrix(u)), 'u must be a numeric vector or matrix.');
    assert(isnumeric(dx) && isscalar(dx) && dx > 0, 'dx must be a positive scalar.');
    assert(isnumeric(tend) && isscalar(tend) && tend > 0, 'tend must be a positive scalar.');
    assert(isa(fhat, 'function_handle'), 'fhat must be a function handle.');
    assert(isa(df, 'function_handle'), 'df must be a function handle.');
    assert(isa(limiter, 'function_handle') || isequal(limiter, false), 'limiter must be a function handle or false.');

    if isequal(limiter, false)
        precessor = @(ul_plus, ur_minus, ul, u, ur, dx) deal(ul_plus, ur_minus);
    else
        precessor = fv_limiter(limiter);
    end

    tnow = 0;

    while tnow < tend
        dt = dx / (2 * max(abs(df(u))));
        dt = min([dt, tend - tnow]);

        u1 = u + dt * L_op(u, dx, fhat, precessor);
        u2 = (3/4) * u + (1/4) * (u1 + dt * L_op(u1, dx, fhat, precessor));
        u3 = (1/3) * u + (2/3) * (u2 + dt * L_op(u2, dx, fhat, precessor));
        u = u3;

        tnow = tnow + dt;
    end

end

function result = L_op(u, dx, fhat, precessor)
    ur = circshift(u, -1);
    ul = circshift(u, 1);

    ul_plus = (1/3 * ul) + (5/6 * u) - (1/6 * ur);
    ur_minus = (-1/6 * ul) + (5/6 * u) + (1/3 * ur);

    [ul_plus, ur_minus] = precessor(ul_plus, ur_minus, ul, u, ur, dx);

    ul_plus_right = circshift(ul_plus, -1);
    ur_minus_left = circshift(ur_minus, 1);

    fhat_left = fhat(ur_minus_left, ul_plus);
    fhat_right = fhat(ur_minus, ul_plus_right);

    result = (-1) * (fhat_right - fhat_left) / dx;
end
