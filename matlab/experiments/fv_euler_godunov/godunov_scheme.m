function u = godunov_scheme(u, dx, tend, fhat, df)
    % godunov_scheme
    %
    % INPUT:
    %   u         - Initial solution, must be a numeric vector or matrix.
    %   dx        - Spatial step size, must be a non-negative scalar.
    %   tend      - Final time, must be a non-negative scalar.
    %   fhat      - Numerical flux function handle, defined as fhat(u_L, u_R, c).
    %   df        - Derivative of the flux function, defined as df(u).
    %
    % OUTPUT:
    %   u         - Numerical solution at final time.

    assert(isnumeric(u) && (isvector(u) || ismatrix(u)), 'u must be a numeric vector or matrix.');
    assert(isscalar(dx) && dx > 0, 'dx must be a positive scalar.');
    assert(isscalar(tend) && tend >= 0, 'tend must be a non-negative scalar.');
    assert(isa(fhat, 'function_handle'), 'fhat must be a function handle.');
    assert(isa(df, 'function_handle'), 'df must be a function handle.');

    tnow = 0;

    while tnow < tend
        dt = dx / (2 * max(abs(df(u))));
        dt = min([dt, tend - tnow]);

        u = u + dt * L_op(u, dx, fhat);
        tnow = tnow + dt;
    end

end

function result = L_op(u, dx, fhat)
    fhat_left = fhat(circshift(u, 1), u);
    fhat_right = fhat(u, circshift(u, -1));
    result = (-1) * (fhat_right - fhat_left) / dx;
end
