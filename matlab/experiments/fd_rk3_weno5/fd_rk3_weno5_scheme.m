function u = fd_rk3_weno5_scheme(u, dx, tend, fl, fr, df)
    % fd_rk3_weno5_scheme
    %
    % INPUT:
    %   u         - Initial solution, must be a numeric vector or matrix.
    %   dx        - Spatial step size, must be a positive scalar.
    %   tend      - Final time, must be a positive scalar.
    %   fl        - Left flux function handle, defined as fl(u).
    %   fr        - Right flux function handle, defined as fr(u).
    %   df        - Derivative of the flux function, defined as df(u).
    %
    % OUTPUT:
    %   u         - Numerical solution at final time.

    assert(isnumeric(u) && (isvector(u) || ismatrix(u)), 'u must be a numeric vector or matrix.');
    assert(isnumeric(dx) && isscalar(dx) && dx > 0, 'dx must be a positive scalar.');
    assert(isnumeric(tend) && isscalar(tend) && tend > 0, 'tend must be a positive scalar.');
    assert(isa(fl, 'function_handle'), 'fl must be a function handle.');
    assert(isa(fr, 'function_handle'), 'fr must be a function handle.');
    assert(isa(df, 'function_handle'), 'df must be a function handle.');

    tnow = 0;

    while tnow < tend
        dt = 1 / (2 * max(abs(df(u)))) * dx ^ (5/3);
        dt = min([dt, tend - tnow]);

        u1 = u + dt * L_op(u, dx, fl, fr, df);
        u2 = (3/4) * u + (1/4) * (u1 + dt * L_op(u1, dx, fl, fr, df));
        u3 = (1/3) * u + (2/3) * (u2 + dt * L_op(u2, dx, fl, fr, df));
        u = u3;

        tnow = tnow + dt;
    end

end

function result = L_op(u, dx, fl, fr, df)
    c = max(abs(df(u)));

    hl = fl(u, c);
    hr = fr(u, c);

    [~, hlr] = weno5(hl);
    [hrl, ~] = weno5(hr);

    result = (-1) * (hlr + circshift(hrl, -1) - circshift(hlr, 1) - hrl) / dx;
end
