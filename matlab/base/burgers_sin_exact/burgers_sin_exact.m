function u = burgers_sin_exact(x, t, alpha, beta)
    % BURGERS_SIN_EXACT Computes the exact solution of Burgers' equation (u_t + u * u_x = 0)
    %       with 2pi periodic boundary condition at time t.
    %
    %       u(x, 0) = alpha + beta * sin(x), x in [-pi, pi]
    %
    % INPUT:
    %   x       - Spatial coordinate (matrix, vector, or scalar), numeric
    %   t       - Time (scalar), numeric and non-negative
    %   alpha   - Scalar, numeric
    %   beta    - Scalar, numeric
    %
    % OUTPUT:
    %   u       - Solution of Burgers' equation at (x, t). size(u) == size(x)

    assert(isnumeric(x), 'x must be a numeric array.');
    assert(isnumeric(t) && isscalar(t) && t >= 0, 't must be a non-negative scalar.');
    assert(isnumeric(alpha) && isscalar(alpha), 'alpha must be a scalar numeric value.');
    assert(isnumeric(beta) && isscalar(beta), 'beta must be a scalar numeric value.');

    u = alpha + beta * burgers_sin_newton(x - alpha * t, beta * t);
end

function u = burgers_sin_newton(x, t)
    % u_t + u u_x = 0, x in [-pi, pi]
    % u(x, 0) = sin(x)

    x = mod(x + pi, 2 * pi) - pi;

    % For t >= Tb, the algorithm becomes highly sensitive
    % to the choice of the initial value for iteration.
    u = x / (pi / 2 + t);

    % Newton iteration
    for iter = 1:1e5
        du = (u - sin(x - u * t)) ./ (1 + cos(x - u * t) * t);
        u = u - du;

        if max(abs(du)) < 1e-10
            break
        end

    end

end
