function errors = dg_errors(uh, u, x, dx, pk, gk, basis)
    % DG_ERRORS Computes the L1, L2, and Linf errors for the DG method.
    %
    % This function computes the errors between the DG approximation `uh` and
    % the exact function `u` using numerical quadrature.
    %
    % INPUT:
    %   uh     - DG coefficients of the numerical solution, size(uh) = [pk+1, nx].
    %   u      - Function handle representing the exact solution.
    %   x      - Array of cell center positions, must be a row vector, size(x) = [1, nx].
    %   dx     - Cell width (scalar, must be positive).
    %   pk     - Polynomial degree of the basis (non-negative integer).
    %   gk     - Number of quadrature points (positive integer).
    %   basis  - An instance of MatBase or its derived class.
    %
    % OUTPUT:
    %   errors - A row vector containing:
    %            errors(1, 1) - L1 norm error.
    %            errors(2, 1) - L2 norm error.
    %            errors(3, 1) - Linf (maximum) norm error.

    arguments
        uh double
        u (1, 1) function_handle
        x (1, :) double
        dx (1, 1) double {mustBePositive}
        pk (1, 1) double {mustBeInteger, mustBeNonnegative}
        gk (1, 1) double {mustBeInteger, mustBePositive}
        basis (1, 1) MatBase
    end

    assert(2 * gk >= basis.funcs_num, 'Numerical quadrature formula requires at least twice the number of basis functions.');

    [points, weights] = gauss_legendre(gk);

    M = basis.eval(points, pk + 1);
    W = diag(weights);

    y = M * uh;
    y_exact = u(x + dx / 2 * points);
    % x is a row vector, points is a column vector

    diff = y - y_exact;
    error_inf = max(abs(diff(:)));

    diff1 = dx / 2 * W * diff;
    error_l1 = sum(abs(diff1(:)));

    diff2 = dx / 2 * W * (diff .^ 2);
    error_l2 = sqrt(sum(diff2(:)));

    errors = [error_l1; error_l2; error_inf];
end
