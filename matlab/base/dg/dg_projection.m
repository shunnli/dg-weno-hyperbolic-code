function u = dg_projection(f, x, dx, pk, gk, basis)
    % DG_PROJECTION Computes the Discontinuous Galerkin (DG) projection of a function.
    %
    % This function projects a given function f onto a polynomial basis using
    % numerical quadrature in the Discontinuous Galerkin (DG) method.
    %
    % INPUT:
    %   f      - Function handle representing the function to be projected.
    %   x      - Array of cell center positions, must be a row vector, size(x) = [1, nx].
    %   dx     - Cell width (scalar, must be positive).
    %   pk     - Polynomial degree of the basis (non-negative integer).
    %   gk     - Number of quadrature points (positive integer).
    %   basis  - An instance of MatBase or its derived class.
    %
    % OUTPUT:
    %   u      - Coefficients of the DG projection, size(u) = [pk+1, nx].

    arguments
        f (1, 1) function_handle
        x (1, :) double
        dx (1, 1) double {mustBePositive}
        pk (1, 1) double {mustBeInteger, mustBeNonnegative}
        gk (1, 1) double {mustBeInteger, mustBePositive}
        basis (1, 1) MatBase
    end

    assert(2 * gk >= basis.funcs_num, 'numerical quadrature formula requires at least twice the number of basis functions.');

    [points, weights] = gauss_legendre(gk);

    M = basis.eval(points, pk + 1);
    W = diag(weights);

    nx = numel(x);
    z = zeros(gk, nx);

    for idx = 1:nx
        y = x(idx) + dx / 2 * points;
        z(:, idx) = f(y);
    end

    u = (M' * W * M) \ (M' * W * z);
end
