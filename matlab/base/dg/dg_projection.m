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

    assert(isa(f, 'function_handle'), 'f must be a function handle.');
    assert(isrow(x) && isnumeric(x), 'x must be a numeric row vector.');
    assert(isscalar(dx) && dx > 0, 'dx must be a positive scalar.');
    assert(isscalar(pk) && pk >= 0 && mod(pk, 1) == 0, 'pk must be a non-negative integer.');
    assert(isscalar(gk) && gk > 0 && mod(gk, 1) == 0, 'gk must be a positive integer.');
    assert(isa(basis, 'MatBase'), 'basis must be an instance of MatBase or its derived class.');
    assert(2 * gk >= basis.funcs_num, ...
    'numerical quadrature formula requires at least twice the number of basis functions.');

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
