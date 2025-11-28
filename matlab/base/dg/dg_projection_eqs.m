function u = dg_projection_eqs(f, x, dx, pk, gk, basis, dim)
    % DG_PROJECTION_EQS Computes the Discontinuous Galerkin (DG) projection for systems of equations.
    %
    % This function projects a system of equations onto a polynomial basis using
    % numerical quadrature in the Discontinuous Galerkin (DG) method. Unlike the
    % scalar function version, this function handles systems with multiple equations.
    %
    % INPUT:
    %   f      - Function handle representing the system of equations.
    %            It must return `dim` outputs, where each output is evaluated at given points.
    %   x      - Array of cell center positions, must be a row vector, size(x) = [1, nx].
    %   dx     - Cell width (scalar, must be positive).
    %   pk     - Polynomial degree of the basis (non-negative integer).
    %   gk     - Number of quadrature points (positive integer).
    %   basis  - An instance of MatBase or its derived class.
    %   dim    - Number of equations in the system (positive integer).
    %
    % OUTPUT:
    %   u      - Coefficients of the DG projection, size(u) = [(pk+1)*dim, nx].
    %
    % NOTE:
    %   dim == nargout(f) (not checked by assert)

    assert(isa(f, 'function_handle'), 'f must be a function handle.');
    assert(isrow(x) && isnumeric(x), 'x must be a numeric row vector.');
    assert(isscalar(dx) && dx > 0, 'dx must be a positive scalar.');
    assert(isscalar(pk) && pk >= 0 && mod(pk, 1) == 0, 'pk must be a non-negative integer.');
    assert(isscalar(gk) && gk > 0 && mod(gk, 1) == 0, 'gk must be a positive integer.');
    assert(isa(basis, 'MatBase'), 'basis must be an instance of MatBase or its derived class.');
    assert(2 * gk >= basis.funcs_num, ...
    'numerical quadrature formula requires at least twice the number of basis functions.');
    assert(isscalar(dim) && dim > 0 && mod(dim, 1) == 0, 'dim must be a positive integer.');

    [points, weights] = gauss_legendre(gk);

    M = basis.eval(points, pk + 1);
    W = diag(weights);

    Md = kron(eye(dim), M);
    Wd = kron(eye(dim), W);

    nx = numel(x);
    z = zeros(dim * gk, nx);

    v_cs = cell(dim, 1);

    for idx = 1:nx
        y = x(idx) + dx / 2 * points;
        [v_cs{:}] = f(y);
        z(:, idx) = cell2mat(v_cs);
    end

    u = (Md' * Wd * Md) \ (Md' * Wd * z);
end
