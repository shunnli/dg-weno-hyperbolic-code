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
    %   dim == nargout(f) (not checked)

    arguments
        f (1, 1) function_handle
        x (1, :) double
        dx (1, 1) double {mustBePositive}
        pk (1, 1) double {mustBeInteger, mustBeNonnegative}
        gk (1, 1) double {mustBeInteger, mustBePositive}
        basis (1, 1) MatBase
        dim (1, 1) double {mustBeInteger, mustBePositive}
    end

    assert(2 * gk >= basis.funcs_num, 'numerical quadrature formula requires at least twice the number of basis functions.');

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
