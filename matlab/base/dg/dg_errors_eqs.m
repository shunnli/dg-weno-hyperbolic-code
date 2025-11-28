function errors = dg_errors_eqs(uh, u, x, dx, pk, gk, basis, dim, trans)
    % DG_ERRORS_EQS Computes the L1, L2, and Linf errors for a system of equations.
    %
    % This function calculates the errors between the DG approximation `uh` and
    % the exact function `u` for a system of equations, using numerical quadrature.
    % The function includes a transformation `trans` before error computation.
    %
    % INPUT:
    %   uh     - DG coefficients of the numerical solution, size(uh) = [dim*(pk+1), nx].
    %   u      - Function handle representing the exact solution, returning `dim` outputs.
    %   x      - Array of cell center positions, must be a row vector, size(x) = [1, nx].
    %   dx     - Cell width (scalar, must be positive).
    %   pk     - Polynomial degree of the basis (non-negative integer).
    %   gk     - Number of quadrature points (positive integer).
    %   basis  - An instance of MatBase or its derived class.
    %   dim    - Number of equations in the system (positive integer).
    %   trans  - Function handle for transformation applied to solutions.
    %
    % OUTPUT:
    %   errors - A row vector containing:
    %            errors(1, 1) - L1 norm error.
    %            errors(2, 1) - L2 norm error.
    %            errors(3, 1) - Linf (maximum) norm error.
    %
    % NOTE:
    %   dim == nargout(trans) == nargin(trans) == nargout(u) (not checked by assert)

    assert(isnumeric(uh) && ismatrix(uh), 'uh must be a numeric matrix.');
    assert(isa(u, 'function_handle'), 'u must be a function handle.');
    assert(isrow(x) && isnumeric(x), 'x must be a numeric row vector.');
    assert(isscalar(dx) && dx > 0, 'dx must be a positive scalar.');
    assert(isscalar(pk) && pk >= 0 && mod(pk, 1) == 0, 'pk must be a non-negative integer.');
    assert(isscalar(gk) && gk > 0 && mod(gk, 1) == 0, 'gk must be a positive integer.');
    assert(isa(basis, 'MatBase'), 'basis must be an instance of MatBase or its derived class.');
    assert(isscalar(dim) && dim > 0 && mod(dim, 1) == 0, 'dim must be a positive integer.');
    assert(isa(trans, 'function_handle'), 'trans must be a function handle.');
    assert(2 * gk >= basis.funcs_num, ...
    'numerical quadrature formula requires at least twice the number of basis functions.');

    [points, weights] = gauss_legendre(gk);

    M = basis.eval(points, pk + 1);
    W = diag(weights);

    % x: row vector, points: column vector
    v_cs = cell(1, dim);
    [v_cs{:}] = u(x + dx / 2 * points); % v_cs{i}: matrix
    u_cs = cell(1, dim);
    [u_cs{:}] = trans(v_cs{:});

    vh_cs = cell(1, dim);

    for ii = 1:dim
        vh_cs{ii} = M * uh((ii - 1) * (pk + 1) + 1:ii * (pk + 1), :);
    end

    uh_cs = cell(1, dim);
    [uh_cs{:}] = trans(vh_cs{:});

    errors_l1 = zeros(1, dim);
    errors_l2 = zeros(1, dim);
    errors_inf = zeros(1, dim);

    for ii = 1:dim
        y = uh_cs{ii};
        y_exact = u_cs{ii};

        diff = y - y_exact;
        errors_inf(ii) = max(abs(diff(:)));

        diff1 = dx / 2 * W * diff;
        errors_l1(ii) = sum(abs(diff1(:)));

        diff2 = dx / 2 * W * (diff .^ 2);
        errors_l2(ii) = sum(diff2(:));
    end

    error_inf = max(errors_inf);
    error_l1 = sum(errors_l1);
    error_l2 = sqrt(sum(errors_l2));

    errors = [error_l1; error_l2; error_inf];
end
