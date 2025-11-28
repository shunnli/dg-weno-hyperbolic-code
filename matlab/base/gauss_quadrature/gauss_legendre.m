function [x, w] = gauss_legendre(n)
    % GAUSS_LEGENDRE Computation of the Nodes and Weights for the Gauss-Legendre Quadrature.
    %
    % INPUT:
    %   n   - Number of nodes. (n >= 2)
    %
    % OUTPUT:
    %   x   - Nodes in [-1, 1], size(x) = [n, 1].
    %   w   - Weights, size(w) = [n, 1].
    %
    % EXAMPLE:
    %   [x, w] = gauss_legendre(3);
    %   % x = [0.7746; 0; -0.7746];
    %   % w = [0.5556; 0.8889; 0.5556];
    %
    % NOTE:
    %   x(i) from 1 to -1

    validateattributes(n, {'numeric'}, {'integer', 'nonnegative', '>=', 2});

    % [p_0(x), p_1(x), ..., p_n(x)]
    L = zeros(n, n + 1);

    % Initial guess.
    y = cos(pi * (2 * (1:n) - 1) / (2 * n))';
    delta = (0.27 / n) * sin(pi * linspace(-1, 1, n) * (n - 1) / (n + 1))';
    x = y + delta;

    % Newton's iteration.
    % Find the roots of p_n(x).
    while true
        L(:, 1) = 1; % p_0(x) = 1
        L(:, 2) = x; % p_1(x) = x

        for k = 2:n
            % p_k(x) = (2k-1)/k * x * p_{k-1}(x) - (k-1)/k * p_{k-2}(x)
            L(:, k + 1) = (2 * k - 1) / k * x .* L(:, k) - (k - 1) / k * L(:, k - 1);
        end

        % p_n'(x) = n/(x^2-1) * [x * p_n(x) - p_{n-1}(x)]
        dLn = (n ./ (x .^ 2 - 1)) .* (x .* L(:, n + 1) - L(:, n));
        dx = L(:, n + 1) ./ dLn;
        x = x - dx;

        if max(abs(dx)) < eps
            break
        end

    end

    % Compute the weights.
    w = 2 ./ ((1 - x .^ 2) .* dLn .^ 2);
end
