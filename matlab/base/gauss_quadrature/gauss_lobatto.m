function [x, w] = gauss_lobatto(n)
    % GAUSS_LOBATTO Computation of the Nodes and Weights for the Gauss-Lobatto Quadrature.
    %
    % INPUT:
    %   n   - Number of nodes. (n >= 2)
    %
    % OUTPUT:
    %   x   - Nodes in [-1, 1], size(x) = [n, 1].
    %   w   - Weights, size(w) = [n, 1].
    %
    % EXAMPLE:
    %   [x, w] = gauss_lobatto(3);
    %   % x = [1; 0; -1];
    %   % w = [0.3333; 1.3333; 0.3333];
    %
    % NOTE:
    %   x(i) from 1 to -1

    validateattributes(n, {'numeric'}, {'integer', 'nonnegative', '>=', 2});

    % [p_0(x), p_1(x), ..., p_{n-1}(x)]
    L = zeros(n, n);

    % Initial guess.
    x = cos(pi * (0:n - 1) / (n - 1))';

    % Newton's iteration.
    % Find the roots of q_n(x) = (1-x^2) p_{n-1}'(x).
    while true
        L(:, 1) = 1; % p_0(x) = 1
        L(:, 2) = x; % p_1(x) = x

        for k = 2:n - 1
            % p_k(x) = (2k-1)/k * x * p_{k-1}(x) - (k-1)/k * p_{k-2}(x)
            L(:, k + 1) = (2 * k - 1) / k * x .* L(:, k) - (k - 1) / k * L(:, k - 1);
        end

        % q_n(x) = (1-x^2) p_{n-1}'(x) = -(n-1) * [x * p_{n-1}(x) - p_{n-2}(x)]
        % q_n'(x) = -(n-1)n * p_{n-1}(x)
        dx = (x .* L(:, n) - L(:, n - 1)) ./ (n * L(:, n));
        x = x - dx;

        if max(abs(dx)) < eps
            break
        end

    end

    % Compute the weights.
    w = 2 ./ ((n - 1) * n * L(:, n) .^ 2);
end
