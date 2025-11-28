function tests = test_gauss_rules()
    % TESTGAUSSRULES Unit tests for Gauss quadrature rules (gauss_legendre, gauss_lobatto)
    tests = functiontests(localfunctions);
end

%% === Test Functions ===

function testGaussLegendre(testCase)
    nRange = 2:15;

    for n = nRange
        [x, w] = gauss_legendre(n);
        [deg_min, deg_max] = deal(1, 2 * n - 1);

        for k = deg_min:(deg_max + 1)
            f = @(x) x .^ k;
            exact = (mod(k, 2) == 0) * 2 / (k + 1);
            approx = sum(f(x) .* w);
            err = abs(approx - exact);

            if k <= deg_max
                verifyLessThan(testCase, err, 10 * eps, ...
                    sprintf('Gauss-Legendre failed for n=%d, k=%d', n, k));
            end

        end

    end

end

function testGaussLobatto(testCase)
    nRange = 2:15;

    for n = nRange
        [x, w] = gauss_lobatto(n);
        [deg_min, deg_max] = deal(2, 2 * n - 3);

        for k = deg_min:(deg_max + 1)
            f = @(x) x .^ k;
            exact = (mod(k, 2) == 0) * 2 / (k + 1);
            approx = sum(f(x) .* w);
            err = abs(approx - exact);

            if k <= deg_max
                verifyLessThan(testCase, err, 10 * eps, ...
                    sprintf('Gauss-Lobatto failed for n=%d, k=%d', n, k));
            end

        end

    end

end
