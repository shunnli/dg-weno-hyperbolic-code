classdef MatLegendreDx < MatBase
    % MatLegendreDx: A class that provides Legendre polynomial derivatives.
    %
    % Constructor:
    %   obj = MatLegendreDx(n) - Generates the first 'n' Legendre polynomial derivatives.
    %
    % EXAMPLE:
    %   basis_dx = MatLegendreDx(5);
    %   x = linspace(-1, 1, 10);
    %   u = basis_dx.eval(x, 5); % size = (5, 10)
    %
    % NOTE:
    %   p0(x) = 1
    %   p1(x) = x
    %   p2(x) = (1/2)*(3*x^2 - 1)
    %   p3(x) = (1/2)*(5*x^3 - 3*x)
    %   p4(x) = (1/8)*(35*x^4 - 30*x^2 + 3)

    methods

        function obj = MatLegendreDx(n)
            % Constructor for MatLegendreDx
            %
            % Generates the first 'n' derivatives of Legendre polynomials.
            %
            % INPUT:
            %   n       - Number of Legendre polynomials to generate.
            %
            % EXAMPLE:
            %   basis_dx = MatLegendreDx(5);

            if ~isscalar(n) || n < 1 || mod(n, 1) ~= 0
                error('MatLegendreDx:InvalidInput', 'n must be a positive integer.');
            end

            % Generate the derivatives of Legendre polynomial function handles
            dfuncs = MatLegendreDx.generateLegendreDerivatives(n);

            % Call parent constructor
            obj@MatBase(dfuncs);
        end

    end

    methods (Static)

        function dfuncs = generateLegendreDerivatives(n)
            % Generate the first 'n' Legendre polynomial derivatives
            %
            % INPUT:
            %   n       - Number of polynomials to generate
            %
            % OUTPUT:
            %   dfuncs  - Cell array of derivative function handles
            %
            % NOTE:
            %   Return NaN if x^2 == 1 because 1/(x^2-1) is used in the computation of the derivative.

            dfuncs = cell(1, n);
            dfuncs{1} = @(x) zeros(size(x));

            for k = 2:n
                dfuncs{k} = @(x) (k - 1) ./ (x .^ 2 - 1) .* (x .* legendreP(k - 1, x) - legendreP(k - 2, x));
            end

        end

    end

end
