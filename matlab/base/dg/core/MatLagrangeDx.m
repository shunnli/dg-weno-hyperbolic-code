classdef MatLagrangeDx < MatBase
    % MatLagrangeDx: A class that provides the derivative of Lagrange polynomial basis functions based on a given set of points.
    %
    % Constructor:
    %   obj = MatLagrangeDx(x_points) - Generates the derivatives of the Lagrange polynomials based on the given points.
    %
    % EXAMPLE:
    %   points = [-1, 0, 1];
    %   basis_dx = MatLagrangeDx(points);
    %   x = linspace(-1, 1, 10);
    %   du = basis_dx.eval(x, numel(points)); % size = (3, 10)

    methods

        function obj = MatLagrangeDx(x_points)
            % Constructor for MatLagrangeDx
            %
            % Generates the derivatives of the Lagrange polynomials based on the given points.
            %
            % INPUT:
            %   x_points    - A vector of points where the Lagrange polynomials are constructed.
            %
            % EXAMPLE:
            %   basis_dx = MatLagrangeDx([-1, 0, 1]);

            if ~isvector(x_points) || numel(x_points) < 2
                error('MatLagrangeDx:InvalidInput', 'x_points must be a vector with at least two elements.');
            end

            % Generate Lagrange polynomial derivatives
            funcs = MatLagrangeDx.generateLagrangeFuncsDerivative(x_points);

            % Call parent constructor
            obj@MatBase(funcs);
        end

    end

    methods (Static)

        function funcs = generateLagrangeFuncsDerivative(x_points)
            % Generate the derivatives of the Lagrange polynomial functions based on the given points
            %
            % INPUT:
            %   x_points    - A vector of points for interpolation
            %
            % OUTPUT:
            %   funcs       - Cell array of Lagrange polynomial derivatives function handles

            n = numel(x_points);
            funcs = cell(1, n);

            % Loop over each point and construct the derivative of the corresponding Lagrange basis polynomial
            for i = 1:n
                % Derivative of Lagrange basis polynomial for point i

                L_prime = @(x) 0;

                for k = 1:n

                    if k == i
                        continue
                    end

                    term = @(x) 1;

                    for j = 1:n

                        if j == k || j == i
                            continue
                        end

                        term = @(x) term(x) .* (x - x_points(j));
                    end

                    L_prime = @(x) L_prime(x) + term(x);
                end

                c0 = 1;

                for j = 1:n

                    if j ~= i
                        c0 = c0 * (x_points(i) - x_points(j));
                    end

                end

                funcs{i} = @(x) L_prime(x) / c0; % Store derivative function handle
            end

        end

    end

end
