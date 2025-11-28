classdef MatLagrange < MatBase
    % MatLagrange: A class that provides Lagrange polynomial basis functions based on a given set of points.
    %
    % Constructor:
    %   obj = MatLagrange(x_points) - Generates the Lagrange polynomials based on the given points.
    %
    % EXAMPLE:
    %   points = [-1, 0, 1];
    %   basis = MatLagrange(points);
    %   x = linspace(-1, 1, 10);
    %   u = basis.eval(x, numel(points));

    methods

        function obj = MatLagrange(x_points)
            % Constructor for MatLagrange
            %
            % Generates the Lagrange polynomials based on the given points.
            %
            % INPUT:
            %   x_points - A vector of points where the Lagrange polynomials are constructed.
            %
            % EXAMPLE:
            %   basis = MatLagrange([-1, 0, 1]);

            if ~isvector(x_points) || numel(x_points) < 2
                error('MatLagrange:InvalidInput', 'x_points must be a vector with at least two elements.');
            end

            % Generate Lagrange polynomial function handles
            funcs = MatLagrange.generateLagrangeFuncs(x_points);

            % Call parent constructor
            obj@MatBase(funcs);
        end

    end

    methods (Static)

        function funcs = generateLagrangeFuncs(x_points)
            % Generate the Lagrange polynomial functions based on the given points
            %
            % INPUT:
            %   x_points    - A vector of points for interpolation
            %
            % OUTPUT:
            %   funcs       - Cell array of Lagrange polynomial function handles

            n = numel(x_points);
            funcs = cell(1, n);

            % Loop over each point and construct the corresponding Lagrange basis polynomial
            for i = 1:n
                % Lagrange basis polynomial for point i
                L = @(x) 1;

                for j = 1:n

                    if i ~= j
                        L = @(x) L(x) .* (x - x_points(j)) / (x_points(i) - x_points(j));
                    end

                end

                funcs{i} = L; % Store function handle in cell array
            end

        end

    end

end
