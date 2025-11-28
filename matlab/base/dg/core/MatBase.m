classdef MatBase < handle
    % MATBASE A base class for Discontinuous Galerkin (DG) basis functions.
    %
    % MatBase Properties:
    %   funcs      - Cell array of function handles representing the basis functions.
    %   funcs_num  - Number of basis functions.
    %
    % MatBase Methods:
    %   eval       - Evaluate the first 'n' basis functions at points 'x'.

    properties
        funcs % Basis function list (cell array of function handles)
        funcs_num % Number of basis functions
    end

    methods

        function obj = MatBase(funcs)
            % Constructor for MatBase
            %
            % INPUT:
            %   funcs   - A cell array of function handles representing basis functions.
            %
            % EXAMPLE:
            %   funcs = {@(x) x, @(x) x.^2, @(x) x.^3, @(x) sin(x)};
            %   basis = MatBase(funcs);

            if ~iscell(funcs) || isempty(funcs) || ~all(cellfun(@(f) isa(f, 'function_handle'), funcs))
                error('MatBase:InvalidInput', 'funcs must be a non-empty cell array of function handles.');
            end

            obj.funcs = funcs;
            obj.funcs_num = numel(funcs);
        end

        function u = eval(obj, x, n)
            % EVAL Evaluate the first 'n' basis functions at points 'x'.
            %
            % INPUT:
            %   x   - A vector of evaluation points in the range [-1, 1].
            %   n   - The number of basis functions to evaluate.
            %
            % OUTPUT:
            %   u   - A matrix containing function values.
            %           - If x is a scalar (1, 1), return a (n, 1) column vector.
            %           - If x is a row vector (1, m), return a (n, m) matrix.
            %           - If x is a column vector (m, 1), return a (m, n) matrix.
            %
            % EXAMPLE:
            %   funcs = {@(x) x, @(x) x.^2, @(x) x.^3};
            %   basis = MatBase(funcs);
            %   x = linspace(-1, 1, 5);
            %   u = basis.eval(x, 3);

            if ~isvector(x)
                error('MatBase:InvalidInput', 'x must be a vector.');
            end

            if n > obj.funcs_num
                error('MatBase:InvalidInput', 'n=%d exceeds available basis functions (%d).', n, obj.funcs_num);
            end

            % Preserve input shape
            is_row_vector = isrow(x);
            x_col = x(:); % Ensure x is a column vector

            % Evaluate basis functions
            u = zeros(numel(x_col), n);

            for i = 1:n
                u(:, i) = obj.funcs{i}(x_col);
            end

            % Restore original shape
            if is_row_vector
                u = u.';
            end

        end

    end

end
