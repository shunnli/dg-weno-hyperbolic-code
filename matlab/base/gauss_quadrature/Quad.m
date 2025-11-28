classdef Quad
    % QUAD Gaussian quadrature integration.
    %
    % This class implements Gaussian quadrature integration methods.
    % It supports manual specification of quadrature nodes and weights,
    % or automatic generation via standard rules.
    %
    % QUAD Properties:
    %   nodes           - Quadrature nodes (column vector).
    %   weights         - Quadrature weights (column vector).
    %
    % QUAD Methods:
    %   integrate       - Numerically integrates a function over specified intervals.
    %
    % NOTE:
    %   Supported quadrature rules:
    %       - 'GaussLegendre': Gauss-Legendre quadrature.
    %       - 'GaussLobatto': Gauss-Lobatto quadrature.

    properties (SetAccess = private)
        nodes % Quadrature nodes (column vector)
        weights % Quadrature weights (column vector)
    end

    methods

        function obj = Quad(varargin)
            % Constructor for Quad class.
            %
            % Key-Value Parameters:
            %   'nodes'     - Quadrature nodes (numeric vector).
            %   'weights'   - Quadrature weights (numeric vector), same size as nodes.
            %   'type'      - Quadrature rule type ('GaussLegendre' or 'GaussLobatto'). (default: 'GaussLegendre')
            %   'gk'        - Number of quadrature points (integer, >= 2). (default: 5)
            %
            % USAGE:
            %   The constructor supports two modes:
            %       (1) Rule-based generation: Provide 'type' and 'gk'.
            %           obj = Quad(type = type, gk = gk)
            %
            %       (2) Manual specification: Provide both 'nodes' and 'weights'.
            %           obj = Quad(nodes = nodes, weights = weights)
            %
            %   Default: `obj = Quad()` is equivalent to `obj = Quad(type = 'Gauss-Legendre', gk = 5)`.
            %
            % EXAMPLE:
            %   q1 = Quad();
            %   q2 = Quad(type = 'Gauss-Legendre', gk = 5); % same as q1
            %   q3 = Quad(nodes = nodes, weights = weights);

            p = inputParser;
            addParameter(p, 'nodes', [], @(x) isnumeric(x) && isvector(x));
            addParameter(p, 'weights', [], @(x) isnumeric(x) && isvector(x));
            addParameter(p, 'type', 'gausslegendre', @(x) ischar(x));
            addParameter(p, 'gk', 5, @(x) isnumeric(x) && isscalar(x) && x >= 2);
            parse(p, varargin{:});
            args = p.Results;

            if ~isempty(args.nodes) && ~isempty(args.weights)
                % Manual nodes & weights mode

                assert(isequal(size(args.nodes), size(args.weights)), ...
                'nodes and weights must have the same size.');

                obj.nodes = args.nodes(:);
                obj.weights = args.weights(:);

            elseif ~isempty(args.type) && ~isempty(args.gk)
                % Rule-based mode

                type = lower(string(args.type));
                gk = args.gk;

                switch type
                    case 'gausslegendre'
                        [nodes, weights] = gauss_legendre(gk);
                    case 'gausslobatto'
                        [nodes, weights] = gauss_lobatto(gk);
                    otherwise
                        error('Quad:UnsupportedType', 'Unsupported quadrature type: %s', args.type);
                end

                obj.nodes = nodes(:);
                obj.weights = weights(:);

            else
                error('Quad:InvalidParameter', 'Invalid parameter combination. Provide either nodes+weights or type+gk.');
            end

        end

        function result = integrate(obj, f, xleft, xright)
            % INTEGRATE Numerically integrates a function over intervals.
            %
            % INPUT:
            %   f           - Function handle @(x) to integrate.
            %   xleft       - Left integration boundary.
            %   xright      - Right integration boundary. same size as xleft
            %
            % OUTPUT:
            %   result      - Numerical integral of f over each interval [xleft(i), xright(i)].
            %                 The result has the same size as xleft.
            %
            % Example:
            %   q = Quad('type', 'GaussLegendre', 'gk', 7);
            %   f = @(x) sin(x);
            %   xleft = [0, pi, 2*pi];
            %   xright = [pi, 2*pi, 3*pi];
            %   result = q.integrate(f, xleft, xright); % size = [1, 3]

            assert(isa(f, 'function_handle'), 'f must be a function handle.');
            assert(isnumeric(xleft) && isnumeric(xright), ...
            'xleft and xright must be numeric vectors.');
            assert(isequal(size(xleft), size(xright)), ...
            'xleft and xright must have the same size.');

            input_size = size(xleft);
            xleft = xleft(:)'; % size = [1, n]
            xright = xright(:)'; % size = [1, n]

            half_width = (xright - xleft) ./ 2; % size = [1, n]
            mid_point = (xleft + xright) ./ 2; % size = [1, n]

            nodes_local = mid_point + half_width .* obj.nodes; % size = [gk, n]
            integral_row = half_width .* sum(obj.weights .* f(nodes_local), 1); % size = [1, n]

            % reshape output
            result = reshape(integral_row, input_size);
        end

    end

end
