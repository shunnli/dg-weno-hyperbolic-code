function result = order(xlist, ylist)
    % ORDER Computes the local order of convergence based on given data points.
    %
    % This function calculates the order of convergence for a given sequence
    % of x-values (xlist) and corresponding y-values (ylist). The order is
    % determined using the formula:
    %
    %   order = -log(y_{k+1} / y_k) / log(x_{k+1} / x_k)
    %
    % Usage:
    %   result = order(ylist)
    %   result = order(xlist, ylist)
    %
    % If only one input is provided, it's interpreted as ylist, and xlist is
    % set to [1, 2, 4, 8, ...].
    %
    % INPUT:
    %   xlist   - (optional) A numeric vector containing x-values.
    %   ylist   - A numeric vector containing y-values, must have the same size as xlist.
    %
    % OUTPUT:
    %   result  - A vector containing the computed order values.
    %             result is a row/column vector if xlist and ylist are row/column vectors.
    %
    % EXAMPLE:
    %   order([10, 20, 40], [0.5, 0.25, 0.125])   => [1, 1]
    %   order([10; 20; 40]', [0.5; 0.25; 0.125])  => [1; 1]
    %   order([0.5, 0.25, 0.125])  => [1, 1]
    %   order([0.5; 0.25; 0.125])  => [1; 1]

    if nargin == 1
        ylist = xlist;
        n = numel(ylist);
        xlist = reshape(2 .^ (0:n - 1), size(ylist)); % default to 1, 2, 4, ..., 2^(n-1)
    end

    assert(isvector(xlist) && isnumeric(xlist), 'xlist must be a numeric vector.');
    assert(isvector(ylist) && isnumeric(ylist), 'ylist must be a numeric vector.');
    assert(isequal(size(xlist), size(ylist)), 'xlist and ylist must have the same size.');
    assert(numel(xlist) >= 2, 'xlist and ylist must have at least two elements.');

    % Use column vector for computation
    is_row_vector = isrow(xlist);
    xlist = xlist(:);
    ylist = ylist(:);
    result = -log(ylist(2:end) ./ ylist(1:end - 1)) ./ log(xlist(2:end) ./ xlist(1:end - 1));

    % Match input orientation
    if is_row_vector
        result = result';
    end

end
