function show_results(nums, error_l1, error_l2, error_inf)
    % SHOW_RESULTS Displays a table of errors and convergence orders.
    %
    % This function prints a LaTeX-style table with error values and
    % corresponding convergence orders for L1, L2, and Linf norms.
    %
    % INPUT:
    %   nums      - A numeric vector containing problem sizes.
    %   error_l1  - A numeric vector containing L1 norm errors.
    %   error_l2  - A numeric vector containing L2 norm errors.
    %   error_inf - A numeric vector containing Linf norm errors.
    %
    %   See also ORDER.

    arguments
        nums double {mustBeVector, mustBePositive}
        error_l1 double {mustBeVector, mustBePositive}
        error_l2 double {mustBeVector, mustBePositive}
        error_inf double {mustBeVector, mustBePositive}
    end

    assert(isequal(size(nums), size(error_l1), size(error_l2), size(error_inf)), 'All input vectors must have the same size.');

    order_l1 = order(nums, error_l1);
    order_l2 = order(nums, error_l2);
    order_inf = order(nums, error_inf);

    fprintf('   n & error_l1 & order & error_l2 & order & error_linf & order \\\\ \n');
    fprintf('%4d & %.2e &   -  & %.2e &   -  & %.2e &   -  \\\\ \n', ...
        nums(1), error_l1(1), error_l2(1), error_inf(1));

    for cnt = 2:numel(nums)
        fprintf('%4d & %.2e & %.2f & %.2e & %.2f & %.2e & %.2f \\\\ \n', ...
            nums(cnt), error_l1(cnt), order_l1(cnt - 1), ...
            error_l2(cnt), order_l2(cnt - 1), error_inf(cnt), order_inf(cnt - 1));
    end

end
