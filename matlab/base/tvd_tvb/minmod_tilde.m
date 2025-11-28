function result = minmod_tilde(a1, a2, a3, tol)
    assert(isnumeric(a1) && isnumeric(a2) && isnumeric(a3), 'a1, a2, and a3 must be numeric arrays.');
    assert(isequal(size(a1), size(a2), size(a3)), 'a1, a2, and a3 must have the same size.');
    assert(isnumeric(tol) && isscalar(tol) && tol >= 0, 'm must be a non-negative scalar.');

    condition = abs(a1) < tol;
    result = zeros(size(a1));

    result(condition) = a1(condition);
    result(~condition) = minmod(a1(~condition), a2(~condition), a3(~condition));
end
