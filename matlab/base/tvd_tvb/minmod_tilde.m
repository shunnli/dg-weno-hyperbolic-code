function result = minmod_tilde(a1, a2, a3, tol)

    arguments
        a1 double
        a2 double
        a3 double
        tol (1, 1) double {mustBeNonnegative}
    end

    assert(isequal(size(a1), size(a2), size(a3)), 'a1, a2, and a3 must have the same size.');

    condition = abs(a1) < tol;
    result = zeros(size(a1));

    result(condition) = a1(condition);
    result(~condition) = minmod(a1(~condition), a2(~condition), a3(~condition));
end
