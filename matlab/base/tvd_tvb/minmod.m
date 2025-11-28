function result = minmod(a1, a2, a3, ~)
    assert(isnumeric(a1) && isnumeric(a2) && isnumeric(a3), 'a1, a2, and a3 must be numeric arrays.');
    assert(isequal(size(a1), size(a2), size(a3)), 'a1, a2, and a3 must have the same size.');

    pos_mask = (a1 > 0) & (a2 > 0) & (a3 > 0);
    neg_mask = (a1 < 0) & (a2 < 0) & (a3 < 0);

    result = zeros(size(a1));
    result(pos_mask) = min(min(a1(pos_mask), a2(pos_mask)), a3(pos_mask));
    result(neg_mask) = max(max(a1(neg_mask), a2(neg_mask)), a3(neg_mask));
end
