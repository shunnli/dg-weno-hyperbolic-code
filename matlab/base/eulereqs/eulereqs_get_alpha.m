function alpha = eulereqs_get_alpha(uh)
    assert(size(uh, 1) == 3);

    [rho, u, p] = eulereqs_trans2raw(uh(1, :), uh(2, :), uh(3, :)); % row vector
    gamma = 1.4;
    c = sqrt(gamma * abs(p ./ rho)); % rho may be negative
    alpha = max([abs(u - c); abs(u); abs(u + c)]); % max on each column
end
