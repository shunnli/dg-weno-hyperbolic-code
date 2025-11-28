function limiter = tvb(m)
    assert(isnumeric(m) && isscalar(m) && m >= 0, 'm must be a non-negative scalar.');

    limiter = @(a1, a2, a3, h) minmod_tilde(a1, a2, a3, m * h ^ 2);
end
