function limiter = tvb(m)

    arguments
        m (1, 1) double {mustBeNonnegative}
    end

    limiter = @(a1, a2, a3, h) minmod_tilde(a1, a2, a3, m * h ^ 2);
end
