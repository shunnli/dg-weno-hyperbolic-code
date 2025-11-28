function lim = fv_limiter(raw_lim)
    lim = @(ul_plus, ur_minus, ul, u, ur, dx) fv_limiter_inner(ul_plus, ur_minus, ul, u, ur, dx, raw_lim);
end

function [ul_plus_mod, ur_minus_mod] = fv_limiter_inner(ul_plus, ur_minus, ul, u, ur, dx, raw_lim)
    ul_plus_mod = u - raw_lim(u - ul_plus, u - ul, ur - u, dx);
    ur_minus_mod = u + raw_lim(ur_minus - u, u - ul, ur - u, dx);
end
