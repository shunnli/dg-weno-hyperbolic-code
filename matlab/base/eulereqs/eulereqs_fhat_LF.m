function f = eulereqs_fhat_LF(u_minus, u_plus)
    alpha1 = eulereqs_get_alpha(u_minus);
    alpha2 = eulereqs_get_alpha(u_plus);
    alpha = max(alpha1, alpha2);

    f = 0.5 * (eulereqs_f(u_minus) + eulereqs_f(u_plus)) - 0.5 * alpha .* (u_plus - u_minus);
end
