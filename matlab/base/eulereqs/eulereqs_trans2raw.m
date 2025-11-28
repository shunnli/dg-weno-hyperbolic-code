function [rho, u, p] = eulereqs_trans2raw(v1, v2, v3)
    % [v1, v2, v3] -> [rho, u, p]

    gamma = 1.4;
    rho = v1;
    u = v2 ./ v1;
    p = (gamma - 1) * (v3 - 0.5 * (v2 .^ 2) ./ v1);
end
