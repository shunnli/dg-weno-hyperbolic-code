function [v1, v2, v3] = eulereqs_trans2cv(rho, u, p)
    % [rho, u, p] -> [v1, v2, v3]

    gamma = 1.4;
    v1 = rho;
    v2 = rho .* u;
    v3 = 1 / (gamma - 1) * p + 0.5 * rho .* u .^ 2;
end
