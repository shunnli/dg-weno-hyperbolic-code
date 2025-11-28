function [L, D, R] = eulereqs_DF(v1, v2, v3)

    [rho, u, p] = eulereqs_trans2raw(v1, v2, v3);

    gamma = 1.4;
    E = 1 / (gamma - 1) * p + 0.5 * rho .* u .^ 2;
    c = sqrt(gamma * abs(p ./ rho));
    H = (E + p) ./ rho;

    R = [1, 1, 1;
         u - c, u, u + c;
         H - u .* c, 0.5 * u .^ 2, H + u .* c;
         ];

    L = inv(R);

    D = diag([u - c, u, u + c]);
end
