function [rho_out, u_out, p_out, more_info] = euler_riemann_exact( ...
        rho_l, u_l, p_l, rho_r, u_r, p_r, gamma, xlist, x_c, t)
    % EULER_RIEMANN_EXACT Computes the exact Riemann solution for the 1D Euler equations.
    %
    % This function calculates the exact solution of the Riemann problem for
    % the one-dimensional Euler equations in gas dynamics. The solution consists
    % of three waves: a left-moving wave, a contact discontinuity, and a right-moving wave.
    %
    % INPUT:
    %   rho_l  - Left state density (must be positive)
    %   u_l    - Left state velocity
    %   p_l    - Left state pressure (must be positive)
    %   rho_r  - Right state density (must be positive)
    %   u_r    - Right state velocity
    %   p_r    - Right state pressure (must be positive)
    %   gamma  - Specific heat ratio (= 1.4)
    %   xlist  - Array of spatial positions
    %   x_c    - Contact discontinuity location
    %   t      - Time at which to compute the solution
    %
    % OUTPUT:
    %   rho_out - Density distribution at given positions, same size as xlist
    %   u_out   - Velocity distribution at given positions, same size as xlist
    %   p_out   - Pressure distribution at given positions, same size as xlist
    %   more_info - Structure containing additional solution details:
    %               - p_s, u_s, rho_s_l, rho_s_r (intermediate states)
    %               - Wave speeds (w_1_l, w_1_r, w_2, w_3_l, w_3_r)
    %               - Wave types (shock or rarefaction)

    assert(isnumeric(rho_l) && isscalar(rho_l) && rho_l > 0, 'rho_l must be positive.');
    assert(isnumeric(u_l) && isscalar(u_l), 'u_l must be a scalar.');
    assert(isnumeric(p_l) && isscalar(p_l) && p_l > 0, 'p_l must be positive.');
    assert(isnumeric(rho_r) && isscalar(rho_r) && rho_r > 0, 'rho_r must be positive.');
    assert(isnumeric(u_r) && isscalar(u_r), 'u_r must be a scalar.');
    assert(isnumeric(p_r) && isscalar(p_r) && p_r > 0, 'p_r must be positive.');
    assert(isnumeric(gamma) && isscalar(gamma) && gamma > 1, 'gamma must be greater than 1.');
    assert(isnumeric(xlist) && isvector(xlist), 'xlist must be a vector.');
    assert(isnumeric(x_c) && isscalar(x_c), 'x_c must be a scalar.');
    assert(isnumeric(t) && isscalar(t) && t >= 0, 't must be non-negative.');

    % Compute the sound speeds in the left and right states
    c_l = sqrt(gamma * p_l / rho_l);
    c_r = sqrt(gamma * p_r / rho_r);

    % Compute constants alpha and beta
    alpha = (gamma + 1.0) / (gamma - 1.0);
    beta = (gamma - 1.0) / (2.0 * gamma);

    % Check for cavitation (not supported)
    if u_l - u_r + 2 * (c_l + c_r) / (gamma - 1.0) < 0
        disp('Cavitation detected! Exiting.');
        rho_out = []; u_out = []; p_out = []; more_info = [];
        return;
    end

    % Define integral curves and Hugoniot locus for 1-wave and 3-wave

    integral_curve_1 = @(p) u_l + 2 * c_l / (gamma - 1.0) * (1.0 - (p / p_l) ^ beta);
    integral_curve_3 = @(p) u_r - 2 * c_r / (gamma - 1.0) * (1.0 - (p / p_r) ^ beta);

    hugoniot_locus_1 = @(p) u_l + 2 * c_l / sqrt(2 * gamma * (gamma - 1.0)) * ...
        ((1 - p / p_l) / sqrt(1 + alpha * p / p_l));
    hugoniot_locus_3 = @(p) u_r - 2 * c_r / sqrt(2 * gamma * (gamma - 1.0)) * ...
        ((1 - p / p_r) / sqrt(1 + alpha * p / p_r));

    phi_l = @(p) (p >= p_l) * hugoniot_locus_1(p) + (p < p_l) * integral_curve_1(p);
    phi_r = @(p) (p >= p_r) * hugoniot_locus_3(p) + (p < p_r) * integral_curve_3(p);

    % Define the function for finding the intersection point in the (p-v) plane
    func = @(p) phi_l(p) - phi_r(p);

    % Initial guess for pressure in the middle state
    p0_PV = (p_l + p_r) / 2.0 -1/8 * (u_r - u_l) * (rho_l + rho_r) * (c_l + c_r);
    p0 = max(p0_PV, 1e-8);

    % Solve for the intersection point to get the intermediate state p and u
    options = optimset('TolFun', 1.0e-12, 'Display', 'off');
    [p_s_tmp, ~, exitflag] = fsolve(func, p0, options);
    p_s = p_s_tmp;
    u_s = 0.5 * (phi_l(p_s) + phi_r(p_s));

    % Warning if the solution fails to converge
    if exitflag <= 0
        disp('Warning: fsolve did not converge.');
    end

    % Compute the density on the left and right side of the contact discontinuity

    if p_s <= p_l
        rho_s_l = (p_s / p_l) ^ (1.0 / gamma) * rho_l; % Rarefaction wave
    else
        rho_s_l = ((1.0 + alpha * p_s / p_l) / ((p_s / p_l) + alpha)) * rho_l; % Shock wave
    end

    if p_s <= p_r
        rho_s_r = (p_s / p_r) ^ (1.0 / gamma) * rho_r; % Rarefaction wave
    else
        rho_s_r = ((1.0 + alpha * p_s / p_r) / ((p_s / p_r) + alpha)) * rho_r; % Shock wave
    end

    % Compute sound speed in the intermediate state
    c_s_l = sqrt(gamma * p_s / rho_s_l);
    c_s_r = sqrt(gamma * p_s / rho_s_r);

    % Compute wave speeds

    % 1-wave
    if p_s > p_l % Shock wave
        w_1_l = (rho_l * u_l - rho_s_l * u_s) / (rho_l - rho_s_l);
        w_1_r = w_1_l;
    else % Rarefaction wave
        w_1_l = u_l - c_l;
        w_1_r = u_s - c_s_l;
    end

    % 2-wave
    w_2 = u_s;

    % 3-wave
    if p_s > p_r % Shock wave
        w_3_l = (rho_r * u_r - rho_s_r * u_s) / (rho_r - rho_s_r);
        w_3_r = w_3_l;
    else % Rarefaction wave
        w_3_l = u_s + c_s_r;
        w_3_r = u_r + c_r;
    end

    w_max = max(abs([w_1_l, w_1_r, w_2, w_3_l, w_3_r]));

    % Warning if the wave is out of range
    if t * w_max > max(abs(xlist - x_c))
        disp('Warning: wave is out of range.');
    end

    % Solve for the state inside the rarefaction wave
    xi = (xlist - x_c) / t;
    u_1_fan = ((gamma - 1.0) * u_l + 2 * (c_l + xi)) / (gamma + 1.0);
    u_3_fan = ((gamma - 1.0) * u_r - 2 * (c_r - xi)) / (gamma + 1.0);
    rho_1_fan = (rho_l ^ gamma * (u_1_fan - xi) .^ 2 / (gamma * p_l)) .^ (1.0 / (gamma - 1.0));
    rho_3_fan = (rho_r ^ gamma * (xi - u_3_fan) .^ 2 / (gamma * p_r)) .^ (1.0 / (gamma - 1.0));
    p_1_fan = p_l * (rho_1_fan / rho_l) .^ gamma;
    p_3_fan = p_r * (rho_3_fan / rho_r) .^ gamma;

    % Default: Right state
    rho_out = rho_r * ones(size(xlist));
    u_out = u_r * ones(size(xlist));
    p_out = p_r * ones(size(xlist));

    % Left of the 1-wave: Left state
    idx_1_l = (xi <= w_1_l);
    rho_out(idx_1_l) = rho_l;
    u_out(idx_1_l) = u_l;
    p_out(idx_1_l) = p_l;

    % Inside the 1-wave (if it's a rarefaction wave)
    idx_1_r = (xi > w_1_l) & (xi <= w_1_r);
    rho_out(idx_1_r) = rho_1_fan(idx_1_r);
    u_out(idx_1_r) = u_1_fan(idx_1_r);
    p_out(idx_1_r) = p_1_fan(idx_1_r);

    % Between the 1-wave and the 2-wave
    idx_2 = (xi > w_1_r) & (xi <= w_2);
    rho_out(idx_2) = rho_s_l;
    u_out(idx_2) = u_s;
    p_out(idx_2) = p_s;

    % Between the 2-wave and the 3-wave
    idx_3_l = (xi > w_2) & (xi <= w_3_l);
    rho_out(idx_3_l) = rho_s_r;
    u_out(idx_3_l) = u_s;
    p_out(idx_3_l) = p_s;

    % Inside the 3-wave (if it's a rarefaction wave)
    idx_3_r = (xi > w_3_l) & (xi <= w_3_r);
    rho_out(idx_3_r) = rho_3_fan(idx_3_r);
    u_out(idx_3_r) = u_3_fan(idx_3_r);
    p_out(idx_3_r) = p_3_fan(idx_3_r);

    % Additional info
    more_info = struct();
    more_info.p_s = p_s;
    more_info.u_s = u_s;
    more_info.rho_s_l = rho_s_l;
    more_info.rho_s_r = rho_s_r;
    more_info.w_1_l = w_1_l;
    more_info.w_1_r = w_1_r;
    more_info.w_2 = w_2;
    more_info.w_3_l = w_3_l;
    more_info.w_3_r = w_3_r;

    if p_s > p_l
        left_type = 'Shock';
    else
        left_type = 'Rarefaction';
    end

    center_type = 'Contract discontinuity';

    if p_s > p_r
        right_type = 'Shock';
    else
        right_type = 'Rarefaction';
    end

    more_info.left_type = left_type;
    more_info.center_type = center_type;
    more_info.right_type = right_type;

    more_info.type_msg = sprintf('[L] %s\n[C] %s\n[R] %s', left_type, center_type, right_type);
    more_info.key_msg = sprintf('%s=%f\t%s=%f\t%s=%f\t%s=%f', 'p_s', p_s, 'u_s', u_s, 'rho_s_l', rho_s_l, 'rho_s_r', rho_s_r);
end
