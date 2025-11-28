function euler_riemann_exact_plot( ...
        rho_l, u_l, p_l, rho_r, u_r, p_r, x_l, x_r, x_c, t)
    % EULER_RIEMANN_EXACT_PLOT Plots the exact solution of the Euler Riemann problem.
    %
    % This function computes and visualizes the exact solution of the Riemann problem
    % for the 1D Euler equations. It generates two plots:
    % 1. Primitive variable profiles (density, velocity, pressure, internal energy).
    % 2. Space-time diagram showing wave propagation.
    %
    % INPUT:
    %   rho_l  - Left state density (must be positive)
    %   u_l    - Left state velocity
    %   p_l    - Left state pressure (must be positive)
    %   rho_r  - Right state density (must be positive)
    %   u_r    - Right state velocity
    %   p_r    - Right state pressure (must be positive)
    %   x_l    - Left boundary of spatial domain
    %   x_r    - Right boundary of spatial domain (must be greater than x_l)
    %   x_c    - Initial location of the contact discontinuity (optional, default = 0.0)
    %   t      - Time at which to compute the solution (must be non-negative)
    %
    % EXAMPLE:
    %   euler_riemann_exact_plot(...
    %        1.0, 0.0, 1.0, ...
    %        0.125, 0.0, 0.1, ...
    %        0.0, 1.0, 0.5, 0.25);

    assert(isnumeric(rho_l) && isscalar(rho_l) && rho_l > 0, 'rho_l must be positive.');
    assert(isnumeric(u_l) && isscalar(u_l), 'u_l must be a scalar.');
    assert(isnumeric(p_l) && isscalar(p_l) && p_l > 0, 'p_l must be positive.');
    assert(isnumeric(rho_r) && isscalar(rho_r) && rho_r > 0, 'rho_r must be positive.');
    assert(isnumeric(u_r) && isscalar(u_r), 'u_r must be a scalar.');
    assert(isnumeric(p_r) && isscalar(p_r) && p_r > 0, 'p_r must be positive.');
    assert(isnumeric(x_l) && isscalar(x_l), 'x_l must be a scalar.');
    assert(isnumeric(x_r) && isscalar(x_r), 'x_r must be a scalar.');
    assert(isnumeric(x_c) && isscalar(x_c), 'x_c must be a scalar.');
    assert(x_l < x_c && x_c < x_r, 'x_l < x_c < x_r.');
    assert(isnumeric(t) && isscalar(t) && t >= 0, 't must be non-negative.');

    xlist = linspace(x_l, x_r, 1000);

    gamma = 1.4;
    [rho, u, p, more_info] = euler_riemann_exact( ...
        rho_l, u_l, p_l, rho_r, u_r, p_r, gamma, xlist, x_c, t);

    % Compute internal energy
    e = 1 / (gamma - 1) * p ./ rho;

    primitive = {rho, u, p, e};
    names = {'Density', 'Velocity', 'Pressure', 'Internal Energy'};

    figure;

    for w = 1:4
        subplot(2, 2, w);
        q = primitive{w};
        plot(xlist, q, 'LineWidth', 2);
        title(names{w});

        qmax = max(q);
        qmin = min(q);
        qdiff = qmax - qmin;
        ylim([qmin - 0.1 * qdiff, qmax + 0.1 * qdiff]);
    end

    % Plot the wave lines
    figure;
    ax = axes;
    hold on;
    tlist = linspace(0, t, 500);

    % Plot left shock or rarefaction
    if strcmp(more_info.left_type, 'Shock')
        w_1_line = x_c + more_info.w_1_l * tlist;
        plot(w_1_line, tlist, 'r', 'DisplayName', 'Left shock');
    else
        w_1_l_line = x_c + more_info.w_1_l * tlist;
        w_1_r_line = x_c + more_info.w_1_r * tlist;
        fill_betweenx(tlist, w_1_l_line, w_1_r_line, 'r', 0.2, 'Left Rarefaction');
    end

    % Contact discontinuity in the middle
    w_2_line = x_c + more_info.w_2 * tlist;
    plot(w_2_line, tlist, 'g', 'DisplayName', 'Contract discontinuity');

    % Plot right shock or rarefaction
    if strcmp(more_info.right_type, 'Shock')
        w_3_line = x_c + more_info.w_3_l * tlist;
        plot(w_3_line, tlist, 'b', 'DisplayName', 'Right shock');
    else
        w_3_l_line = x_c + more_info.w_3_l * tlist;
        w_3_r_line = x_c + more_info.w_3_r * tlist;
        fill_betweenx(tlist, w_3_l_line, w_3_r_line, 'b', 0.2, 'Right Rarefaction');
    end

    % Adjust axes and labels
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.XTick = [x_l, x_c, x_r];
    ax.YTick = t;
    ax.YTickLabel = {'t'};
    xlim([x_l, x_r]);
    ylim([0, t]);
    legend;
    hold off;

    subtitle(sprintf('Riemann Problem (t = %.4f)', t));

end

function fill_betweenx(tlist, y1, y2, color, alpha, label)
    % Fill between two curves, adjusted for x and y
    patch([y1, fliplr(y2)], [tlist, fliplr(tlist)], color, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none', 'DisplayName', label);
end
