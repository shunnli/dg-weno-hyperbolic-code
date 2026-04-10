function [x, dx] = mesh_init_1d(xleft, xright, n)
    % MESH_INIT_1D Generates a uniform 1D mesh with n cells.
    %
    % INPUT:
    %   xleft   - Left boundary coordinate
    %   xright  - Right boundary coordinate (xleft < xright)
    %   n       - Number of cells (integer, must be positive)
    %
    % OUTPUT:
    %   x       - 1D array of cell center coordinates, size [1, n]
    %   dx      - Cell width
    %
    % DESCRIPTION:
    %   x(j) is the center of cell I_j = [x(j)-dx/2, x(j)+dx/2]
    %
    % EXAMPLE:
    %   [x, dx] = mesh_init_1d(0, 1, 4);
    %   % x = [0.125, 0.375, 0.625, 0.875]
    %   % dx = 0.25

    arguments
        xleft (1, 1) double
        xright (1, 1) double
        n (1, 1) double {mustBeInteger, mustBePositive}
    end

    assert(xleft < xright, 'xleft must be less than xright.');

    x = linspace(xleft, xright, n + 1);
    dx = x(2) - x(1);
    x = x(1:end - 1) + dx / 2;
end
