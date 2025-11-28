clc;
close all;

fprintf('Example 1: Using manual nodes & weights\n');
nodes = [-sqrt(1/3); sqrt(1/3)];
weights = [1; 1];

quad1 = Quad('nodes', nodes, 'weights', weights);

f = @(x) x .^ 2;
xleft = [0, 1, 2, 3];
xright = [1, 2, 3, 4];

result = quad1.integrate(f, xleft, xright);
exact = (xright .^ 3 - xleft .^ 3) ./ 3;
fprintf('Approximate: %s\n', mat2str(result));
fprintf('Exact: %s\n', mat2str(exact));
fprintf('Error: %s\n\n', mat2str(abs(result - exact)));

fprintf('Example 2: Using Gauss-Legendre rule (gk = 3)\n');
quad2 = Quad('type', 'GaussLegendre', 'gk', 3);

f = @(x) sin(x);
xleft = 0;
xright = pi;
result = quad2.integrate(f, xleft, xright);
exact = 2;
fprintf('Approximate: %.12f\n', result);
fprintf('Exact: %.12f\n', exact);
fprintf('Error: %.2e\n\n', abs(result - exact));

fprintf('Example 3: Using Gauss-Lobatto rule (gk = 5)\n');
quad3 = Quad('type', 'GaussLobatto', 'gk', 5);

f = @(x) exp(x);
xleft = [-1; 0; 3];
xright = [1; 2; 4];
result = quad3.integrate(f, xleft, xright);
exact = exp(xright) - exp(xleft);
fprintf('Approximate:\n');
disp(result);
fprintf('Exact:\n');
disp(exact);
fprintf('Error:\n');
disp(abs(result - exact));
