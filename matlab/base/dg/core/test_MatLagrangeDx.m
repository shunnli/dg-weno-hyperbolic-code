function tests = test_MatLagrangeDx()
    % TEST_MATLAGRANGEDX Test suite for MatLagrangeDx class (derivatives of Lagrange basis functions).
    tests = functiontests(localfunctions);
end

%% Test MatLagrangeDx Constructor
function testConstructor(testCase)
    % Test case 1: Construct MatLagrangeDx with a set of points
    points = [-1, 0, 1];
    basisDx = MatLagrangeDx(points);

    verifyClass(testCase, basisDx, 'MatLagrangeDx');
    verifyEqual(testCase, numel(basisDx.funcs), numel(points), ...
    'Incorrect number of derivative basis functions stored.');
end

%% Test Evaluation of Lagrange Derivative Polynomials
function testEvalFunction(testCase)
    % Test case 2: Evaluate the derivative of the Lagrange polynomials at test points
    points = [-1, 0, 1];
    basisDx = MatLagrangeDx(points);
    x = linspace(-1, 1, 11)';

    du = basisDx.eval(x, numel(points));

    expected_du = [
                   -1.5, 2., -0.5;
                   -1.3, 1.6, -0.3;
                   -1.1, 1.2, -0.1;
                   -0.9, 0.8, 0.1;
                   -0.7, 0.4, 0.3;
                   -0.5, 0., 0.5;
                   -0.3, -0.4, 0.7;
                   -0.1, -0.8, 0.9;
                   0.1, -1.2, 1.1;
                   0.3, -1.6, 1.3;
                   0.5, -2., 1.5
                   ];

    verifySize(testCase, du, size(expected_du));
    verifyEqual(testCase, du, expected_du, 'AbsTol', 1e-6);
end

%% Test Eval Function with Row Vector Input
function testEvalRowVector(testCase)
    % Test case 3: Verify eval handles row vectors correctly
    points = [-1, 0, 1];
    basisDx = MatLagrangeDx(points);
    x_row = linspace(-1, 1, 11);

    du_row = basisDx.eval(x_row, numel(points));
    expected_du_row = [
                       -1.5, 2., -0.5;
                       -1.3, 1.6, -0.3;
                       -1.1, 1.2, -0.1;
                       -0.9, 0.8, 0.1;
                       -0.7, 0.4, 0.3;
                       -0.5, 0., 0.5;
                       -0.3, -0.4, 0.7;
                       -0.1, -0.8, 0.9;
                       0.1, -1.2, 1.1;
                       0.3, -1.6, 1.3;
                       0.5, -2., 1.5
                       ]';

    verifyEqual(testCase, du_row, expected_du_row, 'AbsTol', 1e-6);
end

%% Test Eval Function with Large Number of Points
function testEvalLargeInput(testCase)
    % Test case 4: Verify eval function with large number of points
    points_large = linspace(-1, 1, 20);
    basisDx_large = MatLagrangeDx(points_large);
    x_large = linspace(-1, 1, 50)';

    du_large = basisDx_large.eval(x_large, numel(points_large));

    verifySize(testCase, du_large, [numel(x_large), numel(points_large)]);
end

%% Test Invalid Constructor Input
function testConstructorInvalidInput(testCase)
    % Test case 5: Verify constructor rejects invalid input
    verifyError(testCase, @() MatLagrangeDx(1), 'MatLagrangeDx:InvalidInput');
end
