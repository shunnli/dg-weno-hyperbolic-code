function tests = test_MatLagrange()
    % TEST_MATLAGRANGE Test suite for MatLagrange class using function-based tests.
    tests = functiontests(localfunctions);
end

%% Test MatLagrange Constructor
function testConstructor(testCase)
    % Test case 1: Construct MatLagrange with a set of points
    points = [-1, 0, 1];
    basis = MatLagrange(points);

    verifyClass(testCase, basis, 'MatLagrange');
    verifyEqual(testCase, numel(basis.funcs), numel(points), ...
    'Incorrect number of basis functions stored.');
end

%% Test Evaluation of Lagrange Polynomials
function testEvalFunction(testCase)
    % Test case 2: Evaluate the Lagrange polynomials at test points
    points = [-1, 0, 1];
    basis = MatLagrange(points);
    x = linspace(-1, 1, 11)';

    u = basis.eval(x, numel(points));

    expected_u = [
                  1., 0., 0.;
                  0.72, 0.36, -0.08;
                  0.48, 0.64, -0.12;
                  0.28, 0.84, -0.12;
                  0.12, 0.96, -0.08;
                  0., 1., 0.;
                  -0.08, 0.96, 0.12;
                  -0.12, 0.84, 0.28;
                  -0.12, 0.64, 0.48;
                  -0.08, 0.36, 0.72;
                  0., 0., 1.
                  ];

    verifySize(testCase, u, size(expected_u));
    verifyEqual(testCase, u, expected_u, 'AbsTol', 1e-6);
end

%% Test Eval Function with Row Vector Input
function testEvalRowVector(testCase)
    % Test case 3: Verify eval handles row vectors correctly
    points = [-1, 0, 1];
    basis = MatLagrange(points);
    x_row = linspace(-1, 1, 11);

    u_row = basis.eval(x_row, numel(points));
    expected_u_row = [
                      1., 0., 0.;
                      0.72, 0.36, -0.08;
                      0.48, 0.64, -0.12;
                      0.28, 0.84, -0.12;
                      0.12, 0.96, -0.08;
                      0., 1., 0.;
                      -0.08, 0.96, 0.12;
                      -0.12, 0.84, 0.28;
                      -0.12, 0.64, 0.48;
                      -0.08, 0.36, 0.72;
                      0., 0., 1.
                      ]';

    verifyEqual(testCase, u_row, expected_u_row, 'AbsTol', 1e-6);
end

%% Test Eval Function with Large Number of Points
function testEvalLargeInput(testCase)
    % Test case 4: Verify eval function with large number of points
    points_large = linspace(-1, 1, 20);
    basis_large = MatLagrange(points_large);
    x_large = linspace(-1, 1, 50)';

    u_large = basis_large.eval(x_large, numel(points_large));

    verifySize(testCase, u_large, [numel(x_large), numel(points_large)]);
end

%% Test Invalid Constructor Input
function testConstructorInvalidInput(testCase)
    % Test case 5: Verify constructor rejects invalid input
    verifyError(testCase, @() MatLagrange(1), 'MatLagrange:InvalidInput');
end
