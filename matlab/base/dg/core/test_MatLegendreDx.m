function tests = test_MatLegendreDx()
    % TEST_MATLEGENDREDX Test suite for MatLegendreDx class (Legendre polynomial derivatives)
    tests = functiontests(localfunctions);
end

%% Test MatLegendreDx Constructor
function testConstructor(testCase)
    % Test case 1: Construct MatLegendreDx with n = 5
    n = 5;
    basis = MatLegendreDx(n);

    verifyClass(testCase, basis, 'MatLegendreDx');
    verifyEqual(testCase, basis.funcs_num, n, ...
    'Incorrect number of basis functions stored.');
end

%% Test Evaluation of Legendre Polynomial Derivatives
function testEvalFunction(testCase)
    % Test case 2: Evaluate the first 5 derivatives of Legendre polynomials at test points
    n = 5;
    basis = MatLegendreDx(n);
    x = linspace(-0.9, 0.9, 5)';

    u = basis.eval(x, n);

    expected_u = [
                  0, 1, -2.7, 4.575, -6.0075,;
                  0, 1, -1.35, 0.01875, 1.7803125,;
                  0, 1, 0, -1.5, 0;
                  0, 1, 1.35, 0.01875 -1.7803125;
                  0, 1, 2.7, 4.575, 6.0075;
                  ];

    verifySize(testCase, u, size(expected_u));
    verifyEqual(testCase, u, expected_u, 'AbsTol', 1e-6);
end

%% Test Eval Function with Row Vector Input
function testEvalRowVector(testCase)
    % Test case 3: Verify eval handles row vectors correctly
    n = 5;
    basis = MatLegendreDx(n);
    x_row = linspace(-0.9, 0.9, 5);

    u_row = basis.eval(x_row, n);
    expected_u_row = [
                      0, 1, -2.7, 4.575, -6.0075,;
                      0, 1, -1.35, 0.01875, 1.7803125,;
                      0, 1, 0, -1.5, 0;
                      0, 1, 1.35, 0.01875 -1.7803125;
                      0, 1, 2.7, 4.575, 6.0075;
                      ]';

    verifyEqual(testCase, u_row, expected_u_row, 'AbsTol', 1e-6);
end

%% Test Invalid Constructor Input
function testConstructorInvalidInput(testCase)
    % Test case 4: Verify constructor rejects invalid inputs

    verifyError(testCase, @() MatLegendreDx(0), 'MatLegendreDx:InvalidInput');

    verifyError(testCase, @() MatLegendreDx(-3), 'MatLegendreDx:InvalidInput');

    verifyError(testCase, @() MatLegendreDx(3.5), 'MatLegendreDx:InvalidInput');
end
