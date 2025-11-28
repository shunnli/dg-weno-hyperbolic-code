function tests = test_MatLegendre()
    % TEST_MATLEGENDRE Test suite for MatLegendre class (Legendre polynomials)
    tests = functiontests(localfunctions);
end

%% Test MatLegendre Constructor
function testConstructor(testCase)
    % Test case 1: Construct MatLegendre with n = 5
    n = 5;
    basis = MatLegendre(n);

    verifyClass(testCase, basis, 'MatLegendre');
    verifyEqual(testCase, basis.funcs_num, n, ...
    'Incorrect number of basis functions stored.');
end

%% Test Evaluation of Legendre Polynomials
function testEvalFunction(testCase)
    % Test case 2: Evaluate the first 5 Legendre polynomials at test points
    n = 5;
    basis = MatLegendre(n);
    x = linspace(-1, 1, 5)';

    u = basis.eval(x, n);

    expected_u = [
                  1, -1, 1, -1, 1;
                  1., -0.5, -0.125, 0.4375, -0.289063;
                  1, 0, - (1/2), 0, 3/8;
                  1., 0.5, -0.125, -0.4375, -0.289063;
                  1, 1, 1, 1, 1
                  ];

    verifySize(testCase, u, size(expected_u));
    verifyEqual(testCase, u, expected_u, 'AbsTol', 1e-6);
end

%% Test Eval Function with Row Vector Input
function testEvalRowVector(testCase)
    % Test case 3: Verify eval handles row vectors correctly
    n = 5;
    basis = MatLegendre(n);
    x_row = linspace(-1, 1, 5);

    u_row = basis.eval(x_row, n);
    expected_u_row = [
                      1, -1, 1, -1, 1;
                      1., -0.5, -0.125, 0.4375, -0.289063;
                      1, 0, - (1/2), 0, 3/8;
                      1., 0.5, -0.125, -0.4375, -0.289063;
                      1, 1, 1, 1, 1
                      ]';

    verifyEqual(testCase, u_row, expected_u_row, 'AbsTol', 1e-6);
end

%% Test Invalid Constructor Input
function testConstructorInvalidInput(testCase)
    % Test case 4: Verify constructor rejects invalid inputs

    verifyError(testCase, @() MatLegendre(0), 'MatLegendre:InvalidInput');

    verifyError(testCase, @() MatLegendre(-3), 'MatLegendre:InvalidInput');

    verifyError(testCase, @() MatLegendre(3.5), 'MatLegendre:InvalidInput');
end
