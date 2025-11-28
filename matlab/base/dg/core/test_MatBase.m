function tests = test_MatBase()
    % TEST_MATBASE Test suite for MatBase class using function-based tests.
    tests = functiontests(localfunctions);
end

%% Test MatBase Constructor
function testConstructor(testCase)
    % Test case 1: Construct MatBase with simple basis functions
    funcs = {@(x) ones(size(x)), @(x) x, @(x) x .^ 2};
    basis = MatBase(funcs);

    verifyClass(testCase, basis, 'MatBase');
    verifyEqual(testCase, basis.funcs_num, numel(funcs), ...
    'Incorrect number of basis functions stored.');
end

%% Test Evaluation of Basis Functions
function testEvalFunction(testCase)
    % Test case 2: Evaluate the basis functions at test points
    funcs = {@(x) ones(size(x)), @(x) x, @(x) x .^ 2};
    basis = MatBase(funcs);

    x = linspace(-1, 1, 5)';
    n = 3;
    u = basis.eval(x, n);

    expected_u = [
                  1, -1, 1;
                  1, -0.5, 0.25;
                  1, 0, 0;
                  1, 0.5, 0.25;
                  1, 1, 1
                  ];

    verifySize(testCase, u, size(expected_u));
    verifyEqual(testCase, u, expected_u, 'AbsTol', 1e-6);
end

%% Test Eval Function with Row Vector Input
function testEvalRowVector(testCase)
    % Test case 3: Verify eval handles row vectors correctly
    funcs = {@(x) ones(size(x)), @(x) x, @(x) x .^ 2};
    basis = MatBase(funcs);

    x_row = linspace(-1, 1, 5);
    u_row = basis.eval(x_row, 3);

    expected_u = [
                  1, -1, 1;
                  1, -0.5, 0.25;
                  1, 0, 0;
                  1, 0.5, 0.25;
                  1, 1, 1
                  ]';

    verifyEqual(testCase, u_row, expected_u, 'AbsTol', 1e-6);
end

%% Test Invalid Input (n > funcs_num)
function testEvalInvalidN(testCase)
    % Test case 4: Invalid input (n > funcs_num)
    funcs = {@(x) ones(size(x)), @(x) x, @(x) x .^ 2};
    basis = MatBase(funcs);

    x = linspace(-1, 1, 5)';

    verifyError(testCase, @() basis.eval(x, 4), 'MatBase:InvalidInput');
end
