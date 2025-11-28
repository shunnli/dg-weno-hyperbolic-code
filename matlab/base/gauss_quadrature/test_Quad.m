function tests = test_Quad()
    % TESTQUAD Unit tests for the Quad class
    tests = functiontests(localfunctions);
end

%% Normal usage tests

function testDefaultConstructor(testCase)
    % Test default constructor
    q = Quad();
    verifyClass(testCase, q, 'Quad');
    verifyEqual(testCase, numel(q.nodes), 5);
    verifyEqual(testCase, numel(q.weights), 5);
end

function testGaussLegendreConstructor(testCase)
    % Test rule-based constructor with parameters
    q = Quad('type', 'GaussLegendre', 'gk', 7);
    verifyClass(testCase, q, 'Quad');
    verifyEqual(testCase, numel(q.nodes), 7);
    verifyEqual(testCase, numel(q.weights), 7);
end

function testManualNodesWeightsConstructor(testCase)
    % Test manual specification of nodes and weights
    nodes = [-0.5; 0.5];
    weights = [1; 1];
    q = Quad('nodes', nodes, 'weights', weights);
    verifyClass(testCase, q, 'Quad');
    verifyEqual(testCase, q.nodes, nodes);
    verifyEqual(testCase, q.weights, weights);
end

function testIntegrationBasic(testCase)
    % Test basic integration functionality
    q = Quad('gk', 3);
    f = @(x) x.^2;
    xleft = 0;
    xright = 1;
    result = q.integrate(f, xleft, xright);
    verifyGreaterThan(testCase, result, 0.3);
    verifyLessThan(testCase, result, 0.4);
end

function testIntegrationVectorInput(testCase)
    % Test vector input for integration intervals
    q = Quad('gk', 4);
    f = @(x) cos(x);
    xleft = [0, pi];
    xright = [pi, 2*pi];

    result = q.integrate(f, xleft, xright);
    verifySize(testCase, result, size(xleft));

    result2 = q.integrate(f, xleft', xright');
    verifySize(testCase, result2, size(xleft'));

    verifyEqual(testCase, result, result2');
end

%% Error handling tests (assert failures)

function testInvalidNodesWeightsLength(testCase)
    % nodes and weights size mismatch
    fcn = @() Quad('nodes', [0, 1], 'weights', 1);
    verifyAssertionFailure(testCase, fcn, 'nodes and weights must have the same size.');
end

function testInvalidIntegrationInputs(testCase)
    % Invalid input arguments in integrate()
    q = Quad();

    fcn1 = @() q.integrate(5, 0, 1);
    verifyAssertionFailure(testCase, fcn1, 'f must be a function handle.');

    fcn2 = @() q.integrate(@(x) x, [0, 1], 1);
    verifyAssertionFailure(testCase, fcn2, 'xleft and xright must have the same size.');
end

%% Helper function for assertion failure check

function verifyAssertionFailure(testCase, fcn, expectedMessage)
    % Helper to verify assertion failures (no identifier)
    try
        fcn();
        verifyFail(testCase, 'Expected an assertion failure but none occurred.');
    catch ME
        verifyMatches(testCase, ME.message, expectedMessage);
    end
end
