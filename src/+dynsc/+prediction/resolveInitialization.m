function [cfg, details] = resolveInitialization(cfg, ~, topology)
%RESOLVEINITIALIZATION Build s0 and P0 from the configured initialization.

method = initializationMethod(cfg);
nEdges = size(topology.B1_full, 2);
nTriangles = size(topology.B2_full, 2);
nState = nEdges + nTriangles;

details = struct();
details.method = method;

switch method
    case "manual"
        cfg = dynsc.prediction.validateParameters(cfg, nState);

    case {"bernoulli", "ground_truth"}
        assert(isfield(topology, 'p'), ...
            'dynsc:MissingInitializationEdgeProbability', ...
            ['Bernoulli and ground-truth initialization require the ' ...
             'Erdos-Renyi edge probability in topology.p.']);
        assert(isfield(topology, 'triangleFraction'), ...
            'dynsc:MissingInitializationTriangleFraction', ...
            ['Bernoulli and ground-truth initialization require the ' ...
             'conditional triangle-selection probability in ' ...
             'topology.triangleFraction.']);

        varianceFloor = cfg.initialization.varianceFloor;
        validateattributes(varianceFloor, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'});

        edgeProbability = double(topology.p);
        triangleSelectionProbability = ...
            double(topology.triangleFraction);
        validateattributes(edgeProbability, {'numeric'}, ...
            {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
        validateattributes(triangleSelectionProbability, {'numeric'}, ...
            {'scalar', 'real', 'finite', '>=', 0, '<=', 1});

        % Under the independent-edge Erdos-Renyi approximation, a candidate
        % triangle is feasible only when all three of its boundary edges are
        % present. Conditional on feasibility, it is filled with probability
        % given by triangleFraction. Hence its marginal probability is
        % P(Y_tau = 1) = P(selected | feasible) P(feasible).
        triangleFeasibilityProbability = edgeProbability ^ 3;
        triangleProbability = triangleSelectionProbability * ...
            triangleFeasibilityProbability;

        edgeVariance = max( ...
            edgeProbability * (1 - edgeProbability), varianceFloor);
        triangleVariance = max( ...
            triangleProbability * (1 - triangleProbability), ...
            varianceFloor);

        if method == "ground_truth"
            assert(isfield(topology, 's1') && isfield(topology, 's2'), ...
                'dynsc:MissingGroundTruthInitialization', ...
                ['Ground-truth initialization requires topology.s1 ' ...
                 'and topology.s2.']);
            assert(size(topology.s1, 1) == nEdges && ...
                size(topology.s2, 1) == nTriangles && ...
                size(topology.s1, 2) >= 1 && size(topology.s2, 2) >= 1, ...
                'dynsc:InvalidGroundTruthInitialization', ...
                ['The first columns of topology.s1 and topology.s2 must ' ...
                 'match the candidate edge and triangle dimensions.']);

            initialState = double([topology.s1(:, 1); topology.s2(:, 1)]);
            assert(all(isfinite(initialState)) && ...
                all(initialState == 0 | initialState == 1), ...
                'dynsc:InvalidGroundTruthInitialization', ...
                'The ground-truth initial topology must be binary.');
            details.stateSource = "ground_truth";
        else
            initialState = [
                edgeProbability * ones(nEdges, 1)
                triangleProbability * ones(nTriangles, 1)
            ];
            details.stateSource = "bernoulli_mean";
        end

        cfg.initialState = initialState;
        cfg.initialCovariance = diag([
            edgeVariance * ones(nEdges, 1)
            triangleVariance * ones(nTriangles, 1)
        ]);

        details.numCandidateEdges = nEdges;
        details.numCandidateTriangles = nTriangles;
        details.edgeProbability = edgeProbability;
        details.triangleFeasibilityProbability = ...
            triangleFeasibilityProbability;
        details.triangleSelectionProbability = ...
            triangleSelectionProbability;
        details.triangleProbability = triangleProbability;
        details.expectedC1 = nEdges * edgeProbability;
        details.expectedC2 = nTriangles * triangleProbability;
        details.edgeVariance = edgeVariance;
        details.triangleVariance = triangleVariance;
        details.varianceFloor = varianceFloor;

        cfg = dynsc.prediction.validateParameters(cfg, nState);

    otherwise
        error('dynsc:InvalidInitializationMethod', ...
            ['cfg.prediction.initialization.method must be "bernoulli" ' ...
             '"ground_truth", or "manual".']);
end
end

function method = initializationMethod(cfg)
if isfield(cfg, 'initialization') && ...
        isfield(cfg.initialization, 'method')
    method = string(cfg.initialization.method);
else
    method = "manual";
end

assert(isscalar(method), 'dynsc:InvalidInitializationMethod', ...
    'The prediction initialization method must be scalar text.');
end
