function signals = generateFromTopology(topology, cfg)
%GENERATEFROMTOPOLOGY Sample fresh signals from each topology in the sequence.

signalParameters = struct();
signalParameters.M0 = cfg.numNodeSignals;
signalParameters.M1 = cfg.numEdgeSignals;
signalParameters.sigma = cfg.sigma;
signalParameters.norm_noise = cfg.normalizeNoise;
signalParameters.verbose = cfg.verbose;
signalParameters.sig_type = "smooth";

nTime = cfg.numTimeSteps;
nNodes = topology.N;
nFullEdges = size(topology.B2_full, 1);

assert(topology.numTimeSteps == nTime && ...
    size(topology.s1, 2) == nTime && ...
    size(topology.s2, 2) == nTime, ...
    'dynsc:SignalTopologyTimeMismatch', ...
    ['The topology support and signal configuration must contain the same ' ...
     'number of time steps.']);

signals = struct();
signals.time = 1:nTime;
signals.smoothnessType = string(cfg.smoothnessType);
signals.X_node = zeros(nNodes, cfg.numNodeSignals, nTime);
signals.Y_node = zeros(nNodes, cfg.numNodeSignals, nTime);
signals.X_edge = zeros(nFullEdges, cfg.numEdgeSignals, nTime);
signals.Y_edge = zeros(nFullEdges, cfg.numEdgeSignals, nTime);
% Retain the initial/static field and add the time-indexed sequence needed
% when the edge support evolves.
signals.nodeSmoothnessLaplacian = topology.L0;
signals.nodeSmoothnessLaplacianSequence = cell(1, nTime);
signals.edgeSmoothnessLaplacian = cell(1, nTime);
signals.topologyChangeTimes = topology.changeTimes;
boundaryEdges = dynsc.topology.resolveBoundaryEdges(topology);

for iTime = 1:nTime
    selectedEdges = logical(topology.s1(:, iTime));
    selectedTriangles = logical(topology.s2(:, iTime));
    B1AtTime = topology.B1_full(:, selectedEdges);
    nodeSmoothnessLaplacian = B1AtTime * B1AtTime';

    switch string(cfg.smoothnessType)
        case "low_curl"
            B2 = topology.B2_full(:, selectedTriangles);
            edgeSmoothnessLaplacian = B2 * B2';
        case "new_smoothness"
            edgeSmoothnessLaplacian = ...
                buildTriangleSmoothnessLaplacian( ...
                    selectedTriangles, boundaryEdges, nFullEdges);
        otherwise
            error('dynsc:UnknownSmoothnessType', ...
                'Unknown smoothness type: %s', cfg.smoothnessType);
    end

    signals.edgeSmoothnessLaplacian{iTime} = ...
        edgeSmoothnessLaplacian;
    signals.nodeSmoothnessLaplacianSequence{iTime} = ...
        nodeSmoothnessLaplacian;

    nodeOutput = generate_graph_signals0( ...
        nodeSmoothnessLaplacian, signalParameters);

    if string(cfg.smoothnessType) == "low_curl"
        edgeOutput = generate_graph_signals1( ...
            edgeSmoothnessLaplacian, signalParameters);
    else
        edgeOutput = sampleNewSmoothnessSignals( ...
            edgeSmoothnessLaplacian, cfg.numEdgeSignals);
    end

    signals.X_node(:, :, iTime) = nodeOutput.X_0;
    signals.Y_node(:, :, iTime) = nodeOutput.X;
    signals.X_edge(:, :, iTime) = edgeOutput.X_0;
    signals.Y_edge(:, :, iTime) = edgeOutput.X;
end
end

function out = sampleNewSmoothnessSignals(L, numSignals)
% Match the new_smoothness sampling block in generate_simplicial_complex.m.

alpha = 1;
nEdges = size(L, 1);
[eigenvectors, eigenvalues] = eig(full(L));
spectralCovariance = inv(alpha * eigenvalues + eye(nEdges));
spectralSignals = mvnrnd( ...
    zeros(nEdges, 1), spectralCovariance, numSignals)';
cleanSignals = eigenvectors * spectralSignals;

out = struct();
out.X_0 = cleanSignals;
out.X = cleanSignals;
end

function L = buildTriangleSmoothnessLaplacian( ...
    s2, boundaryEdges, nFullEdges)
% This is the local three-edge construction used by the original generator.

selectedTriangles = find(s2(:));
localLaplacian = [2, -1, -1; -1, 2, -1; -1, -1, 2];
L = sparse(nFullEdges, nFullEdges);

for iTriangle = 1:numel(selectedTriangles)
    triangleIndex = selectedTriangles(iTriangle);
    edgeIndices = boundaryEdges(triangleIndex, :);

    selector = sparse(1:3, edgeIndices, 1, 3, nFullEdges);
    L = L + selector' * localLaplacian * selector;
end
end
