function metrics = evaluateRecoveryOverTime(filterResult, topology, cfg)
%EVALUATERECOVERYOVERTIME Compute requested recovery metrics over time.

if isfield(filterResult, 'updatedState')
    estimatedState = filterResult.updatedState;
    estimatorName = "Dynamic";
elseif isfield(filterResult, 'estimatedState')
    estimatedState = filterResult.estimatedState;
    estimatorName = "Static";
else
    error('dynsc:MissingEstimatedState', ...
        ['The estimator result must contain updatedState or ' ...
         'estimatedState.']);
end
if isfield(filterResult, 'name') && strlength(string(filterResult.name)) > 0
    estimatorName = string(filterResult.name);
end

nEdges = size(topology.s1, 1);
nTriangles = size(topology.s2, 1);
nTime = size(estimatedState, 2);

assert(size(estimatedState, 1) == nEdges + nTriangles, ...
    'dynsc:MetricsStateSizeMismatch', ...
    'The updated state does not match the ground-truth topology.');

assert(size(topology.s1, 2) == nTime && ...
    size(topology.s2, 2) == nTime, ...
    'dynsc:MetricsTopologyTimeMismatch', ...
    ['The ground-truth topology and filter result must contain the same ' ...
     'number of time steps.']);

trueS1 = double(topology.s1);
trueS2 = double(topology.s2);

edgeF1 = zeros(1, nTime);
triangleF1 = zeros(1, nTime);
edgeNormalizedSHD = zeros(1, nTime);
triangleNormalizedSHD = zeros(1, nTime);
edgeRelativeError = zeros(1, nTime);
triangleRelativeError = zeros(1, nTime);

edgeRows = 1:nEdges;
triangleRows = nEdges + (1:nTriangles);

for iTime = 1:nTime
    estimatedS1 = estimatedState(edgeRows, iTime);
    estimatedS2 = estimatedState(triangleRows, iTime);
    trueS1AtTime = trueS1(:, iTime);
    trueS2AtTime = trueS2(:, iTime);

    estimatedEdgeSupport = estimatedS1 >= cfg.supportThreshold;
    estimatedTriangleSupport = estimatedS2 >= cfg.supportThreshold;

    edgeF1(iTime) = supportF1( ...
        trueS1AtTime > 0, estimatedEdgeSupport);
    triangleF1(iTime) = supportF1( ...
        trueS2AtTime > 0, estimatedTriangleSupport);
    edgeNormalizedSHD(iTime) = normalizedSHD( ...
        trueS1AtTime > 0, estimatedEdgeSupport);
    triangleNormalizedSHD(iTime) = normalizedSHD( ...
        trueS2AtTime > 0, estimatedTriangleSupport);
    edgeRelativeError(iTime) = ...
        relativeError(trueS1AtTime, estimatedS1);
    triangleRelativeError(iTime) = ...
        relativeError(trueS2AtTime, estimatedS2);
end

metrics = struct();
metrics.estimatorName = estimatorName;
metrics.supportThreshold = cfg.supportThreshold;
metrics.time = 1:nTime;
metrics.topologyChangeTimes = topology.changeTimes(:)';
metrics.topologyChanged = ismember(metrics.time, ...
    metrics.topologyChangeTimes);
metrics.topologyChangeLog = topology.changeLog;
if isfield(filterResult, 'priorResidual')
    assert(isequal(size(filterResult.priorResidual), ...
        size(estimatedState)), ...
        'dynsc:MetricsPriorResidualSizeMismatch', ...
        ['The prior residual must have the same size as the estimated ' ...
         'state.']);
    metrics.priorResidualNorm = ...
        vecnorm(filterResult.priorResidual, 2, 1);
else
    metrics.priorResidualNorm = zeros(1, 0);
end
metrics.edge.supportF1 = edgeF1;
metrics.edge.normalizedSHD = edgeNormalizedSHD;
metrics.edge.relativeError = edgeRelativeError;
metrics.triangle.supportF1 = triangleF1;
metrics.triangle.normalizedSHD = triangleNormalizedSHD;
metrics.triangle.relativeError = triangleRelativeError;
end

function distance = normalizedSHD(truth, estimate)
% Fraction of candidate simplices whose binary support is incorrect.
truth = logical(truth(:));
estimate = logical(estimate(:));
assert(numel(truth) == numel(estimate), ...
    'dynsc:NormalizedSHDSizeMismatch', ...
    'Truth and estimated supports must have the same size.');

if isempty(truth)
    distance = 0;
else
    distance = nnz(xor(truth, estimate)) / numel(truth);
end
end

function f1 = supportF1(truth, estimate)
truth = logical(truth(:));
estimate = logical(estimate(:));

truePositive = nnz(truth & estimate);
falsePositive = nnz(~truth & estimate);
falseNegative = nnz(truth & ~estimate);
denominator = 2 * truePositive + falsePositive + falseNegative;

if denominator == 0
    f1 = 1;
else
    f1 = 2 * truePositive / denominator;
end
end

function error = relativeError(truth, estimate)
denominator = norm(truth, 2);

if denominator > 0
    error = norm(estimate - truth, 2) / denominator;
elseif norm(estimate, 2) == 0
    error = 0;
else
    error = Inf;
end
end
