function diagnostics = buildPrimarySchoolStateDiagnostics( ...
    filterResult, topology, metrics)
%BUILDPRIMARYSCHOOLSTATEDIAGNOSTICS Summarize each empirical-SC update.

state = filterResult.updatedState;
nEdges = size(topology.s1, 1);
nTriangles = size(topology.s2, 1);
nTime = size(state, 2);
assert(isequal(size(state), [nEdges + nTriangles, nTime]) && ...
    size(topology.s1, 2) == nTime && ...
    size(topology.s2, 2) == nTime, ...
    'dynsc:PrimarySchoolDiagnosticSizeMismatch', ...
    'The estimate and empirical topology must have matching sizes.');
assert(numel(metrics.time) == nTime, ...
    'dynsc:PrimarySchoolDiagnosticMetricSizeMismatch', ...
    'The recovery metrics must match the empirical topology sequence.');

edgeState = state(1:nEdges, :);
triangleState = state(nEdges + 1:end, :);
threshold = metrics.supportThreshold;

k = (1:nTime).';
if isfield(topology, 'windowStart') && ...
        numel(topology.windowStart) == nTime
    windowStart = topology.windowStart(:);
else
    windowStart = NaT(nTime, 1);
end
trueC1 = sum(topology.s1, 1).';
trueC2 = sum(topology.s2, 1).';
trueEdgeSupportEmpty = trueC1 == 0;
trueTriangleSupportEmpty = trueC2 == 0;

estimatedEdgeMass = sum(edgeState, 1).';
estimatedTriangleMass = sum(triangleState, 1).';
estimatedEdgeSupport = sum(edgeState >= threshold, 1).';
estimatedTriangleSupport = sum(triangleState >= threshold, 1).';
maximumEdgeWeight = max(edgeState, [], 1).';
maximumTriangleWeight = max(triangleState, [], 1).';

edgeSupportF1 = metrics.edge.supportF1(:);
triangleSupportF1 = metrics.triangle.supportF1(:);
edgeNormalizedSHD = metrics.edge.normalizedSHD(:);
triangleNormalizedSHD = metrics.triangle.normalizedSHD(:);
edgeRelativeError = metrics.edge.relativeError(:);
triangleRelativeError = metrics.triangle.relativeError(:);
edgeRelativeErrorFinite = isfinite(edgeRelativeError);
triangleRelativeErrorFinite = isfinite(triangleRelativeError);

priorResidualNorm = nan(nTime, 1);
if isfield(filterResult, 'priorResidual') && ...
        isequal(size(filterResult.priorResidual), size(state))
    priorResidualNorm = vecnorm(filterResult.priorResidual, 2, 1).';
end

smoothnessResidual = nan(nTime, 1);
edgeCountResidual = nan(nTime, 1);
triangleCountResidual = nan(nTime, 1);
if isfield(filterResult, 'measurementResidual') && ...
        size(filterResult.measurementResidual, 2) == nTime && ...
        size(filterResult.measurementResidual, 1) >= 3
    smoothnessResidual = filterResult.measurementResidual(1, :).';
    edgeCountResidual = filterResult.measurementResidual(2, :).';
    triangleCountResidual = filterResult.measurementResidual(3, :).';
end

diagnostics = table(k, windowStart, trueC1, trueC2, ...
    trueEdgeSupportEmpty, trueTriangleSupportEmpty, ...
    estimatedEdgeMass, estimatedTriangleMass, estimatedEdgeSupport, ...
    estimatedTriangleSupport, maximumEdgeWeight, maximumTriangleWeight, ...
    edgeSupportF1, triangleSupportF1, edgeNormalizedSHD, ...
    triangleNormalizedSHD, edgeRelativeError, triangleRelativeError, ...
    edgeRelativeErrorFinite, triangleRelativeErrorFinite, ...
    priorResidualNorm, smoothnessResidual, edgeCountResidual, ...
    triangleCountResidual);
end
