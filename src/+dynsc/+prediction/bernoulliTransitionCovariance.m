function [Q, details] = bernoulliTransitionCovariance( ...
    state, boundaryEdges, parameters)
%BERNOULLITRANSITIONCOVARIANCE Exact conditional covariance of one hard step.
%
% The returned sparse matrix is Cov([X_k;Y_k] | state), where edges are
% sampled independently and triangle proposals are independently sampled
% before being multiplied by their hard boundary gates.

state = state(:);
nTriangles = size(boundaryEdges, 1);
nEdges = numel(state) - nTriangles;
assert(nEdges > 0, 'dynsc:InvalidClosureAwareStateSize', ...
    'The state must contain an edge block followed by a triangle block.');
validateattributes(boundaryEdges, {'numeric'}, ...
    {'integer', 'size', [nTriangles, 3], '>=', 1, '<=', nEdges});

[~, transitionDetails] = dynsc.prediction.closureAwareTransition( ...
    state, boundaryEdges, parameters);
p = transitionDetails.edgeProbability;
q = transitionDetails.triangleProposalProbability;

probabilityTolerance = 1e-10;
assert(all(p >= -probabilityTolerance & ...
    p <= 1 + probabilityTolerance) && ...
    all(q >= -probabilityTolerance & ...
    q <= 1 + probabilityTolerance), ...
    'dynsc:InvalidBernoulliProbability', ...
    ['The previous state and transition parameters must produce ' ...
     'probabilities in [0,1].']);
p = min(max(p, 0), 1);
q = min(max(q, 0), 1);

boundaryProbability = reshape( ...
    p(boundaryEdges(:)), size(boundaryEdges));
g = prod(boundaryProbability, 2);
m = q .* g;

edgeVariance = p .* (1 - p);
triangleVariance = m .* (1 - m);
Q11 = spdiags(edgeVariance, 0, nEdges, nEdges);

triangleIndex = repmat((1:nTriangles)', 3, 1);
edgeIndex = boundaryEdges(:);
edgeTriangleCovariance = ...
    m(triangleIndex) .* (1 - p(edgeIndex));
Q12 = sparse(edgeIndex, triangleIndex, ...
    edgeTriangleCovariance, nEdges, nTriangles);

[trianglePair1, trianglePair2, sharedEdge] = ...
    trianglePairsSharingEdges(boundaryEdges, nEdges);
pairCovariance = zeros(numel(sharedEdge), 1);
if ~isempty(sharedEdge)
    otherProduct1 = boundaryProductExcludingSharedEdge( ...
        boundaryEdges(trianglePair1, :), sharedEdge, p);
    otherProduct2 = boundaryProductExcludingSharedEdge( ...
        boundaryEdges(trianglePair2, :), sharedEdge, p);
    pairCovariance = ...
        q(trianglePair1) .* q(trianglePair2) .* ...
        p(sharedEdge) .* (1 - p(sharedEdge)) .* ...
        otherProduct1 .* otherProduct2;
end

Q22 = spdiags(triangleVariance, 0, nTriangles, nTriangles) + ...
    sparse([trianglePair1; trianglePair2], ...
        [trianglePair2; trianglePair1], ...
        [pairCovariance; pairCovariance], ...
        nTriangles, nTriangles);

Q = [Q11, Q12; Q12', Q22];
Q = 0.5 * (Q + Q');

details = struct();
details.conditionalMean = [p; m];
details.edgeProbability = p;
details.triangleProposalProbability = q;
details.softBoundaryGate = g;
details.triangleProbability = m;
details.edgeVariance = edgeVariance;
details.triangleVariance = triangleVariance;
details.numTrianglePairsSharingEdges = numel(sharedEdge);
end

function [triangle1, triangle2, sharedEdge] = ...
    trianglePairsSharingEdges(boundaryEdges, nEdges)
nTriangles = size(boundaryEdges, 1);
triangleIndex = repmat((1:nTriangles)', 3, 1);
edgeIndex = boundaryEdges(:);
incidence = sparse(edgeIndex, triangleIndex, 1, nEdges, nTriangles);

overlap = triu(incidence' * incidence, 1);
assert(isempty(nonzeros(overlap)) || all(nonzeros(overlap) == 1), ...
    'dynsc:InvalidTriangleBoundaryOverlap', ...
    'Distinct candidate triangles may share at most one boundary edge.');

incidentCount = full(sum(incidence, 2));
nPairs = sum(incidentCount .* (incidentCount - 1) / 2);
triangle1 = zeros(nPairs, 1);
triangle2 = zeros(nPairs, 1);
sharedEdge = zeros(nPairs, 1);

nextPair = 1;
for iEdge = 1:nEdges
    incidentTriangles = find(incidence(iEdge, :));
    if numel(incidentTriangles) < 2
        continue;
    end
    pairs = nchoosek(incidentTriangles, 2);
    pairIndices = nextPair:(nextPair + size(pairs, 1) - 1);
    triangle1(pairIndices) = pairs(:, 1);
    triangle2(pairIndices) = pairs(:, 2);
    sharedEdge(pairIndices) = iEdge;
    nextPair = pairIndices(end) + 1;
end
end

function product = boundaryProductExcludingSharedEdge( ...
    pairBoundaries, sharedEdge, p)
nPairs = size(pairBoundaries, 1);
product = ones(nPairs, 1);
for iBoundary = 1:3
    currentEdge = pairBoundaries(:, iBoundary);
    useEdge = currentEdge ~= sharedEdge;
    product(useEdge) = product(useEdge) .* p(currentEdge(useEdge));
end
end
