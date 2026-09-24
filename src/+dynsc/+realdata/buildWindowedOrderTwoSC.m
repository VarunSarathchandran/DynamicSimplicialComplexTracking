function [topology, diagnostics] = buildWindowedOrderTwoSC( ...
    eventTime, simplices, selection, windowMinutes, dayIndex)
%BUILDWINDOWEDORDERTWOSC Aggregate 20-second simplices into SC snapshots.

validateattributes(windowMinutes, {'numeric'}, ...
    {'scalar', 'positive', 'finite'});
assert(~isempty(eventTime) && numel(eventTime) == numel(simplices), ...
    'dynsc:InvalidDayEvents', ...
    'Each simplex must have one event timestamp.');

nNodes = height(selection.nodes);
edgeNodes = nchoosek(1:nNodes, 2);
triangleNodes = nchoosek(1:nNodes, 3);
nEdges = size(edgeNodes, 1);
nTriangles = size(triangleNodes, 1);

edgeLookup = zeros(nNodes, nNodes);
for iEdge = 1:nEdges
    edgeLookup(edgeNodes(iEdge, 1), edgeNodes(iEdge, 2)) = iEdge;
    edgeLookup(edgeNodes(iEdge, 2), edgeNodes(iEdge, 1)) = iEdge;
end

triangleLookup = zeros(nNodes, nNodes, nNodes);
for iTriangle = 1:nTriangles
    nodes = triangleNodes(iTriangle, :);
    permutations = perms(nodes);
    for iPermutation = 1:size(permutations, 1)
        triangleLookup( ...
            permutations(iPermutation, 1), ...
            permutations(iPermutation, 2), ...
            permutations(iPermutation, 3)) = iTriangle;
    end
end

maxAhornId = max(selection.nodes.AhornId);
sourceToLocal = zeros(maxAhornId, 1);
sourceToLocal(selection.nodes.AhornId) = selection.nodes.LocalId;

dayStart = dateshift(min(eventTime), 'start', 'day');
firstMinute = hour(min(eventTime)) * 60 + minute(min(eventTime));
alignedMinute = windowMinutes * floor(firstMinute / windowMinutes);
windowOrigin = dayStart + minutes(alignedMinute);
windowDuration = minutes(windowMinutes);
nWindows = floor(seconds(max(eventTime) - windowOrigin) / ...
    seconds(windowDuration)) + 1;

windowStart = windowOrigin + (0:nWindows-1).' * windowDuration;
windowEnd = windowStart + windowDuration;
s1 = false(nEdges, nWindows);
s2 = false(nTriangles, nWindows);
numRelevantSimplexEvents = zeros(nWindows, 1);

for iEvent = 1:numel(simplices)
    sourceNodes = simplices{iEvent};
    sourceNodes = sourceNodes(sourceNodes <= maxAhornId);
    localNodes = sourceToLocal(sourceNodes);
    localNodes = unique(localNodes(localNodes > 0), 'stable');
    if numel(localNodes) < 2
        continue;
    end

    iWindow = floor(seconds(eventTime(iEvent) - windowOrigin) / ...
        seconds(windowDuration)) + 1;
    assert(iWindow >= 1 && iWindow <= nWindows, ...
        'dynsc:InvalidPrimarySchoolWindowIndex', ...
        'An event was assigned outside the daily window range.');
    numRelevantSimplexEvents(iWindow) = ...
        numRelevantSimplexEvents(iWindow) + 1;

    pairs = nchoosek(localNodes, 2);
    pairIndices = arrayfun( ...
        @(i) edgeLookup(pairs(i, 1), pairs(i, 2)), ...
        (1:size(pairs, 1)).');
    s1(pairIndices, iWindow) = true;

    if numel(localNodes) >= 3
        triples = nchoosek(localNodes, 3);
        triangleIndices = arrayfun( ...
            @(i) triangleLookup( ...
                triples(i, 1), triples(i, 2), triples(i, 3)), ...
            (1:size(triples, 1)).');
        s2(triangleIndices, iWindow) = true;
    end
end

boundaryEdges = zeros(nTriangles, 3);
for iTriangle = 1:nTriangles
    nodes = triangleNodes(iTriangle, :);
    boundaryEdges(iTriangle, :) = [ ...
        edgeLookup(nodes(1), nodes(2)), ...
        edgeLookup(nodes(1), nodes(3)), ...
        edgeLookup(nodes(2), nodes(3))];
end

boundaryGate = false(nTriangles, nWindows);
for iWindow = 1:nWindows
    boundaryState = reshape( ...
        s1(boundaryEdges(:), iWindow), size(boundaryEdges));
    boundaryGate(:, iWindow) = all(boundaryState, 2);
end
assert(all(~s2(:) | boundaryGate(:)), ...
    'dynsc:EmpiricalInclusionViolation', ...
    'A constructed triangle is missing at least one boundary edge.');

fullIncidence = gen_B12(nNodes);
assert(size(fullIncidence.B1, 2) == nEdges && ...
    size(fullIncidence.B2, 2) == nTriangles, ...
    'dynsc:CandidateOrderingMismatch', ...
    'Candidate incidence matrices have unexpected dimensions.');
for iTriangle = 1:nTriangles
    assert(isequal(sort(find(fullIncidence.B2(:, iTriangle))).', ...
        sort(boundaryEdges(iTriangle, :))), ...
        'dynsc:CandidateOrderingMismatch', ...
        'Triangle ordering does not agree with the incidence matrix.');
end

initialEdges = s1(:, 1);
initialTriangles = s2(:, 1);
initialAdjacency = zeros(nNodes);
initialPairs = edgeNodes(initialEdges, :);
for iEdge = 1:size(initialPairs, 1)
    initialAdjacency(initialPairs(iEdge, 1), initialPairs(iEdge, 2)) = 1;
    initialAdjacency(initialPairs(iEdge, 2), initialPairs(iEdge, 1)) = 1;
end

changed = false(1, nWindows);
if nWindows > 1
    changed(2:end) = any(diff(s1, 1, 2) ~= 0, 1) | ...
        any(diff(s2, 1, 2) ~= 0, 1);
end

topology = struct();
topology.N = nNodes;
topology.numTimeSteps = nWindows;
topology.isDynamic = true;
topology.dynamicsType = "empirical_primary_school";
topology.dayIndex = dayIndex;
topology.windowStart = windowStart;
topology.windowEnd = windowEnd;
topology.windowMinutes = windowMinutes;
topology.s1 = s1;
topology.s2 = s2;
topology.C1_true = sum(s1, 1);
topology.C2_true = sum(s2, 1);
topology.num_edges = topology.C1_true;
topology.num_triangles = topology.C2_true;
topology.boundaryGate = boundaryGate;
topology.numFeasibleTriangles = sum(boundaryGate, 1);
topology.numOpenTriangles = sum(boundaryGate & ~s2, 1);
topology.edgeNodes = edgeNodes;
topology.triangleNodes = triangleNodes;
topology.boundaryEdges = boundaryEdges;
topology.B1_full = fullIncidence.B1;
topology.B2_full = fullIncidence.B2;
topology.B1 = topology.B1_full(:, initialEdges);
topology.B2_initial = topology.B2_full(:, initialTriangles);
topology.A = initialAdjacency;
topology.L0 = diag(sum(initialAdjacency, 2)) - initialAdjacency;
topology.changeTimes = find(changed);
topology.changeLog = struct([]);
topology.simulated_flag = false;
topology.selectedNodes = selection.nodes;

diagnostics = table( ...
    repmat(dayIndex, nWindows, 1), (1:nWindows).', ...
    windowStart, windowEnd, topology.C1_true.', topology.C2_true.', ...
    topology.numFeasibleTriangles.', topology.numOpenTriangles.', ...
    numRelevantSimplexEvents, ...
    'VariableNames', { ...
        'Day', 'TimeIndex', 'WindowStart', 'WindowEnd', ...
        'NumEdges', 'NumTriangles', 'NumFeasibleTriangles', ...
        'NumOpenTriangles', 'NumRelevantSimplexEvents'});
end
