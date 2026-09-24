function topology = extractTopologySegment(dayTopology, startClock, durationHours)
%EXTRACTTOPOLOGYSEGMENT Extract one fixed clock interval from a daily SC.

validateattributes(durationHours, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
startClock = string(startClock);
assert(isscalar(startClock), 'dynsc:InvalidSegmentStart', ...
    'startClock must be scalar text in HH:mm format.');

parts = sscanf(char(startClock), '%d:%d');
assert(numel(parts) == 2 && parts(1) >= 0 && parts(1) <= 23 && ...
    parts(2) >= 0 && parts(2) <= 59, ...
    'dynsc:InvalidSegmentStart', ...
    'startClock must use a valid HH:mm time.');

dayOrigin = dateshift(dayTopology.windowStart(1), 'start', 'day');
segmentStart = dayOrigin + hours(parts(1)) + minutes(parts(2));
segmentEnd = segmentStart + hours(durationHours);
keep = dayTopology.windowStart >= segmentStart & ...
    dayTopology.windowEnd <= segmentEnd;
indices = find(keep);

expectedStates = durationHours * 60 / dayTopology.windowMinutes;
assert(abs(expectedStates - round(expectedStates)) < 1e-12, ...
    'dynsc:NonIntegralSegmentLength', ...
    'The segment must contain an integer number of topology windows.');
assert(numel(indices) == round(expectedStates), ...
    'dynsc:IncompleteTopologySegment', ...
    ['Requested %d states in %s--%s, but only %d complete windows ' ...
     'were available.'], round(expectedStates), string(segmentStart), ...
    string(segmentEnd), numel(indices));
assert(all(diff(indices) == 1), 'dynsc:NonContiguousTopologySegment', ...
    'The selected topology windows must be contiguous.');

topology = dayTopology;
topology.s1 = dayTopology.s1(:, indices);
topology.s2 = dayTopology.s2(:, indices);
topology.boundaryGate = dayTopology.boundaryGate(:, indices);
topology.C1_true = sum(topology.s1, 1);
topology.C2_true = sum(topology.s2, 1);
topology.num_edges = topology.C1_true;
topology.num_triangles = topology.C2_true;
topology.numFeasibleTriangles = sum(topology.boundaryGate, 1);
topology.numOpenTriangles = sum(topology.boundaryGate & ~topology.s2, 1);
topology.windowStart = dayTopology.windowStart(indices);
topology.windowEnd = dayTopology.windowEnd(indices);
topology.numTimeSteps = numel(indices);

changed = false(1, topology.numTimeSteps);
if topology.numTimeSteps > 1
    changed(2:end) = any(diff(topology.s1, 1, 2) ~= 0, 1) | ...
        any(diff(topology.s2, 1, 2) ~= 0, 1);
end
topology.changeTimes = find(changed);
topology.changeLog = struct([]);

initialEdges = logical(topology.s1(:, 1));
initialTriangles = logical(topology.s2(:, 1));
topology.B1 = topology.B1_full(:, initialEdges);
topology.B2_initial = topology.B2_full(:, initialTriangles);
initialAdjacency = zeros(topology.N);
initialPairs = topology.edgeNodes(initialEdges, :);
for iEdge = 1:size(initialPairs, 1)
    initialAdjacency(initialPairs(iEdge, 1), initialPairs(iEdge, 2)) = 1;
    initialAdjacency(initialPairs(iEdge, 2), initialPairs(iEdge, 1)) = 1;
end
topology.A = initialAdjacency;
topology.L0 = diag(sum(initialAdjacency, 2)) - initialAdjacency;
topology.segmentStart = segmentStart;
topology.segmentEnd = segmentEnd;
topology.sourceWindowIndices = indices(:).';
end
