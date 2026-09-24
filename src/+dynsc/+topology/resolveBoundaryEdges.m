function boundaryEdges = resolveBoundaryEdges(topology)
%RESOLVEBOUNDARYEDGES Return one N2-by-3 triangle-boundary index matrix.
%
% New topology generators store boundaryEdges directly. The B2 fallback
% keeps older topology structures and standalone solver calls compatible.

if isfield(topology, 'boundaryEdges') && ...
        ~isempty(topology.boundaryEdges)
    boundaryEdges = double(topology.boundaryEdges);
else
    assert(isfield(topology, 'B2_full'), ...
        'dynsc:MissingTriangleBoundaryData', ...
        'The topology must contain boundaryEdges or B2_full.');
    boundaryEdges = fromBoundaryMatrix(topology.B2_full);
end

nEdges = size(topology.B2_full, 1);
nTriangles = size(topology.B2_full, 2);
validateattributes(boundaryEdges, {'numeric'}, ...
    {'real', 'finite', 'integer', 'size', [nTriangles, 3], ...
     '>=', 1, '<=', nEdges});

assert(all(arrayfun(@(iRow) ...
    numel(unique(boundaryEdges(iRow, :))) == 3, ...
    (1:nTriangles)')), 'dynsc:InvalidTriangleBoundary', ...
    'Every candidate triangle must have three distinct boundary edges.');
end

function boundaryEdges = fromBoundaryMatrix(B2)
nTriangles = size(B2, 2);
[edgeIndex, triangleIndex] = find(abs(B2) > 0);
counts = accumarray(triangleIndex, 1, [nTriangles, 1]);
assert(all(counts == 3), 'dynsc:InvalidTriangleBoundary', ...
    'Every column of B2_full must contain exactly three boundary edges.');

indexedEdges = sortrows([triangleIndex, edgeIndex], [1, 2]);
boundaryEdges = reshape(indexedEdges(:, 2), 3, nTriangles)';
end
