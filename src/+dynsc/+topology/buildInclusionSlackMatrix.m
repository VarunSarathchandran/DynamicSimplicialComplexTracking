function A = buildInclusionSlackMatrix(boundaryEdges, nEdges, nTriangles)
%BUILDINCLUSIONSLACKMATRIX Construct all edge-minus-triangle slacks.
%
% Each row corresponds to one pair (e,t) with e on the boundary of t and
% satisfies
%
%   (A*s)_(e,t) = s1(e) - s2(t).

validateattributes(boundaryEdges, {'numeric'}, ...
    {'integer', 'size', [nTriangles, 3], '>=', 1, '<=', nEdges});

nInclusions = 3 * nTriangles;
triangleIndex = repelem((1:nTriangles)', 3, 1);
edgeIndex = reshape(boundaryEdges', [], 1);
rows = [(1:nInclusions)'; (1:nInclusions)'];
columns = [edgeIndex; nEdges + triangleIndex];
values = [ones(nInclusions, 1); -ones(nInclusions, 1)];
A = sparse(rows, columns, values, ...
    nInclusions, nEdges + nTriangles);
end
