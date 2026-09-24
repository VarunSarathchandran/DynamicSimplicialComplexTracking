function F = closureAwareJacobian(state, boundaryEdges, parameters)
%CLOSUREAWAREJACOBIAN Sparse Jacobian of closureAwareTransition.

state = state(:);
nTriangles = size(boundaryEdges, 1);
nEdges = numel(state) - nTriangles;
assert(nEdges > 0, 'dynsc:InvalidClosureAwareStateSize', ...
    'The state must contain an edge block followed by a triangle block.');
validateattributes(boundaryEdges, {'numeric'}, ...
    {'integer', 'size', [nTriangles, 3], '>=', 1, '<=', nEdges});

parameters = dynsc.transition.resolveParameters( ...
    parameters, nEdges, nTriangles);
x = state(1:nEdges);
y = state(nEdges + (1:nTriangles));
p = parameters.beta1 + parameters.phi1 .* x;
q = parameters.beta2 + parameters.phi2 .* y;

edge1 = boundaryEdges(:, 1);
edge2 = boundaryEdges(:, 2);
edge3 = boundaryEdges(:, 3);
p1 = p(edge1);
p2 = p(edge2);
p3 = p(edge3);

% Compute the three products directly. In particular, do not divide the
% complete gate by p_e because valid parameter settings can make p_e zero.
value1 = parameters.phi1(edge1) .* q .* p2 .* p3;
value2 = parameters.phi1(edge2) .* q .* p1 .* p3;
value3 = parameters.phi1(edge3) .* q .* p1 .* p2;

triangleRows = repmat((1:nTriangles)', 3, 1);
edgeColumns = [edge1; edge2; edge3];
F21 = sparse(triangleRows, edgeColumns, ...
    [value1; value2; value3], nTriangles, nEdges);

softBoundaryGate = p1 .* p2 .* p3;
F11 = spdiags(parameters.phi1, 0, nEdges, nEdges);
F12 = sparse(nEdges, nTriangles);
F22 = spdiags(parameters.phi2 .* softBoundaryGate, ...
    0, nTriangles, nTriangles);
F = [F11, F12; F21, F22];
end
