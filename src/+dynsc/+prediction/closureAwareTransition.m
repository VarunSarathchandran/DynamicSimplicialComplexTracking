function [nextState, details] = closureAwareTransition( ...
    state, boundaryEdges, parameters)
%CLOSUREAWARETRANSITION Evaluate the relaxed conditional-mean transition.
%
%   f([x;y]) = [p; q .* g],
%   p_e = beta1_e + phi1_e*x_e,
%   q_t = beta2_t + phi2_t*y_t,
%   g_t = prod_{e in boundary(t)} p_e.

state = state(:);
nTriangles = size(boundaryEdges, 1);
nEdges = numel(state) - nTriangles;
assert(nEdges > 0, 'dynsc:InvalidClosureAwareStateSize', ...
    'The state must contain an edge block followed by a triangle block.');
validateattributes(state, {'numeric'}, ...
    {'real', 'finite', 'numel', nEdges + nTriangles});
validateattributes(boundaryEdges, {'numeric'}, ...
    {'integer', 'size', [nTriangles, 3], '>=', 1, '<=', nEdges});

parameters = dynsc.transition.resolveParameters( ...
    parameters, nEdges, nTriangles);
x = state(1:nEdges);
y = state(nEdges + (1:nTriangles));

p = parameters.beta1 + parameters.phi1 .* x;
q = parameters.beta2 + parameters.phi2 .* y;
boundaryProbability = reshape( ...
    p(boundaryEdges(:)), size(boundaryEdges));
softBoundaryGate = prod(boundaryProbability, 2);
nextState = [p; q .* softBoundaryGate];

details = struct();
details.edgeProbability = p;
details.triangleProposalProbability = q;
details.softBoundaryGate = softBoundaryGate;
end
