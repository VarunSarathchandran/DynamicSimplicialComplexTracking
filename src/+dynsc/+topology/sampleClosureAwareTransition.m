function [newEdges, newTriangles, details] = ...
    sampleClosureAwareTransition( ...
        previousEdges, previousTriangles, boundaryEdges, parameters)
%SAMPLECLOSUREAWARETRANSITION Sample one hard binary Markov transition.

previousEdges = previousEdges(:);
previousTriangles = previousTriangles(:);
nEdges = numel(previousEdges);
nTriangles = numel(previousTriangles);

assert(all(previousEdges == 0 | previousEdges == 1) && ...
    all(previousTriangles == 0 | previousTriangles == 1), ...
    'dynsc:NonBinaryPreviousTopology', ...
    'The hard topology transition requires binary previous states.');
validateattributes(boundaryEdges, {'numeric'}, ...
    {'integer', 'size', [nTriangles, 3], '>=', 1, '<=', nEdges});

parameters = dynsc.transition.resolveParameters( ...
    parameters, nEdges, nTriangles);

edgeProbability = parameters.beta1 + ...
    parameters.phi1 .* double(previousEdges);
newEdges = rand(nEdges, 1) < edgeProbability;

triangleProposalProbability = parameters.beta2 + ...
    parameters.phi2 .* double(previousTriangles);
triangleProposal = rand(nTriangles, 1) < triangleProposalProbability;
boundaryState = reshape( ...
    newEdges(boundaryEdges(:)), size(boundaryEdges));
boundaryGate = all(boundaryState, 2);
newTriangles = boundaryGate & triangleProposal;

details = struct();
details.edgeProbability = edgeProbability;
details.triangleProposalProbability = triangleProposalProbability;
details.triangleProposal = triangleProposal;
details.boundaryGate = boundaryGate;
end
