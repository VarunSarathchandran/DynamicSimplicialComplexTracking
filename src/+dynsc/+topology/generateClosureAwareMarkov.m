function topology = generateClosureAwareMarkov(cfg, nTime)
%GENERATECLOSUREAWAREMARKOV Generate a binary closure-aware dynamic SC.
%
% Time index one is the valid ER-based SC produced by generateStaticER.
% Subsequent states use independent Bernoulli edge transitions followed by
% Bernoulli triangle proposals gated by the newly sampled edge boundaries.

validateattributes(nTime, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
assert(isfield(cfg, 'transition'), 'dynsc:MissingTransitionParameters', ...
    'Closure-aware generation requires cfg.transition.');

topology = dynsc.topology.generateStaticER(cfg, 1);
boundaryEdges = topology.boundaryEdges;
nEdges = size(topology.B1_full, 2);
nTriangles = size(topology.B2_full, 2);
parameters = dynsc.transition.resolveParameters( ...
    cfg.transition, nEdges, nTriangles);

s1 = false(nEdges, nTime);
s2 = false(nTriangles, nTime);
s1(:, 1) = logical(topology.s1(:, 1));
s2(:, 1) = logical(topology.s2(:, 1));

edgeProbability = nan(nEdges, nTime);
triangleProposalProbability = nan(nTriangles, nTime);
triangleProposal = false(nTriangles, nTime);
boundaryGate = false(nTriangles, nTime);
initialBoundaryState = reshape( ...
    s1(boundaryEdges(:), 1), size(boundaryEdges));
boundaryGate(:, 1) = all(initialBoundaryState, 2);

for iTime = 2:nTime
    [s1(:, iTime), s2(:, iTime), step] = ...
        dynsc.topology.sampleClosureAwareTransition( ...
            s1(:, iTime - 1), s2(:, iTime - 1), ...
            boundaryEdges, parameters);
    edgeProbability(:, iTime) = step.edgeProbability;
    triangleProposalProbability(:, iTime) = ...
        step.triangleProposalProbability;
    triangleProposal(:, iTime) = step.triangleProposal;
    boundaryGate(:, iTime) = step.boundaryGate;
end

changed = false(1, nTime);
if nTime > 1
    changed(2:end) = any(diff(s1, 1, 2) ~= 0, 1) | ...
        any(diff(s2, 1, 2) ~= 0, 1);
end
changeTimes = find(changed);
changeLog = buildChangeLog(s1, s2, changeTimes);

topology.numTimeSteps = nTime;
topology.isDynamic = true;
topology.dynamicsType = "closure_aware_markov";
topology.transitionParameters = parameters;
topology.changeTimes = changeTimes;
topology.changeLog = changeLog;
topology.s1 = s1;
topology.s2 = s2;
topology.C1_true = sum(s1, 1);
topology.C2_true = sum(s2, 1);
topology.num_edges = topology.C1_true;
topology.num_triangles = topology.C2_true;
topology.edgeTransitionProbability = edgeProbability;
topology.triangleProposalProbability = triangleProposalProbability;
topology.triangleProposal = triangleProposal;
topology.boundaryGate = boundaryGate;
end

function changeLog = buildChangeLog(s1, s2, changeTimes)
nChanges = numel(changeTimes);
changeLog = repmat(struct( ...
    'time', [], ...
    'addedEdgeIndices', [], ...
    'removedEdgeIndices', [], ...
    'addedTriangleIndices', [], ...
    'removedTriangleIndices', []), 1, nChanges);

for iChange = 1:nChanges
    iTime = changeTimes(iChange);
    previousEdges = s1(:, iTime - 1);
    currentEdges = s1(:, iTime);
    previousTriangles = s2(:, iTime - 1);
    currentTriangles = s2(:, iTime);

    changeLog(iChange).time = iTime;
    changeLog(iChange).addedEdgeIndices = ...
        find(~previousEdges & currentEdges);
    changeLog(iChange).removedEdgeIndices = ...
        find(previousEdges & ~currentEdges);
    changeLog(iChange).addedTriangleIndices = ...
        find(~previousTriangles & currentTriangles);
    changeLog(iChange).removedTriangleIndices = ...
        find(previousTriangles & ~currentTriangles);
end
end
