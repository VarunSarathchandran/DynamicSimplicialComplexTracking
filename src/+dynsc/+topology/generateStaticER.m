function topology = generateStaticER(cfg, nTime)
%GENERATESTATICER Generate a static binary SC over all requested time steps.
%
% The construction follows generate_simplicial_complex.m:
%   1. draw a connected Erdos-Renyi graph;
%   2. enumerate its 3-cliques;
%   3. retain a uniform random fraction of them as filled triangles.

if nargin < 2
    nTime = 1;
end

validateattributes(nTime, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});

nNodes = cfg.numNodes;
edgeProbability = cfg.edgeProbability;

fullIncidence = gen_B12(nNodes);
fullB1 = fullIncidence.B1;
fullB2 = fullIncidence.B2;

A = generate_connected_ER(nNodes, edgeProbability);
L0 = diag(sum(A)) - A;

edgeOutput = compute_B1(A);
s1Initial = edgeOutput.w1(:);
B1 = edgeOutput.B1;

triangleOutput = compute_B2(A, B1, s1Initial);
cliqueIndicator = logical(triangleOutput.w2(:));
cliqueIndices = find(cliqueIndicator);
numCliques = numel(cliqueIndices);
numSelectedTriangles = round(cfg.triangleFraction * numCliques);

s2Initial = false(size(cliqueIndicator));
if numSelectedTriangles > 0
    selectedLocalIndices = randperm(numCliques, numSelectedTriangles);
    selectedCliqueIndices = cliqueIndices(selectedLocalIndices);
    s2Initial(selectedCliqueIndices) = true;
end

s1 = repmat(s1Initial, 1, nTime);
s2 = repmat(s2Initial, 1, nTime);
B2Initial = fullB2(:, s2Initial);
LB1 = B1' * B1;
LB2Initial = B2Initial * B2Initial';
L1Initial = LB1 + LB2Initial;

topology = struct();
topology.N = nNodes;
topology.p = edgeProbability;
topology.triangleFraction = cfg.triangleFraction;
topology.numTimeSteps = nTime;
topology.isDynamic = false;
topology.changeTimes = zeros(1, 0);
topology.changeLog = emptyChangeLog();
topology.A = A;
topology.L0 = L0;
topology.s1 = s1;
topology.s2 = s2;
topology.B1 = B1;
topology.B1_full = fullB1;
topology.B2_full = fullB2;
topology.boundaryEdges = dynsc.topology.resolveBoundaryEdges(topology);
topology.B2_initial = B2Initial;
topology.L1_initial = L1Initial;
topology.LB1 = LB1;
topology.LB2_initial = LB2Initial;
topology.edgeNodes = nchoosek(1:nNodes, 2);
topology.triangleNodes = nchoosek(1:nNodes, 3);
topology.feasibleTriangleMask = cliqueIndicator;
topology.num_3cliques = numCliques;
topology.num_edges = sum(s1 > 0, 1);
topology.num_triangles = sum(s2 > 0, 1);
topology.C1_true = sum(s1, 1);
topology.C2_true = sum(s2, 1);
topology.simulated_flag = true;
end

function changeLog = emptyChangeLog()
changeLog = struct( ...
    'time', {}, ...
    'removedTriangleIndex', {}, ...
    'addedTriangleIndex', {}, ...
    'removedTriangleNodes', {}, ...
    'addedTriangleNodes', {});
end
