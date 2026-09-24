function estimate = solveMeasurementOnlyGurobi( ...
    y, H, noiseCovariance, topology, gurobiParameters)
%SOLVEMEASUREMENTONLYGUROBI Solve one static weighted-residual problem.
%
%   minimize_s (y - H*s)' R^{-1} (y - H*s)
%
% subject to 0 <= s <= 1 and s2(t) <= s1(e) for every edge e in the
% boundary of candidate triangle t.

if nargin < 5 || isempty(gurobiParameters)
    gurobiParameters = struct();
    gurobiParameters.outputflag = 0;
end

y = y(:);
nMeasurements = numel(y);
nEdges = size(topology.B1_full, 2);
nTriangles = size(topology.B2_full, 2);
nState = nEdges + nTriangles;

validateattributes(H, {'numeric'}, ...
    {'real', 'finite', 'size', [nMeasurements, nState]});
validateattributes(noiseCovariance, {'numeric'}, ...
    {'real', 'finite', 'size', [nMeasurements, nMeasurements]});

noiseCovariance = (noiseCovariance + noiseCovariance') / 2;
[noiseFactor, noiseFlag] = chol(noiseCovariance, 'lower');
assert(noiseFlag == 0, ...
    'dynsc:NonPositiveDefiniteMeasurementCovariance', ...
    'The measurement-noise covariance must be positive definite.');
noisePrecision = noiseFactor' \ ...
    (noiseFactor \ eye(nMeasurements));

quadraticMatrix = H' * noisePrecision * H;
quadraticMatrix = (quadraticMatrix + quadraticMatrix') / 2;
linearVector = H' * noisePrecision * y;

% Gurobi minimizes x'Qx + obj'x + objcon, with no factor of 1/2.
model = struct();
model.Q = sparse(quadraticMatrix);
model.obj = -2 * linearVector;
model.objcon = y' * noisePrecision * y;
model.modelsense = 'min';
model.lb = zeros(nState, 1);
model.ub = ones(nState, 1);
model.vtype = repmat('C', 1, nState);

idxS1 = 1:nEdges;
idxS2 = nEdges + (1:nTriangles);
[edgeIndex, triangleIndex] = find(abs(topology.B2_full));
nInclusions = numel(edgeIndex);

rows = [(1:nInclusions)'; (1:nInclusions)'];
edgeColumns = idxS1(edgeIndex);
triangleColumns = idxS2(triangleIndex);
columns = [edgeColumns(:); triangleColumns(:)];
values = [-ones(nInclusions, 1); ones(nInclusions, 1)];

model.A = sparse(rows, columns, values, nInclusions, nState);
model.rhs = zeros(nInclusions, 1);
model.sense = repmat('<', nInclusions, 1);

result = gurobi(model, gurobiParameters);
assert(isfield(result, 'x'), 'dynsc:GurobiStaticBaselineFailed', ...
    'Gurobi returned status %s without a state estimate.', result.status);

state = result.x(:);
measurementResidual = y - H * state;
objectiveValue = measurementResidual' * ...
    noisePrecision * measurementResidual;

offDiagonalNoise = noiseCovariance - diag(diag(noiseCovariance));
noiseScale = max(1, norm(noiseCovariance, 'fro'));
assert(norm(offDiagonalNoise, 'fro') <= 1e-12 * noiseScale, ...
    'dynsc:NonDiagonalMeasurementCovariance', ...
    ['The static objective decomposition requires independent virtual ' ...
     'measurement errors.']);

measurementVariances = diag(noiseCovariance);
objectiveTerms = struct();
objectiveTerms.smoothness = ...
    measurementResidual(1)^2 / measurementVariances(1);
objectiveTerms.edgeTotal = ...
    measurementResidual(2)^2 / measurementVariances(2);
objectiveTerms.triangleTotal = ...
    measurementResidual(3)^2 / measurementVariances(3);
objectiveTerms.total = objectiveTerms.smoothness + ...
    objectiveTerms.edgeTotal + objectiveTerms.triangleTotal;

estimate = struct();
estimate.state = state;
estimate.measurementResidual = measurementResidual;
estimate.objectiveValue = objectiveValue;
estimate.objectiveTerms = objectiveTerms;
estimate.gurobiResult = result;
estimate.gurobiModel = model;
end
