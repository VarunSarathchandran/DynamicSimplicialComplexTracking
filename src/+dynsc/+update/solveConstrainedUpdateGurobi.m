function update = solveConstrainedUpdateGurobi( ...
    predictedState, predictedCovariance, y, H, noiseCovariance, ...
    topology, gurobiParameters)
%SOLVECONSTRAINEDUPDATEGUROBI Constrained quadratic state update.
%
% The optimized objective is
%
%   (s - sPred)' / PPred * (s - sPred) ...
%       + (y - H*s)' / R * (y - H*s),
%
% subject to 0 <= s <= 1 and s2(t) <= s1(e) for every edge e on
% the boundary of triangle t.

if nargin < 7 || isempty(gurobiParameters)
    gurobiParameters = struct();
    gurobiParameters.outputflag = 0;
end

predictedState = predictedState(:);
y = y(:);
nState = numel(predictedState);
nMeasurements = numel(y);
nEdges = size(topology.B1_full, 2);
nTriangles = size(topology.B2_full, 2);

assert(nState == nEdges + nTriangles, ...
    'dynsc:UpdateStateSizeMismatch', ...
    'The state must contain %d edge and %d triangle entries.', ...
    nEdges, nTriangles);
validateattributes(predictedCovariance, {'numeric'}, ...
    {'real', 'finite', 'size', [nState, nState]});
validateattributes(H, {'numeric'}, ...
    {'real', 'finite', 'size', [nMeasurements, nState]});
validateattributes(noiseCovariance, {'numeric'}, ...
    {'real', 'finite', 'size', [nMeasurements, nMeasurements]});

predictedCovariance = (predictedCovariance + predictedCovariance') / 2;
noiseCovariance = (noiseCovariance + noiseCovariance') / 2;

[predictedFactor, predictedFlag] = chol( ...
    predictedCovariance, 'lower');
assert(predictedFlag == 0, 'dynsc:NonPositiveDefinitePredictionCovariance', ...
    'The predicted covariance must be symmetric positive definite.');

[noiseFactor, noiseFlag] = chol(noiseCovariance, 'lower');
assert(noiseFlag == 0, 'dynsc:NonPositiveDefiniteMeasurementCovariance', ...
    'The measurement-noise covariance must be symmetric positive definite.');

priorPrecision = predictedFactor' \ ...
    (predictedFactor \ eye(nState));
noisePrecision = noiseFactor' \ ...
    (noiseFactor \ eye(nMeasurements));

informationMatrix = priorPrecision + H' * noisePrecision * H;
informationMatrix = (informationMatrix + informationMatrix') / 2;
informationVector = priorPrecision * predictedState + ...
    H' * noisePrecision * y;

% Gurobi's MATLAB model uses x'Qx + obj'x + objcon. Therefore Q is the
% information matrix itself and the linear coefficient is -2 times the
% information vector.
model = struct();
model.Q = sparse(informationMatrix);
model.obj = -2 * informationVector;
model.objcon = predictedState' * priorPrecision * predictedState + ...
    y' * noisePrecision * y;
model.modelsense = 'min';
model.lb = zeros(nState, 1);
model.ub = ones(nState, 1);
model.vtype = repmat('C', 1, nState);

idx_s1 = 1:nEdges;
idx_s2 = nEdges + (1:nTriangles);
B2plus = abs(topology.B2_full);
[e_idx, t_idx] = find(B2plus);
m = numel(e_idx);

r = [(1:m)'; (1:m)'];
edgeColumns = idx_s1(e_idx);
triangleColumns = idx_s2(t_idx);
c = [edgeColumns(:); triangleColumns(:)];
v = [-ones(m, 1); ones(m, 1)];

model.A = sparse(r, c, v, m, nState);
model.rhs = zeros(m, 1);
model.sense = repmat('<', m, 1);

result = gurobi(model, gurobiParameters);
assert(isfield(result, 'x'), 'dynsc:GurobiUpdateFailed', ...
    'Gurobi returned status %s without a state estimate.', result.status);

updatedState = result.x(:);
priorResidual = updatedState - predictedState;
measurementResidual = y - H * updatedState;
objectiveValue = priorResidual' * priorPrecision * priorResidual + ...
    measurementResidual' * noisePrecision * measurementResidual;

offDiagonalNoise = noiseCovariance - ...
    diag(diag(noiseCovariance));
noiseScale = max(1, norm(noiseCovariance, 'fro'));
assert(norm(offDiagonalNoise, 'fro') <= 1e-12 * noiseScale, ...
    'dynsc:NonDiagonalMeasurementCovariance', ...
    ['The four-term objective decomposition requires independent virtual ' ...
     'measurement errors, so the noise covariance must be diagonal.']);

measurementVariances = diag(noiseCovariance);
objectiveTerms = struct();
objectiveTerms.prior = ...
    priorResidual' * priorPrecision * priorResidual;
objectiveTerms.smoothness = ...
    measurementResidual(1)^2 / measurementVariances(1);
objectiveTerms.edgeTotal = ...
    measurementResidual(2)^2 / measurementVariances(2);
objectiveTerms.triangleTotal = ...
    measurementResidual(3)^2 / measurementVariances(3);
objectiveTerms.measurementTotal = objectiveTerms.smoothness + ...
    objectiveTerms.edgeTotal + objectiveTerms.triangleTotal;
objectiveTerms.total = objectiveTerms.prior + ...
    objectiveTerms.measurementTotal;

[informationFactor, informationFlag] = chol( ...
    informationMatrix, 'lower');
assert(informationFlag == 0, 'dynsc:NonPositiveDefiniteInformationMatrix', ...
    'The posterior information matrix must be positive definite.');
updatedCovariance = informationFactor' \ ...
    (informationFactor \ eye(nState));
updatedCovariance = (updatedCovariance + updatedCovariance') / 2;

update = struct();
update.state = updatedState;
update.covariance = updatedCovariance;
update.informationMatrix = informationMatrix;
update.priorResidual = priorResidual;
update.measurementResidual = measurementResidual;
update.objectiveValue = objectiveValue;
update.objectiveTerms = objectiveTerms;
update.gurobiResult = result;
update.gurobiModel = model;
end
