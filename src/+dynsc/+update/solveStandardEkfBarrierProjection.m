function update = solveStandardEkfBarrierProjection( ...
    predictedState, predictedCovariance, y, H, noiseCovariance, ...
    topology, barrierConfig, gurobiParameters)
%SOLVESTANDARDEKFBARRIERPROJECTION EKF barriers plus mean projection.
%
% The nonlinear virtual measurements are
%
%   0 = h_inc(s) + v_inc,
%   0 = h_box(s) + v_box,
%
% where h_inc(s)^2 = -sum(log(A_inc*s)) and
% h_box(s)^2 = -sum(log(s)+log(1-s)). They are linearized once at the
% predicted mean. A standard augmented Kalman update produces
% the unconstrained mean and covariance. The mean is then projected onto the
% closed SC feasible set in the PPost^{-1} metric; PPost itself is retained.

if nargin < 8 || isempty(gurobiParameters)
    gurobiParameters = struct();
    gurobiParameters.outputflag = 0;
end
predictedState = predictedState(:);
y = y(:);
linearizationState = predictedState;
nState = numel(predictedState);
nMeasurements = numel(y);
nEdges = size(topology.B1_full, 2);
nTriangles = size(topology.B2_full, 2);

assert(nState == nEdges + nTriangles, ...
    'dynsc:UpdateStateSizeMismatch', ...
    'The state must contain %d edge and %d triangle entries.', ...
    nEdges, nTriangles);
assert(isfield(topology, 'boundaryEdges'), ...
    'dynsc:MissingBoundaryEdges', ...
    'The standard-EKF barrier update requires topology.boundaryEdges.');
validateattributes(predictedState, {'numeric'}, ...
    {'column', 'real', 'finite'});
validateattributes(predictedCovariance, {'numeric'}, ...
    {'real', 'finite', 'size', [nState, nState]});
validateattributes(H, {'numeric'}, ...
    {'real', 'finite', 'size', [nMeasurements, nState]});
validateattributes(noiseCovariance, {'numeric'}, ...
    {'real', 'finite', 'size', [nMeasurements, nMeasurements]});
validateBarrierConfig(barrierConfig);

lambdaInc = barrierConfig.lambda_inc;
lambdaBox = barrierConfig.lambda_box;
sigmaSquaredInc = 1 / (2 * lambdaInc);
sigmaSquaredBox = 1 / (2 * lambdaBox);

predictedCovariance = 0.5 * ...
    (predictedCovariance + predictedCovariance');
noiseCovariance = 0.5 * ...
    (noiseCovariance + noiseCovariance');
assertPositiveDefinite(predictedCovariance, ...
    'dynsc:NonPositiveDefinitePredictionCovariance', ...
    'The predicted covariance must be symmetric positive definite.');
assertPositiveDefinite(noiseCovariance, ...
    'dynsc:NonPositiveDefiniteMeasurementCovariance', ...
    'The measurement-noise covariance must be symmetric positive definite.');

inclusionSlackMatrix = ...
    dynsc.topology.buildInclusionSlackMatrix( ...
        topology.boundaryEdges, nEdges, nTriangles);
linearizationState = makeInteriorLinearizationState( ...
    linearizationState, topology.boundaryEdges, nEdges, nTriangles, ...
    inclusionSlackMatrix, barrierConfig.interiorMargin);
barrierMeasurement = dynsc.update.barrierMeasurements( ...
    linearizationState, inclusionSlackMatrix);

jInc = barrierMeasurement.inclusionJacobian;
jBox = barrierMeasurement.boxJacobian;
pseudoTargetInc = -barrierMeasurement.inclusionValue + ...
    jInc' * linearizationState;
pseudoTargetBox = -barrierMeasurement.boxValue + ...
    jBox' * linearizationState;

augmentedH = [H; jInc'; jBox'];
augmentedY = [y; pseudoTargetInc; pseudoTargetBox];
augmentedNoiseCovariance = blkdiag( ...
    noiseCovariance, sigmaSquaredInc, sigmaSquaredBox);

innovation = augmentedY - augmentedH * predictedState;
innovationCovariance = augmentedH * predictedCovariance * ...
    augmentedH' + augmentedNoiseCovariance;
innovationCovariance = 0.5 * ...
    (innovationCovariance + innovationCovariance');
assertPositiveDefinite(innovationCovariance, ...
    'dynsc:NonPositiveDefiniteInnovationCovariance', ...
    'The augmented innovation covariance must be positive definite.');

kalmanGain = (predictedCovariance * augmentedH') / ...
    innovationCovariance;
unprojectedState = predictedState + kalmanGain * innovation;

identity = eye(nState);
josephFactor = identity - kalmanGain * augmentedH;
updatedCovariance = josephFactor * predictedCovariance * ...
    josephFactor' + kalmanGain * augmentedNoiseCovariance * ...
    kalmanGain';
updatedCovariance = 0.5 * ...
    (updatedCovariance + updatedCovariance');
assertNumericallyPsd(updatedCovariance);

[updatedState, projectionResult, projectionModel, ...
    posteriorPrecision] = covarianceWeightedProjection( ...
        unprojectedState, updatedCovariance, inclusionSlackMatrix, ...
        gurobiParameters);

priorPrecision = symmetricInverse(predictedCovariance, ...
    'dynsc:NonPositiveDefinitePredictionCovariance');
noisePrecision = symmetricInverse(noiseCovariance, ...
    'dynsc:NonPositiveDefiniteMeasurementCovariance');
linearMeasurementInformation = H' * noisePrecision * H;
LambdaNormal = priorPrecision + linearMeasurementInformation;
LambdaNormal = 0.5 * (LambdaNormal + LambdaNormal');
inclusionInformation = (1 / sigmaSquaredInc) * (jInc * jInc');
boxInformation = (1 / sigmaSquaredBox) * (jBox * jBox');
informationMatrix = LambdaNormal + ...
    inclusionInformation + boxInformation;
informationMatrix = 0.5 * ...
    (informationMatrix + informationMatrix');

% The Joseph covariance and information-form covariance describe the same
% one-shot augmented EKF update in exact arithmetic. Retain the Joseph form
% as the algorithmic covariance because it is numerically stable. Near a
% barrier boundary, the Jacobians and hence informationMatrix can be very
% ill-conditioned, so disagreement with its explicit inverse is a
% diagnostic rather than a reason to stop the filter.
informationCovariance = symmetricInverse(informationMatrix, ...
    'dynsc:NonPositiveDefiniteInformationMatrix');
covarianceConsistencyError = norm( ...
    updatedCovariance - informationCovariance, 'fro');
covarianceConsistencyScale = max([ ...
    norm(updatedCovariance, 'fro'), ...
    norm(informationCovariance, 'fro'), realmin]);
covarianceConsistencyRelativeError = ...
    covarianceConsistencyError / covarianceConsistencyScale;
informationReciprocalCondition = rcond(full(informationMatrix));
informationIdentityResidual = norm( ...
    informationMatrix * updatedCovariance - identity, 'fro') / ...
    max(1, norm(informationMatrix, 'fro') * ...
        norm(updatedCovariance, 'fro'));

priorResidual = updatedState - predictedState;
measurementResidual = y - H * updatedState;
linearizedIncResidual = pseudoTargetInc - jInc' * updatedState;
linearizedBoxResidual = pseudoTargetBox - jBox' * updatedState;
objectiveTerms = decomposeLinearizedObjective( ...
    priorResidual, priorPrecision, measurementResidual, ...
    noiseCovariance, linearizedIncResidual, linearizedBoxResidual, ...
    sigmaSquaredInc, sigmaSquaredBox);

projectionResidual = updatedState - unprojectedState;
projectionObjective = projectionResidual' * ...
    posteriorPrecision * projectionResidual;

update = struct();
update.method = "barrier_standard_ekf_projection";
update.solverName = "Kalman + Gurobi projection";
update.solverStatus = string(projectionResult.status);
update.state = updatedState;
update.unprojectedState = unprojectedState;
update.linearizationState = linearizationState;
update.covariance = updatedCovariance;
update.priorResidual = priorResidual;
update.measurementResidual = measurementResidual;
update.linearizedBarrierResidual = [ ...
    linearizedIncResidual; linearizedBoxResidual];
update.objectiveValue = objectiveTerms.total;
update.objectiveTerms = objectiveTerms;
update.priorPrecision = priorPrecision;
update.linearMeasurementInformation = linearMeasurementInformation;
update.LambdaNormal = LambdaNormal;
update.informationMatrix = informationMatrix;
update.informationTerms.inclusion = inclusionInformation;
update.informationTerms.box = boxInformation;
update.informationTerms.totalBarrier = ...
    inclusionInformation + boxInformation;
update.barrierMeasurement = barrierMeasurement;
update.inclusionSlackMatrix = inclusionSlackMatrix;
update.augmentedMeasurementMatrix = augmentedH;
update.augmentedMeasurementVector = augmentedY;
update.augmentedNoiseCovariance = augmentedNoiseCovariance;
update.innovation = innovation;
update.innovationCovariance = innovationCovariance;
update.kalmanGain = kalmanGain;
update.covarianceConsistencyError = covarianceConsistencyError;
update.covarianceConsistencyRelativeError = ...
    covarianceConsistencyRelativeError;
update.informationReciprocalCondition = ...
    informationReciprocalCondition;
update.informationIdentityResidual = informationIdentityResidual;
update.impliedNoiseVariance.inclusion = sigmaSquaredInc;
update.impliedNoiseVariance.box = sigmaSquaredBox;
update.projectionObjective = projectionObjective;
update.projectionPrecision = posteriorPrecision;
update.gurobiResult = projectionResult;
update.gurobiModel = projectionModel;
end

function validateBarrierConfig(cfg)
required = {'lambda_inc', 'lambda_box', 'interiorMargin'};
for iField = 1:numel(required)
    assert(isfield(cfg, required{iField}), ...
        'dynsc:MissingBarrierUpdateConfig', ...
        'Missing cfg.update.barrier.%s.', required{iField});
end
validateattributes(cfg.lambda_inc, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.lambda_box, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.interiorMargin, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', '<', 1e-3});
end

function state = makeInteriorLinearizationState( ...
    state, boundaryEdges, nEdges, nTriangles, inclusionSlackMatrix, ...
    interiorMargin)
isInterior = all(state > interiorMargin) && ...
    all(state < 1 - interiorMargin) && ...
    all(inclusionSlackMatrix * state > interiorMargin);
if isInterior
    return;
end

% A projected estimate may lie on the boundary, where log barriers are not
% differentiable. Move it only by a numerical amount toward a fixed strict
% SC interior point before evaluating the next Jacobian.
candidate = min(max(state, 0), 1);
for iTriangle = 1:nTriangles
    triangleIndex = nEdges + iTriangle;
    candidate(triangleIndex) = min(candidate(triangleIndex), ...
        min(candidate(boundaryEdges(iTriangle, :))));
end
reference = [ ...
    0.75 * ones(nEdges, 1)
    0.25 * ones(nTriangles, 1)
];
mixingWeight = max(5 * interiorMargin, sqrt(eps));
state = (1 - mixingWeight) * candidate + mixingWeight * reference;

assert(all(state > 0 & state < 1) && ...
    all(inclusionSlackMatrix * state > 0), ...
    'dynsc:NoStrictBarrierLinearizationState', ...
    'Could not construct a strict-interior barrier linearization state.');
end

function [projectedState, result, model, precision] = ...
        covarianceWeightedProjection( ...
            unprojectedState, covariance, inclusionSlackMatrix, ...
            gurobiParameters)
nState = numel(unprojectedState);
precision = symmetricInverse(covariance, ...
    'dynsc:NonPositiveDefinitePosteriorCovariance');
if ~isfield(gurobiParameters, 'FeasibilityTol') && ...
        ~isfield(gurobiParameters, 'feasibilitytol')
    gurobiParameters.FeasibilityTol = 1e-9;
end

% Gurobi minimizes x'Qx + obj'x. These coefficients therefore implement
% (x-unprojectedState)'*precision*(x-unprojectedState), up to its constant.
model = struct();
model.Q = sparse(precision);
model.obj = -2 * precision * unprojectedState;
model.objcon = unprojectedState' * precision * unprojectedState;
model.modelsense = 'min';
model.lb = zeros(nState, 1);
model.ub = ones(nState, 1);
model.vtype = repmat('C', 1, nState);
model.A = sparse(-inclusionSlackMatrix);
model.rhs = zeros(size(inclusionSlackMatrix, 1), 1);
model.sense = repmat('<', size(inclusionSlackMatrix, 1), 1);

result = gurobi(model, gurobiParameters);
assert(isfield(result, 'x'), 'dynsc:GurobiProjectionFailed', ...
    'Gurobi returned status %s without a projected state.', result.status);
projectedState = result.x(:);

feasibilityTolerance = 1e-8;
assert(all(projectedState >= -feasibilityTolerance & ...
    projectedState <= 1 + feasibilityTolerance), ...
    'dynsc:ProjectedStateOutsideBox', ...
    'The covariance-weighted projection violated the unit box.');
assert(all(inclusionSlackMatrix * projectedState >= ...
    -feasibilityTolerance), 'dynsc:ProjectedStateViolatesInclusion', ...
    'The covariance-weighted projection violated simplex inclusion.');
end

function inverse = symmetricInverse(matrix, errorId)
matrix = 0.5 * (matrix + matrix');
[factor, flag] = chol(matrix, 'lower');
assert(flag == 0, errorId, ...
    'The matrix must be symmetric positive definite.');
inverse = factor' \ (factor \ eye(size(matrix, 1)));
inverse = 0.5 * (inverse + inverse');
end

function assertPositiveDefinite(matrix, errorId, message)
[~, flag] = chol(0.5 * (matrix + matrix'), 'lower');
assert(flag == 0, errorId, message);
end

function assertNumericallyPsd(covariance)
symmetryScale = max(1, norm(covariance, 'fro'));
assert(norm(covariance - covariance', 'fro') <= ...
    1e-12 * symmetryScale, 'dynsc:AsymmetricPosteriorCovariance', ...
    'PPost must be numerically symmetric.');
minimumEigenvalue = min(eig(covariance));
psdTolerance = 1e-10 * max(1, norm(covariance, 2));
assert(minimumEigenvalue >= -psdTolerance, ...
    'dynsc:IndefinitePosteriorCovariance', ...
    'PPost must be numerically positive semidefinite.');
end

function terms = decomposeLinearizedObjective( ...
    priorResidual, priorPrecision, measurementResidual, ...
    noiseCovariance, linearizedIncResidual, linearizedBoxResidual, ...
    sigmaSquaredInc, sigmaSquaredBox)
offDiagonalNoise = noiseCovariance - diag(diag(noiseCovariance));
noiseScale = max(1, norm(noiseCovariance, 'fro'));
assert(norm(offDiagonalNoise, 'fro') <= 1e-12 * noiseScale, ...
    'dynsc:NonDiagonalMeasurementCovariance', ...
    ['The objective decomposition requires independent linear virtual-' ...
     'measurement errors.']);

measurementVariances = diag(noiseCovariance);
terms = struct();
terms.prior = 0.5 * ...
    (priorResidual' * priorPrecision * priorResidual);
terms.smoothness = 0.5 * ...
    measurementResidual(1)^2 / measurementVariances(1);
terms.edgeTotal = 0.5 * ...
    measurementResidual(2)^2 / measurementVariances(2);
terms.triangleTotal = 0.5 * ...
    measurementResidual(3)^2 / measurementVariances(3);
terms.measurementTotal = terms.smoothness + ...
    terms.edgeTotal + terms.triangleTotal;
terms.normalTotal = terms.prior + terms.measurementTotal;
terms.inclusionLinearized = 0.5 * ...
    linearizedIncResidual^2 / sigmaSquaredInc;
terms.boxLinearized = 0.5 * ...
    linearizedBoxResidual^2 / sigmaSquaredBox;
% Retain the common field names used by the existing diagnostics pipeline.
terms.inclusionBarrier = terms.inclusionLinearized;
terms.boxBarrier = terms.boxLinearized;
terms.barrierTotal = terms.inclusionLinearized + terms.boxLinearized;
terms.total = terms.normalTotal + terms.barrierTotal;
end
