function filterResult = runPredictionCorrectionLoop( ...
    predictionConfig, measurement, topology, gurobiParameters, updateConfig)
%RUNPREDICTIONCORRECTIONLOOP Alternate prediction and constrained correction.

if nargin < 4
    gurobiParameters = [];
end
if nargin < 5 || isempty(updateConfig)
    updateConfig.method = "constrained_map";
end

predictionConfig = dynsc.prediction.validateParameters(predictionConfig);
nState = numel(predictionConfig.initialState);
nTime = size(measurement.H, 3);

assert(size(measurement.y, 2) == nTime, ...
    'dynsc:FilterTimeDimensionMismatch', ...
    'measurement.H and measurement.y must contain the same time indices.');
assert(size(measurement.H, 2) == nState, ...
    'dynsc:FilterStateDimensionMismatch', ...
    'The measurement matrix must have one column per state variable.');

predictedState = zeros(nState, nTime);
predictedCovariance = zeros(nState, nState, nTime);
updatedState = zeros(nState, nTime);
updatedCovariance = zeros(nState, nState, nTime);
priorResidual = zeros(nState, nTime);
measurementResidual = zeros(size(measurement.y));
objectiveValue = zeros(1, nTime);
gurobiStatus = strings(1, nTime);
solverStatus = strings(1, nTime);
solverName = strings(1, nTime);
priorObjective = zeros(1, nTime);
smoothnessObjective = zeros(1, nTime);
edgeTotalObjective = zeros(1, nTime);
triangleTotalObjective = zeros(1, nTime);
normalObjective = zeros(1, nTime);
inclusionBarrierObjective = zeros(1, nTime);
boxBarrierObjective = zeros(1, nTime);
barrierObjective = zeros(1, nTime);
LambdaNormal = cell(1, nTime);
inclusionInformation = cell(1, nTime);
boxInformation = cell(1, nTime);
barrierMeasurement = cell(1, nTime);
barrierLinearizationState = cell(1, nTime);
unprojectedUpdateState = cell(1, nTime);
transitionJacobian = cell(1, nTime);
processNoiseCovariance = cell(1, nTime);
transitionApplied = false(1, nTime);

previousUpdatedState = predictionConfig.initialState;
previousUpdatedCovariance = predictionConfig.initialCovariance;

for iTime = 1:nTime
    skipPrediction = iTime == 1 && ...
        predictionConfig.skipFirstPrediction;
    if skipPrediction
        % The initial pair is the prior for the first observed ER topology.
        % Assimilate y_1 before applying the Markov model, which starts at k=2.
        predictedState(:, iTime) = previousUpdatedState;
        predictedCovariance(:, :, iTime) = ...
            previousUpdatedCovariance;
    else
        [predictedState(:, iTime), predictedCovariance(:, :, iTime), ...
            transitionJacobian{iTime}, processNoiseCovariance{iTime}] = ...
            dynsc.prediction.predictConfiguredState( ...
                previousUpdatedState, previousUpdatedCovariance, ...
                predictionConfig, topology);
        transitionApplied(iTime) = true;
    end

    update = dynsc.update.solveConfiguredUpdate( ...
        predictedState(:, iTime), ...
        predictedCovariance(:, :, iTime), ...
        measurement.y(:, iTime), ...
        measurement.H(:, :, iTime), ...
        measurement.noiseCovariance, topology, updateConfig, ...
        gurobiParameters);

    updatedState(:, iTime) = update.state;
    updatedCovariance(:, :, iTime) = update.covariance;
    priorResidual(:, iTime) = update.priorResidual;
    measurementResidual(:, iTime) = update.measurementResidual;
    objectiveValue(iTime) = update.objectiveValue;
    solverStatus(iTime) = string(update.solverStatus);
    solverName(iTime) = string(update.solverName);
    if isfield(update, 'gurobiResult') && ~isempty(update.gurobiResult)
        gurobiStatus(iTime) = string(update.gurobiResult.status);
    end
    priorObjective(iTime) = update.objectiveTerms.prior;
    smoothnessObjective(iTime) = update.objectiveTerms.smoothness;
    edgeTotalObjective(iTime) = update.objectiveTerms.edgeTotal;
    triangleTotalObjective(iTime) = ...
        update.objectiveTerms.triangleTotal;
    normalObjective(iTime) = getObjectiveTerm( ...
        update.objectiveTerms, 'normalTotal', ...
        update.objectiveTerms.prior + ...
        update.objectiveTerms.measurementTotal);
    inclusionBarrierObjective(iTime) = getObjectiveTerm( ...
        update.objectiveTerms, 'inclusionBarrier', 0);
    boxBarrierObjective(iTime) = getObjectiveTerm( ...
        update.objectiveTerms, 'boxBarrier', 0);
    barrierObjective(iTime) = getObjectiveTerm( ...
        update.objectiveTerms, 'barrierTotal', 0);
    if isfield(update, 'LambdaNormal')
        LambdaNormal{iTime} = update.LambdaNormal;
    else
        LambdaNormal{iTime} = update.informationMatrix;
    end
    if isfield(update, 'informationTerms')
        inclusionInformation{iTime} = ...
            update.informationTerms.inclusion;
        boxInformation{iTime} = update.informationTerms.box;
    end
    if isfield(update, 'barrierMeasurement')
        barrierMeasurement{iTime} = update.barrierMeasurement;
    end
    if isfield(update, 'linearizationState')
        barrierLinearizationState{iTime} = ...
            update.linearizationState;
    end
    if isfield(update, 'unprojectedState')
        unprojectedUpdateState{iTime} = update.unprojectedState;
    end

    previousUpdatedState = update.state;
    previousUpdatedCovariance = update.covariance;
end

filterResult = struct();
filterResult.initialState = predictionConfig.initialState;
filterResult.initialCovariance = predictionConfig.initialCovariance;
filterResult.predictionModel = ...
    dynsc.prediction.predictionModel(predictionConfig);
filterResult.skipFirstPrediction = ...
    predictionConfig.skipFirstPrediction;
filterResult.transitionApplied = transitionApplied;
filterResult.transitionJacobian = transitionJacobian;
filterResult.processNoiseCovariance = processNoiseCovariance;
filterResult.predictedState = predictedState;
filterResult.predictedCovariance = predictedCovariance;
filterResult.updatedState = updatedState;
filterResult.updatedCovariance = updatedCovariance;
filterResult.priorResidual = priorResidual;
filterResult.measurementResidual = measurementResidual;
filterResult.objectiveValue = objectiveValue;
filterResult.gurobiStatus = gurobiStatus;
filterResult.solverStatus = solverStatus;
filterResult.solverName = solverName;
filterResult.updateMethod = string(updateConfig.method);
filterResult.LambdaNormal = LambdaNormal;
filterResult.inclusionInformation = inclusionInformation;
filterResult.boxInformation = boxInformation;
filterResult.barrierMeasurement = barrierMeasurement;
filterResult.barrierLinearizationState = barrierLinearizationState;
filterResult.unprojectedUpdateState = unprojectedUpdateState;
filterResult.objectiveTerms = struct();
filterResult.objectiveTerms.prior = priorObjective;
filterResult.objectiveTerms.smoothness = smoothnessObjective;
filterResult.objectiveTerms.edgeTotal = edgeTotalObjective;
filterResult.objectiveTerms.triangleTotal = triangleTotalObjective;
filterResult.objectiveTerms.measurementTotal = smoothnessObjective + ...
    edgeTotalObjective + triangleTotalObjective;
filterResult.objectiveTerms.normalTotal = normalObjective;
filterResult.objectiveTerms.inclusionBarrier = ...
    inclusionBarrierObjective;
filterResult.objectiveTerms.boxBarrier = boxBarrierObjective;
filterResult.objectiveTerms.barrierTotal = barrierObjective;
filterResult.objectiveTerms.total = objectiveValue;
end

function value = getObjectiveTerm(terms, name, defaultValue)
if isfield(terms, name)
    value = terms.(name);
else
    value = defaultValue;
end
end
