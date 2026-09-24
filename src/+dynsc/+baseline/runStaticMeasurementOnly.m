function baselineResult = runStaticMeasurementOnly( ...
    measurement, topology, gurobiParameters)
%RUNSTATICMEASUREMENTONLY Solve each time index independently.

if nargin < 3
    gurobiParameters = [];
end

nState = size(measurement.H, 2);
nTime = size(measurement.H, 3);

assert(size(measurement.y, 2) == nTime, ...
    'dynsc:StaticBaselineTimeDimensionMismatch', ...
    'measurement.H and measurement.y must have matching time dimensions.');

estimatedState = zeros(nState, nTime);
measurementResidual = zeros(size(measurement.y));
objectiveValue = zeros(1, nTime);
gurobiStatus = strings(1, nTime);
smoothnessObjective = zeros(1, nTime);
edgeTotalObjective = zeros(1, nTime);
triangleTotalObjective = zeros(1, nTime);

for iTime = 1:nTime
    estimate = dynsc.baseline.solveMeasurementOnlyGurobi( ...
        measurement.y(:, iTime), ...
        measurement.H(:, :, iTime), ...
        measurement.noiseCovariance, topology, gurobiParameters);

    estimatedState(:, iTime) = estimate.state;
    measurementResidual(:, iTime) = estimate.measurementResidual;
    objectiveValue(iTime) = estimate.objectiveValue;
    gurobiStatus(iTime) = string(estimate.gurobiResult.status);
    smoothnessObjective(iTime) = estimate.objectiveTerms.smoothness;
    edgeTotalObjective(iTime) = estimate.objectiveTerms.edgeTotal;
    triangleTotalObjective(iTime) = ...
        estimate.objectiveTerms.triangleTotal;
end

baselineResult = struct();
baselineResult.name = "Static";
baselineResult.estimatedState = estimatedState;
baselineResult.measurementResidual = measurementResidual;
baselineResult.objectiveValue = objectiveValue;
baselineResult.gurobiStatus = gurobiStatus;
baselineResult.objectiveTerms = struct();
baselineResult.objectiveTerms.smoothness = smoothnessObjective;
baselineResult.objectiveTerms.edgeTotal = edgeTotalObjective;
baselineResult.objectiveTerms.triangleTotal = triangleTotalObjective;
baselineResult.objectiveTerms.total = smoothnessObjective + ...
    edgeTotalObjective + triangleTotalObjective;
end
