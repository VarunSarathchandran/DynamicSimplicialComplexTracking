function prediction = runSequence(cfg, numTimeSteps, topology)
%RUNSEQUENCE Propagate the initialized state through prediction-only steps.
%
% Until a measurement-update step is added, each predicted state and
% covariance become the inputs to the following prediction.

cfg = dynsc.prediction.validateParameters(cfg);
if nargin < 3
    topology = struct();
end
validateattributes(numTimeSteps, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});

nState = numel(cfg.initialState);
predictedState = zeros(nState, numTimeSteps);
predictedCovariance = zeros(nState, nState, numTimeSteps);
transitionJacobian = cell(1, numTimeSteps);
processNoiseCovariance = cell(1, numTimeSteps);

previousState = cfg.initialState;
previousCovariance = cfg.initialCovariance;

for iTime = 1:numTimeSteps
    [predictedState(:, iTime), predictedCovariance(:, :, iTime), ...
        transitionJacobian{iTime}, processNoiseCovariance{iTime}] = ...
        dynsc.prediction.predictConfiguredState( ...
            previousState, previousCovariance, cfg, topology);

    previousState = predictedState(:, iTime);
    previousCovariance = predictedCovariance(:, :, iTime);
end

prediction = struct();
prediction.initialState = cfg.initialState;
prediction.initialCovariance = cfg.initialCovariance;
prediction.predictionModel = dynsc.prediction.predictionModel(cfg);
prediction.transitionJacobian = transitionJacobian;
prediction.processNoiseCovariance = processNoiseCovariance;
prediction.predictedState = predictedState;
prediction.predictedCovariance = predictedCovariance;
end
