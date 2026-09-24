function [predictedState, predictedCovariance, F, Q] = ...
    predictConfiguredState( ...
        previousState, previousCovariance, cfg, topology)
%PREDICTCONFIGUREDSTATE Dispatch to the selected prediction model.

model = dynsc.prediction.predictionModel(cfg);
usesTextProcessNoise = isstring(cfg.Q) || ischar(cfg.Q);
needsBoundaryEdges = model == "closure_aware" || usesTextProcessNoise;
if needsBoundaryEdges
    boundaryEdges = dynsc.topology.resolveBoundaryEdges(topology);
else
    boundaryEdges = [];
end
[Q, ~] = dynsc.prediction.processNoiseAtState( ...
    cfg, previousState, boundaryEdges);

switch model
    case "linear"
        F = sparse(cfg.F);
        [predictedState, predictedCovariance] = ...
            dynsc.prediction.predictState( ...
                previousState, previousCovariance, cfg.F, Q);

    case "closure_aware"
        [predictedState, predictedCovariance, F] = ...
            dynsc.prediction.predictClosureAwareState( ...
                previousState, previousCovariance, boundaryEdges, ...
                cfg.transition, Q);
end
end
