function [selected, trials] = resolvePrimarySchoolFixedHyperparameters( ...
    fixed, calibratedSmoothnessVariance, supportThreshold)
%RESOLVEPRIMARYSCHOOLFIXEDHYPERPARAMETERS Record one fixed tracking profile.

required = {'smoothnessVariance', 'lambdaInc', 'lambdaBox'};
for iField = 1:numel(required)
    name = required{iField};
    assert(isfield(fixed, name), ...
        'dynsc:MissingFixedPrimarySchoolHyperparameters', ...
        'Missing tuning.fixed.%s.', name);
    validateattributes(fixed.(name), {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
end
validateattributes(calibratedSmoothnessVariance, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(supportThreshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>=', 0, '<=', 1});

smoothnessVariance = fixed.smoothnessVariance;
lambdaInc = fixed.lambdaInc;
lambdaBox = fixed.lambdaBox;
multiplier = smoothnessVariance / calibratedSmoothnessVariance;

selected = struct();
selected.trial = 1;
selected.smoothnessVariance = smoothnessVariance;
selected.smoothnessVarianceMultiplier = multiplier;
selected.lambdaInc = lambdaInc;
selected.lambdaBox = lambdaBox;
selected.selectionScore = NaN;
selected.selectionEstimator = "fixed";
selected.supportThreshold = supportThreshold;
selected.mode = "fixed";

trials = table(1, multiplier, smoothnessVariance, lambdaInc, lambdaBox, ...
    NaN, NaN, NaN, 'VariableNames', ...
    {'Trial', 'SmoothnessVarianceMultiplier', 'SmoothnessVariance', ...
     'LambdaInc', 'LambdaBox', 'PooledEdgeF1', 'PooledTriangleF1', ...
     'SelectionScore'});
end
