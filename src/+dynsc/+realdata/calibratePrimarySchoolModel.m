function calibration = calibratePrimarySchoolModel( ...
    topology, signals, measurementConfig, calibrationConfig)
%CALIBRATEPRIMARYSCHOOLMODEL Fit transition and measurement quantities.
%
% Every quantity returned here uses only the supplied calibration topology
% and signals. No test topology is accepted by this function.

alpha = calibrationConfig.betaPriorAlpha;
beta = calibrationConfig.betaPriorBeta;
varianceFloor = calibrationConfig.varianceFloor;
% alpha = beta = 0 is permitted and reduces the estimator below to the
% plain maximum-likelihood ratio successes/trials. It is only safe while
% every transition category supplies at least one trial, which is asserted
% once the counts are known.
validateattributes([alpha, beta], {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'finite', 'nonnegative'});
validateattributes(varianceFloor, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});

X = logical(topology.s1);
Y = logical(topology.s2);
G = logical(topology.boundaryGate);
assert(size(X, 2) >= 2, 'dynsc:InsufficientCalibrationTransitions', ...
    'At least two topology states are required for calibration.');

edgePrevious = X(:, 1:end-1);
edgeCurrent = X(:, 2:end);
trianglePrevious = Y(:, 1:end-1);
triangleCurrent = Y(:, 2:end);
currentGate = G(:, 2:end);

counts = struct();
counts.edgeBirthSuccess = nnz(~edgePrevious & edgeCurrent);
counts.edgeBirthTrials = nnz(~edgePrevious);
counts.edgeSurvivalSuccess = nnz(edgePrevious & edgeCurrent);
counts.edgeSurvivalTrials = nnz(edgePrevious);
counts.triangleBirthSuccess = nnz( ...
    ~trianglePrevious & currentGate & triangleCurrent);
counts.triangleBirthTrials = nnz(~trianglePrevious & currentGate);
counts.triangleSurvivalSuccess = nnz( ...
    trianglePrevious & currentGate & triangleCurrent);
counts.triangleSurvivalTrials = nnz(trianglePrevious & currentGate);

if alpha + beta == 0
    trialCounts = [counts.edgeBirthTrials, counts.edgeSurvivalTrials, ...
        counts.triangleBirthTrials, counts.triangleSurvivalTrials];
    assert(all(trialCounts > 0), ...
        'dynsc:EmptyTransitionCategory', ...
        ['An unsmoothed transition estimate requires at least one trial ' ...
         'in every category; observed trials were [%s].'], ...
        num2str(trialCounts));
end

posteriorMean = @(successes, trials) ...
    (successes + alpha) / (trials + alpha + beta);
transition.beta1 = posteriorMean( ...
    counts.edgeBirthSuccess, counts.edgeBirthTrials);
transition.rho1 = posteriorMean( ...
    counts.edgeSurvivalSuccess, counts.edgeSurvivalTrials);
transition.beta2 = posteriorMean( ...
    counts.triangleBirthSuccess, counts.triangleBirthTrials);
transition.rho2 = posteriorMean( ...
    counts.triangleSurvivalSuccess, counts.triangleSurvivalTrials);

C1 = double(topology.C1_true(:));
C2 = double(topology.C2_true(:));
countMean = [mean(C1); mean(C2)];
countVariance = [sampleVariance(C1); sampleVariance(C2)];
countVariance = max(countVariance, varianceFloor);

temporaryMeasurementConfig = measurementConfig;
temporaryMeasurementConfig.C1 = countMean(1);
temporaryMeasurementConfig.C2 = countMean(2);
temporaryMeasurementConfig.noiseVariances = ones(3, 1);
measurement = dynsc.measurement.buildModel( ...
    topology, signals, temporaryMeasurementConfig);
trueState = double([topology.s1; topology.s2]);
smoothnessResidual = sum(measurement.h .* trueState, 1);
% The virtual observation is fixed at zero, so the zero-mean Gaussian MLE
% is its second moment about zero, rather than variance about the sample mean.
smoothnessVariance = mean(smoothnessResidual .^ 2);
smoothnessVariance = max(smoothnessVariance, varianceFloor);

nEdges = size(topology.s1, 1);
nTriangles = size(topology.s2, 1);
edgeProbability = countMean(1) / nEdges;
triangleProbability = countMean(2) / nTriangles;
assert(triangleProbability <= edgeProbability + 1e-12, ...
    'dynsc:InfeasibleEmpiricalInitialization', ...
    ['The uniform triangle initialization exceeds the uniform edge ' ...
     'initialization and would violate relaxed inclusion.']);

initialState = [
    edgeProbability * ones(nEdges, 1)
    triangleProbability * ones(nTriangles, 1)
];
initialVariance = [
    edgeProbability * (1 - edgeProbability) * ones(nEdges, 1)
    triangleProbability * (1 - triangleProbability) * ...
        ones(nTriangles, 1)
];
initialVariance = max(initialVariance, varianceFloor);

calibration = struct();
calibration.transition = transition;
calibration.transitionCounts = counts;
calibration.transitionPrior = [alpha, beta];
calibration.C1 = countMean(1);
calibration.C2 = countMean(2);
% The full calibration-day count sequences, retained so that a test segment
% covering the same clock interval can use the time-matched counts instead
% of their mean. These are observed on the calibration day only.
calibration.C1Sequence = reshape(C1, 1, []);
calibration.C2Sequence = reshape(C2, 1, []);
calibration.countVariances = countVariance;
calibration.smoothnessResidual = smoothnessResidual;
calibration.smoothnessVariance = smoothnessVariance;
calibration.initialEdgeProbability = edgeProbability;
calibration.initialTriangleProbability = triangleProbability;
calibration.initialState = initialState;
calibration.initialCovariance = spdiags( ...
    initialVariance, 0, nEdges + nTriangles, nEdges + nTriangles);
calibration.numStates = topology.numTimeSteps;
calibration.sourceDay = topology.dayIndex;
calibration.segmentStart = topology.segmentStart;
calibration.segmentEnd = topology.segmentEnd;
end

function value = sampleVariance(x)
if numel(x) > 1
    value = var(x, 0);
else
    value = 0;
end
end
