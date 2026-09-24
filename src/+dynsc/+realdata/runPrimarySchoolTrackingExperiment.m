function result = runPrimarySchoolTrackingExperiment(cfg)
%RUNPRIMARYSCHOOLTRACKINGEXPERIMENT Calibrate and test empirical SC tracking.
%
% This adapter intentionally treats all existing prediction, measurement,
% correction, metric, and plotting implementations as immutable functions.

validateExperimentConfig(cfg);
if ~isfolder(cfg.output.root)
    mkdir(cfg.output.root);
end
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
runDirectory = fullfile(cfg.output.root, ...
    sprintf('%s_%s', matlab.lang.makeValidName(char(cfg.experiment.name)), ...
    timestamp));
mkdir(runDirectory);
save(fullfile(runDirectory, 'config.mat'), 'cfg');

previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng));
nRealizations = cfg.experiment.numRealizations;
realizations = cell(nRealizations, 1);

fprintf('Running primary-school calibration/test experiment.\n');
fprintf('Calibration day %d, test day %d, segment %s for %g hours.\n', ...
    cfg.segment.calibrationDay, cfg.segment.testDay, ...
    cfg.segment.startClock, cfg.segment.durationHours);

for iRealization = 1:nRealizations
    realizationDirectory = fullfile(runDirectory, ...
        sprintf('realization_%03d', iRealization));
    mkdir(realizationDirectory);

    realizationCfg = cfg;
    realizationCfg.selection.seed = cfg.selection.seed + iRealization - 1;
    realizationSeed = cfg.experiment.seed + iRealization - 1;
    rng(realizationSeed, 'twister');

    dataset = buildConfiguredDataset(realizationCfg);
    calibrationTopology = dynsc.realdata.extractTopologySegment( ...
        dataset.days(cfg.segment.calibrationDay).topology, ...
        cfg.segment.startClock, cfg.segment.durationHours);
    testTopology = dynsc.realdata.extractTopologySegment( ...
        dataset.days(cfg.segment.testDay).topology, ...
        cfg.segment.startClock, cfg.segment.durationHours);

    signalConfig = cfg.signals;
    signalConfig.numTimeSteps = calibrationTopology.numTimeSteps;
    calibrationSignals = signalConfig.generator( ...
        calibrationTopology, signalConfig);
    signalConfig.numTimeSteps = testTopology.numTimeSteps;
    testSignals = signalConfig.generator(testTopology, signalConfig);

    calibration = dynsc.realdata.calibratePrimarySchoolModel( ...
        calibrationTopology, calibrationSignals, cfg.measurement, ...
        cfg.calibration);
    predictionConfig = resolvedPredictionConfig( ...
        cfg.prediction, calibration, calibrationTopology);

    [selected, tuningTrials] = resolveHyperparameters( ...
        calibrationTopology, calibrationSignals, predictionConfig, ...
        calibration, cfg);
    updateConfig = cfg.update;
    updateConfig.barrier.lambda_inc = selected.lambdaInc;
    updateConfig.barrier.lambda_box = selected.lambdaBox;
    measurementConfig = dynsc.realdata.configurePrimarySchoolMeasurement( ...
        cfg.measurement, calibration, selected.smoothnessVariance);

    calibrationMeasurement = measurementConfig.builder( ...
        calibrationTopology, calibrationSignals, measurementConfig);
    testMeasurement = measurementConfig.builder( ...
        testTopology, testSignals, measurementConfig);

    calibrationMethods = runMethodSuite( ...
        predictionConfig, calibrationMeasurement, ...
        calibrationTopology, updateConfig, cfg);
    testMethods = runMethodSuite( ...
        predictionConfig, testMeasurement, testTopology, ...
        updateConfig, cfg);

    saveMetricFigure(calibrationMethods, false, ...
        fullfile(realizationDirectory, 'calibration_recovery_metrics.pdf'));
    saveMetricFigure(testMethods, false, ...
        fullfile(realizationDirectory, 'test_recovery_metrics.pdf'));

    diagnosticFiles = struct();
    diagnosticFiles.calibration = saveMethodDiagnostics( ...
        calibrationMethods, "calibration", realizationDirectory);
    diagnosticFiles.test = saveMethodDiagnostics( ...
        testMethods, "test", realizationDirectory);
    testSnapshots = saveConfiguredTestSnapshots( ...
        testMethods, testTopology, cfg, realizationDirectory);

    tuningFile = fullfile(realizationDirectory, ...
        'calibration_hyperparameter_trials.csv');
    writetable(tuningTrials, tuningFile);
    selectedNodesFile = fullfile(realizationDirectory, ...
        'selected_nodes.csv');
    writetable(dataset.selection.nodes, selectedNodesFile);

    realization = struct();
    realization.index = iRealization;
    realization.seed = realizationSeed;
    realization.selection = dataset.selection;
    realization.calibrationTopology = calibrationTopology;
    realization.testTopology = testTopology;
    realization.calibrationSignals = calibrationSignals;
    realization.testSignals = testSignals;
    realization.calibration = calibration;
    realization.countSource = measurementConfig.countSource;
    realization.selectedHyperparameters = selected;
    realization.tuningTrials = tuningTrials;
    realization.calibrationMeasurement = calibrationMeasurement;
    realization.testMeasurement = testMeasurement;
    realization.calibrationMethods = calibrationMethods;
    realization.testMethods = testMethods;
    realization.testSnapshots = testSnapshots;
    realization.files.tuning = tuningFile;
    realization.files.selectedNodes = selectedNodesFile;
    realization.files.diagnostics = diagnosticFiles;

    realizationFile = fullfile(realizationDirectory, 'result.mat');
    save(realizationFile, 'realization', '-v7.3');
    realizations{iRealization} = realization;

    fprintf(['  Realization %d: beta=[%.4f %.4f], rho=[%.4f %.4f], ' ...
        'C=[%.2f %.2f], selected trial %d.\n'], iRealization, ...
        calibration.transition.beta1, calibration.transition.beta2, ...
        calibration.transition.rho1, calibration.transition.rho2, ...
        calibration.C1, calibration.C2, selected.trial);
end

result = struct();
result.experimentName = cfg.experiment.name;
result.outputDirectory = runDirectory;
result.realizations = realizations;
save(fullfile(runDirectory, 'result.mat'), 'result', 'cfg', '-v7.3');
fprintf('Finished. Results saved in:\n%s\n', runDirectory);
end

function dataset = buildConfiguredDataset(cfg)
selectionMode = "two_class";
if isfield(cfg.selection, 'mode')
    selectionMode = string(cfg.selection.mode);
end
switch selectionMode
    case "single_class"
        dataset = dynsc.realdata.buildPrimarySchoolSingleClassDataset(cfg);
    case "two_class"
        dataset = dynsc.realdata.buildPrimarySchoolDataset(cfg);
    otherwise
        error('dynsc:InvalidCohortMode', ...
            'selection.mode must be "single_class" or "two_class".');
end
end

function predictionConfig = resolvedPredictionConfig(base, calibration, topology)
predictionConfig = base;
nState = size(topology.s1, 1) + size(topology.s2, 1);
predictionConfig.model = "closure_aware";
predictionConfig.transition = calibration.transition;
predictionConfig.F = speye(nState);
predictionConfig.Q = "bernoulli_second_moment";
predictionConfig.initialization.method = "manual";
predictionConfig.initialState = calibration.initialState;
predictionConfig.initialCovariance = full(calibration.initialCovariance);
predictionConfig = dynsc.prediction.validateParameters( ...
    predictionConfig, nState);
end

function [selected, trials] = resolveHyperparameters( ...
    topology, signals, predictionConfig, calibration, cfg)
mode = trackingHyperparameterMode(cfg.tuning);
switch mode
    case "fixed"
        [selected, trials] = ...
            dynsc.realdata.resolvePrimarySchoolFixedHyperparameters( ...
                cfg.tuning.fixed, calibration.smoothnessVariance, ...
                cfg.metrics.supportThreshold);
        fprintf(['  Using fixed hyperparameters: sigma_s^2=%.4g, ' ...
            'lambda=[%.4g %.4g].\n'], ...
            selected.smoothnessVariance, selected.lambdaInc, ...
            selected.lambdaBox);
    case "calibrated_screen"
        [selected, trials] = tuneOnCalibration( ...
            topology, signals, predictionConfig, calibration, cfg);
    otherwise
        error('dynsc:InvalidPrimarySchoolTuningMode', ...
            ['tuning.mode must be "fixed" or ' ...
             '"calibrated_screen".']);
end
end

function [selected, trials] = tuneOnCalibration( ...
    topology, signals, predictionConfig, calibration, cfg)
multipliers = cfg.tuning.smoothnessVarianceMultiplier(:);
lambdaInc = cfg.tuning.lambdaInc(:);
lambdaBox = cfg.tuning.lambdaBox(:);
assert(numel(multipliers) == numel(lambdaInc) && ...
    numel(multipliers) == numel(lambdaBox), ...
    'dynsc:InconsistentTuningCandidates', ...
    'All hyperparameter candidate vectors must have equal length.');

if ~cfg.tuning.enabled
    multipliers = multipliers(1);
    lambdaInc = lambdaInc(1);
    lambdaBox = lambdaBox(1);
end
nTrials = numel(multipliers);
score = zeros(nTrials, 1);
edgeF1 = zeros(nTrials, 1);
triangleF1 = zeros(nTrials, 1);
smoothnessVariance = calibration.smoothnessVariance * multipliers;

fprintf('  Screening %d day-1 hyperparameter profile(s).\n', nTrials);
for iTrial = 1:nTrials
    measurementConfig = dynsc.realdata.configurePrimarySchoolMeasurement( ...
        cfg.measurement, calibration, smoothnessVariance(iTrial));
    measurement = measurementConfig.builder( ...
        topology, signals, measurementConfig);
    updateConfig = cfg.update;
    updateConfig.method = updateMethodForEstimator(cfg.tuning.estimator);
    updateConfig.barrier.lambda_inc = lambdaInc(iTrial);
    updateConfig.barrier.lambda_box = lambdaBox(iTrial);

    filterResult = dynsc.filter.runPredictionCorrectionLoop( ...
        predictionConfig, measurement, topology, [], updateConfig);
    [edgeF1(iTrial), triangleF1(iTrial)] = pooledSupportF1( ...
        filterResult.updatedState, topology, cfg.metrics.supportThreshold);
    score(iTrial) = 0.5 * (edgeF1(iTrial) + triangleF1(iTrial));
    fprintf(['    Trial %d/%d: sigma_s^2=%.4g, lambda=[%.4g %.4g], ' ...
        'pooled F1=[%.3f %.3f], score=%.3f.\n'], ...
        iTrial, nTrials, smoothnessVariance(iTrial), ...
        lambdaInc(iTrial), lambdaBox(iTrial), ...
        edgeF1(iTrial), triangleF1(iTrial), score(iTrial));
end

trials = table((1:nTrials).', multipliers, smoothnessVariance, ...
    lambdaInc, lambdaBox, edgeF1, triangleF1, score, ...
    'VariableNames', {'Trial', 'SmoothnessVarianceMultiplier', ...
    'SmoothnessVariance', 'LambdaInc', 'LambdaBox', ...
    'PooledEdgeF1', 'PooledTriangleF1', 'SelectionScore'});
[~, best] = max(score);
selected = struct();
selected.trial = best;
selected.smoothnessVariance = smoothnessVariance(best);
selected.smoothnessVarianceMultiplier = multipliers(best);
selected.lambdaInc = lambdaInc(best);
selected.lambdaBox = lambdaBox(best);
selected.selectionScore = score(best);
selected.selectionEstimator = string(cfg.tuning.estimator);
selected.supportThreshold = cfg.metrics.supportThreshold;
selected.mode = "calibrated_screen";
end

function mode = trackingHyperparameterMode(tuning)
% Older configurations without a mode retain the original behavior.
if isfield(tuning, 'mode') && ~isempty(tuning.mode)
    mode = string(tuning.mode);
else
    mode = "calibrated_screen";
end
assert(isscalar(mode), 'dynsc:InvalidPrimarySchoolTuningMode', ...
    'tuning.mode must be scalar text.');
end

function method = updateMethodForEstimator(estimator)
switch string(estimator)
    case "ekf"
        method = "barrier_standard_ekf_projection";
    case "iekf"
        method = "barrier_virtual_measurement";
    case "laplace"
        method = "barrier_exact_hessian";
    otherwise
        error('dynsc:InvalidTuningEstimator', ...
            'Tuning estimator must be "ekf", "iekf", or "laplace".');
end
end

function methods = runMethodSuite(predictionConfig, measurement, topology, updateConfig, cfg)
methods = struct();
if cfg.methods.constrainedRls
    rlsPrediction = predictionConfig;
    nState = numel(rlsPrediction.initialState);
    rlsPrediction.model = "linear";
    rlsPrediction.F = speye(nState);
    rlsPrediction.Q = sparse(nState, nState);
    methods.constrainedRls = runOne( ...
        "Constrained RLS", rlsPrediction, "constrained_map", ...
        measurement, topology, updateConfig, cfg.metrics, ...
        cfg.output.saveDetailedFilters);
else
    methods.constrainedRls = [];
end
if cfg.methods.constrainedMap
    methods.constrainedMap = runOne( ...
        "Constrained MAP", predictionConfig, "constrained_map", ...
        measurement, topology, updateConfig, cfg.metrics, ...
        cfg.output.saveDetailedFilters);
else
    methods.constrainedMap = [];
end
if cfg.methods.ekf
    methods.ekf = runOne( ...
        "EKF", predictionConfig, "barrier_standard_ekf_projection", ...
        measurement, topology, updateConfig, cfg.metrics, ...
        cfg.output.saveDetailedFilters);
else
    methods.ekf = [];
end
if cfg.methods.iekf
    methods.iekf = runOne( ...
        "IEKF", predictionConfig, "barrier_virtual_measurement", ...
        measurement, topology, updateConfig, cfg.metrics, ...
        cfg.output.saveDetailedFilters);
else
    methods.iekf = [];
end
if cfg.methods.laplace
    methods.laplace = runOne( ...
        "Laplace", predictionConfig, "barrier_exact_hessian", ...
        measurement, topology, updateConfig, cfg.metrics, ...
        cfg.output.saveDetailedFilters);
else
    methods.laplace = [];
end
end

function method = runOne(name, prediction, updateMethod, measurement, topology, updateConfig, metricsConfig, keepDetails)
updateConfig.method = updateMethod;
filterResult = dynsc.filter.runPredictionCorrectionLoop( ...
    prediction, measurement, topology, [], updateConfig);
filterResult.name = name;
metrics = dynsc.metrics.evaluateRecoveryOverTime( ...
    filterResult, topology, metricsConfig);
method = struct();
method.name = name;
method.metrics = metrics;
method.diagnostics = ...
    dynsc.realdata.buildPrimarySchoolStateDiagnostics( ...
        filterResult, topology, metrics);
if keepDetails
    method.filter = filterResult;
else
    method.filter = compactFilter(filterResult);
end
end

function compact = compactFilter(filterResult)
compact = struct();
compact.name = filterResult.name;
compact.predictionModel = filterResult.predictionModel;
compact.updateMethod = filterResult.updateMethod;
compact.initialState = filterResult.initialState;
compact.initialCovariance = filterResult.initialCovariance;
compact.predictedState = filterResult.predictedState;
compact.updatedState = filterResult.updatedState;
compact.priorResidual = filterResult.priorResidual;
compact.measurementResidual = filterResult.measurementResidual;
compact.objectiveValue = filterResult.objectiveValue;
compact.solverStatus = filterResult.solverStatus;
end

function [edgeF1, triangleF1] = pooledSupportF1(state, topology, threshold)
nEdges = size(topology.s1, 1);
edgeF1 = binaryF1(logical(topology.s1(:)), ...
    state(1:nEdges, :) >= threshold);
triangleF1 = binaryF1(logical(topology.s2(:)), ...
    state(nEdges + 1:end, :) >= threshold);
end

function value = binaryF1(truth, estimate)
truth = truth(:);
estimate = estimate(:);
tp = nnz(truth & estimate);
fp = nnz(~truth & estimate);
fn = nnz(truth & ~estimate);
denominator = 2 * tp + fp + fn;
if denominator == 0
    value = 1;
else
    value = 2 * tp / denominator;
end
end

function saveMetricFigure(methods, visible, path)
assert(~isempty(methods.iekf), 'dynsc:MissingIekfForPlot', ...
    'The paper plotting pipeline requires IEKF metrics as its reference.');
aggregate = singleRealizationAggregate(methods);
fig = dynsc.plot.plotAverageRecoveryMetrics(aggregate, visible, false);
dynsc.plot.exportFigurePdf(fig, path);
close(fig);
end

function aggregate = singleRealizationAggregate(methods)
realization = struct();
realization.metrics = methods.iekf.metrics;
if ~isempty(methods.constrainedRls)
    realization.constrainedRlsBaseline = metricWrapper( ...
        methods.constrainedRls.metrics);
end
if ~isempty(methods.constrainedMap)
    realization.previousUpdateBaseline = metricWrapper( ...
        methods.constrainedMap.metrics);
end
if ~isempty(methods.ekf)
    realization.standardEkfBarrierBaseline = metricWrapper( ...
        methods.ekf.metrics);
end
if ~isempty(methods.laplace)
    realization.exactHessianBaseline = metricWrapper( ...
        methods.laplace.metrics);
end
aggregate = dynsc.metrics.aggregateRecoveryMetrics({realization});
end

function wrapped = metricWrapper(metrics)
wrapped = struct('metrics', metrics);
end

function files = saveMethodDiagnostics(methods, prefix, directory)
files = struct();
methodFields = fieldnames(methods);
for iMethod = 1:numel(methodFields)
    fieldName = methodFields{iMethod};
    method = methods.(fieldName);
    if isempty(method)
        continue;
    end
    path = fullfile(directory, sprintf('%s_%s_diagnostics.csv', ...
        char(prefix), fieldName));
    writetable(method.diagnostics, path);
    files.(fieldName) = string(path);
end
end

function snapshot = saveConfiguredTestSnapshots( ...
    methods, topology, cfg, directory)
snapshot = struct('enabled', logical(cfg.testSnapshots.enabled));
if ~snapshot.enabled
    return;
end

methodField = char(string(cfg.testSnapshots.method));
assert(isfield(methods, methodField) && ...
    ~isempty(methods.(methodField)), ...
    'dynsc:UnavailablePrimarySchoolSnapshotMethod', ...
    'The requested test snapshot method "%s" is not enabled.', ...
    methodField);
method = methods.(methodField);
[figures, indices] = dynsc.realdata.plotSingleStateSnapshots( ...
    method.filter, topology, method.name, ...
    string(cfg.plot.visible) == "on", cfg.testSnapshots.interval);

pdfFiles = strings(1, numel(indices));
figFiles = strings(1, numel(indices));
for iSnapshot = 1:numel(indices)
    baseName = sprintf('test_%s_state_snapshot_k%04d', ...
        methodField, indices(iSnapshot));
    pdfFiles(iSnapshot) = fullfile(directory, baseName + ".pdf");
    figFiles(iSnapshot) = fullfile(directory, baseName + ".fig");
    exportgraphics(figures(iSnapshot), pdfFiles(iSnapshot), ...
        'ContentType', 'vector');
    savefig(figures(iSnapshot), figFiles(iSnapshot));
    close(figures(iSnapshot));
end

snapshot.method = string(methodField);
snapshot.methodName = method.name;
snapshot.indices = indices;
snapshot.pdf = pdfFiles;
snapshot.fig = figFiles;
end

function validateExperimentConfig(cfg)
validateattributes(cfg.experiment.numRealizations, {'numeric'}, ...
    {'scalar', 'integer', 'positive'});
assert(cfg.segment.calibrationDay ~= cfg.segment.testDay, ...
    'dynsc:CalibrationTestDayOverlap', ...
    'Calibration and test days must be distinct.');
% The aggregation width is configurable, but the construction rule is not:
% edges are pair contacts anywhere in the window and a triangle still
% requires its three nodes in one raw maximal simplex, so a wider window
% never fills an aggregate 3-clique. A width below 10 minutes has not been
% vetted, and the width must tile the segment exactly.
validateattributes(cfg.window.durationMinutes, {'numeric'}, ...
    {'scalar', 'positive', 'finite', 'integer', '>=', 10}, ...
    '', 'cfg.window.durationMinutes');
assert(mod(cfg.window.durationMinutes, cfg.window.rawResolutionSeconds / 60) == 0, ...
    'dynsc:UnexpectedTopologyWindow', ...
    'The window width must be a whole number of raw sampling intervals.');
assert(mod(cfg.segment.durationHours * 60, cfg.window.durationMinutes) == 0, ...
    'dynsc:UnexpectedTopologyWindow', ...
    ['The segment duration must be an exact multiple of the topology ' ...
     'window width.']);
assert(cfg.signals.smoothnessType == cfg.measurement.smoothnessType, ...
    'dynsc:SignalMeasurementPriorMismatch', ...
    'Signal generation and measurement must use the same smoothness type.');
mode = trackingHyperparameterMode(cfg.tuning);
assert(any(mode == ["fixed", "calibrated_screen"]), ...
    'dynsc:InvalidPrimarySchoolTuningMode', ...
    'tuning.mode must be "fixed" or "calibrated_screen".');
if mode == "fixed"
    assert(isfield(cfg.tuning, 'fixed') && isstruct(cfg.tuning.fixed), ...
        'dynsc:MissingFixedPrimarySchoolHyperparameters', ...
        'tuning.fixed must contain the fixed hyperparameters.');
    required = {'smoothnessVariance', 'lambdaInc', 'lambdaBox'};
    for iField = 1:numel(required)
        name = required{iField};
        assert(isfield(cfg.tuning.fixed, name), ...
            'dynsc:MissingFixedPrimarySchoolHyperparameters', ...
            'Missing tuning.fixed.%s.', name);
        validateattributes(cfg.tuning.fixed.(name), {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'});
    end
end
assert(isfield(cfg, 'testSnapshots') && ...
    isfield(cfg.testSnapshots, 'enabled'), ...
    'dynsc:MissingPrimarySchoolSnapshotConfig', ...
    'Supply cfg.testSnapshots.enabled for the tracking experiment.');
validateattributes(cfg.testSnapshots.enabled, {'logical', 'numeric'}, ...
    {'scalar'});
if logical(cfg.testSnapshots.enabled)
    assert(isfield(cfg.testSnapshots, 'method') && ...
        isfield(cfg.testSnapshots, 'interval'), ...
        'dynsc:MissingPrimarySchoolSnapshotConfig', ...
        ['Enabled test snapshots require testSnapshots.method and ' ...
         'testSnapshots.interval.']);
    method = string(cfg.testSnapshots.method);
    assert(isscalar(method) && any(method == ...
        ["constrainedRls", "constrainedMap", "ekf", "iekf", ...
         "laplace"]), ...
        'dynsc:InvalidPrimarySchoolSnapshotMethod', ...
        'Invalid testSnapshots.method.');
    validateattributes(cfg.testSnapshots.interval, {'numeric'}, ...
        {'scalar', 'integer', 'positive', 'finite'});
end
end
