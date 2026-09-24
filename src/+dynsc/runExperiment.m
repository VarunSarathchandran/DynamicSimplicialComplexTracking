function result = runExperiment(cfg)
%RUNEXPERIMENT Generate, track, evaluate, and save SC realizations.

cfg = dynsc.validateConfig(cfg);
runDirectory = createRunDirectory(cfg);
save(fullfile(runDirectory, 'config.mat'), 'cfg');

nRealizations = cfg.experiment.numRealizations;
realizations = cell(nRealizations, 1);
previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>

fprintf('Running %s with %d realization(s).\n', ...
    char(cfg.experiment.name), nRealizations);

for iRealization = 1:nRealizations
    realizationSeed = cfg.experiment.seed + iRealization - 1;
    rng(realizationSeed, 'twister');

    topology = cfg.topology.generator( ...
        cfg.topology, cfg.signals.numTimeSteps);
    signals = cfg.signals.generator(topology, cfg.signals);

    realization = struct();
    realization.index = iRealization;
    realization.seed = realizationSeed;
    realization.topology = topology;
    realization.signals = signals;

    if dynsc.measurement.isConfigured(cfg.measurement)
        realization.measurement = cfg.measurement.builder( ...
            topology, signals, cfg.measurement);
    else
        realization.measurement = [];
    end

    predictionConfigured = dynsc.prediction.isConfigured(cfg.prediction);
    measurementConfigured = ...
        dynsc.measurement.isConfigured(cfg.measurement);

    predictionConfig = cfg.prediction;
    realization.predictionInitialization = [];
    realization.identityBaseline = [];
    realization.constrainedRlsBaseline = [];
    realization.previousUpdateBaseline = [];
    realization.exactHessianBaseline = [];
    realization.standardEkfBarrierBaseline = [];
    if predictionConfigured
        [predictionConfig, realization.predictionInitialization] = ...
            dynsc.prediction.resolveInitialization( ...
                predictionConfig, realization.measurement, topology);
    end

    if measurementConfigured && staticBaselineEnabled(cfg)
        realization.staticBaseline = struct();
        realization.staticBaseline.result = cfg.baseline.solver( ...
            realization.measurement, topology);
        realization.staticBaseline.result.name = "Static";
        realization.staticBaseline.metrics = cfg.metrics.evaluator( ...
            realization.staticBaseline.result, topology, cfg.metrics);
        realization.staticBaseline.metrics.estimatorName = "Static";
    else
        realization.staticBaseline = [];
    end

    if predictionConfigured && measurementConfigured
        realization.filter = ...
            dynsc.filter.runPredictionCorrectionLoop( ...
                predictionConfig, realization.measurement, topology, ...
                [], cfg.update);
        realization.filter.name = "IEKF";
        realization.prediction = [];
        realization.metrics = cfg.metrics.evaluator( ...
            realization.filter, topology, cfg.metrics);
        realization.metrics.estimatorName = "IEKF";

        if previousUpdateBaselineEnabled(cfg)
            previousUpdateConfig = cfg.update;
            previousUpdateConfig.method = "constrained_map";
            realization.previousUpdateBaseline = struct();
            realization.previousUpdateBaseline.result = ...
                dynsc.filter.runPredictionCorrectionLoop( ...
                    predictionConfig, realization.measurement, topology, ...
                    [], previousUpdateConfig);
            realization.previousUpdateBaseline.result.name = ...
                "Constrained MAP";
            realization.previousUpdateBaseline.metrics = ...
                cfg.metrics.evaluator( ...
                    realization.previousUpdateBaseline.result, ...
                    topology, cfg.metrics);
            realization.previousUpdateBaseline.metrics.estimatorName = ...
                "Constrained MAP";
        else
            realization.previousUpdateBaseline = [];
        end

        if identityTransitionBaselineEnabled(cfg)
            identityConfig = predictionConfig;
            identityConfig.model = "linear";
            identityConfig.F = speye(numel(identityConfig.initialState));
            % Persistence prior: retain the previous corrected mean and add
            % fixed, simplex-level-specific uncertainty. In particular, do
            % not inherit the primary filter's state-dependent Bernoulli Q_k.
            identityConfig.Q = "block_diagonal";

            % Keep the persistence ablation incremental: it changes the
            % transition prior but retains the historical correction.
            persistenceUpdateConfig = cfg.update;
            persistenceUpdateConfig.method = "constrained_map";

            realization.identityBaseline = struct();
            realization.identityBaseline.result = ...
                dynsc.filter.runPredictionCorrectionLoop( ...
                    identityConfig, realization.measurement, topology, ...
                    [], persistenceUpdateConfig);
            realization.identityBaseline.result.name = ...
                "Persistence";
            realization.identityBaseline.metrics = cfg.metrics.evaluator( ...
                realization.identityBaseline.result, topology, cfg.metrics);
            realization.identityBaseline.metrics.estimatorName = ...
                "Persistence";
        else
            realization.identityBaseline = [];
        end

        if constrainedRlsBaselineEnabled(cfg)
            rlsPredictionConfig = predictionConfig;
            rlsPredictionConfig.model = "linear";
            nRlsState = numel(rlsPredictionConfig.initialState);
            rlsPredictionConfig.F = speye(nRlsState);
            rlsPredictionConfig.Q = sparse(nRlsState, nRlsState);

            rlsUpdateConfig = cfg.update;
            rlsUpdateConfig.method = "constrained_map";

            realization.constrainedRlsBaseline = struct();
            realization.constrainedRlsBaseline.result = ...
                dynsc.filter.runPredictionCorrectionLoop( ...
                    rlsPredictionConfig, realization.measurement, ...
                    topology, [], rlsUpdateConfig);
            realization.constrainedRlsBaseline.result.name = ...
                "Constrained RLS";
            realization.constrainedRlsBaseline.metrics = ...
                cfg.metrics.evaluator( ...
                    realization.constrainedRlsBaseline.result, ...
                    topology, cfg.metrics);
            realization.constrainedRlsBaseline.metrics.estimatorName = ...
                "Constrained RLS";
        else
            realization.constrainedRlsBaseline = [];
        end

        if exactHessianBaselineEnabled(cfg)
            exactHessianUpdateConfig = cfg.update;
            exactHessianUpdateConfig.method = "barrier_exact_hessian";
            realization.exactHessianBaseline = struct();
            realization.exactHessianBaseline.result = ...
                dynsc.filter.runPredictionCorrectionLoop( ...
                    predictionConfig, realization.measurement, topology, ...
                    [], exactHessianUpdateConfig);
            realization.exactHessianBaseline.result.name = ...
                "Laplace";
            realization.exactHessianBaseline.metrics = ...
                cfg.metrics.evaluator( ...
                    realization.exactHessianBaseline.result, ...
                    topology, cfg.metrics);
            realization.exactHessianBaseline.metrics.estimatorName = ...
                "Laplace";
        else
            realization.exactHessianBaseline = [];
        end

        if standardEkfBarrierBaselineEnabled(cfg)
            standardEkfUpdateConfig = cfg.update;
            standardEkfUpdateConfig.method = ...
                "barrier_standard_ekf_projection";
            realization.standardEkfBarrierBaseline = struct();
            realization.standardEkfBarrierBaseline.result = ...
                dynsc.filter.runPredictionCorrectionLoop( ...
                    predictionConfig, realization.measurement, topology, ...
                    [], standardEkfUpdateConfig);
            realization.standardEkfBarrierBaseline.result.name = ...
                "EKF";
            realization.standardEkfBarrierBaseline.metrics = ...
                cfg.metrics.evaluator( ...
                    realization.standardEkfBarrierBaseline.result, ...
                    topology, cfg.metrics);
            realization.standardEkfBarrierBaseline.metrics.estimatorName = ...
                "EKF";
        else
            realization.standardEkfBarrierBaseline = [];
        end

        metricsFigure = dynsc.plot.plotRecoveryMetrics( ...
            realization.metrics, ...
            selectStaticMetrics(realization.staticBaseline), false, ...
            selectIdentityMetrics(realization.identityBaseline), ...
            selectPreviousUpdateMetrics( ...
                realization.previousUpdateBaseline), ...
            selectExactHessianMetrics( ...
                realization.exactHessianBaseline), ...
            selectStandardEkfMetrics( ...
                realization.standardEkfBarrierBaseline), ...
            selectConstrainedRlsMetrics( ...
                realization.constrainedRlsBaseline));
        metricsBaseName = sprintf( ...
            'recovery_metrics_%03d', iRealization);
        metricsPdf = fullfile(runDirectory, ...
            [metricsBaseName, '.pdf']);
        metricsFig = fullfile(runDirectory, ...
            [metricsBaseName, '.fig']);
        exportgraphics(metricsFigure, metricsPdf, ...
            'ContentType', 'vector');
        savefig(metricsFigure, metricsFig);
        close(metricsFigure);

        realization.figureFiles.recoveryMetricsPdf = metricsPdf;
        realization.figureFiles.recoveryMetricsFig = metricsFig;

        [snapshotsEnabled, snapshotInterval] = ...
            stateSnapshotSettings(cfg);
        if snapshotsEnabled
            [snapshotFigures, snapshotTimes] = ...
                dynsc.plot.plotUpdatedStateSnapshots( ...
                    realization.filter, topology, false, ...
                    snapshotInterval, ...
                    selectStaticResult(realization.staticBaseline), ...
                    selectIdentityResult(realization.identityBaseline), ...
                    selectPreviousUpdateResult( ...
                        realization.previousUpdateBaseline), ...
                    selectExactHessianResult( ...
                        realization.exactHessianBaseline), ...
                    selectStandardEkfResult( ...
                        realization.standardEkfBarrierBaseline), ...
                    selectConstrainedRlsResult( ...
                        realization.constrainedRlsBaseline));
            nSnapshots = numel(snapshotTimes);
            snapshotPdf = strings(1, nSnapshots);
            snapshotFig = strings(1, nSnapshots);

            for iSnapshot = 1:nSnapshots
                snapshotBaseName = sprintf( ...
                    'state_snapshot_%03d_k%04d', ...
                    iRealization, snapshotTimes(iSnapshot));
                snapshotPdf(iSnapshot) = fullfile( ...
                    runDirectory, snapshotBaseName + ".pdf");
                snapshotFig(iSnapshot) = fullfile( ...
                    runDirectory, snapshotBaseName + ".fig");
                exportgraphics(snapshotFigures(iSnapshot), ...
                    snapshotPdf(iSnapshot), 'ContentType', 'vector');
                savefig(snapshotFigures(iSnapshot), ...
                    snapshotFig(iSnapshot));
                close(snapshotFigures(iSnapshot));
            end

            realization.figureFiles.stateSnapshotTimes = snapshotTimes;
            realization.figureFiles.stateSnapshotsPdf = snapshotPdf;
            realization.figureFiles.stateSnapshotsFig = snapshotFig;
        end

    elseif predictionConfigured
        realization.prediction = dynsc.prediction.runSequence( ...
            predictionConfig, cfg.signals.numTimeSteps, topology);
        realization.filter = [];
        realization.metrics = [];
    else
        realization.prediction = [];
        realization.filter = [];
        realization.metrics = [];
    end

    realizations{iRealization} = realization;

    if cfg.output.saveRealizations
        outputPath = fullfile(runDirectory, ...
            sprintf('realization_%03d.mat', iRealization));
        save(outputPath, 'realization', '-v7.3');
    end

    fprintf(['  Realization %d/%d: %d nodes, %d edges and ' ...
        '%d triangles at k=1, %d topology changes.\n'], ...
        iRealization, nRealizations, topology.N, ...
        topology.num_edges(1), topology.num_triangles(1), ...
        numel(topology.changeTimes));
end

result = struct();
result.experimentName = cfg.experiment.name;
result.realizations = realizations;
result.outputDirectory = runDirectory;
result.aggregateMetrics = [];
result.figureFiles = struct();

metricsAvailable = cellfun(@(r) ...
    isfield(r, 'metrics') && ~isempty(r.metrics), realizations);
if all(metricsAvailable)
    result.aggregateMetrics = ...
        dynsc.metrics.aggregateRecoveryMetrics(realizations);

    showStandardDeviation = averageStandardDeviationSetting(cfg);
    averageFigure = dynsc.plot.plotAverageRecoveryMetrics( ...
        result.aggregateMetrics, false, showStandardDeviation);
    averagePdf = fullfile(runDirectory, ...
        'average_recovery_metrics.pdf');
    averageFig = fullfile(runDirectory, ...
        'average_recovery_metrics.fig');
    dynsc.plot.exportFigurePdf(averageFigure, averagePdf);
    savefig(averageFigure, averageFig);
    close(averageFigure);

    result.figureFiles.averageRecoveryMetricsPdf = averagePdf;
    result.figureFiles.averageRecoveryMetricsFig = averageFig;

    [edgePdf, edgeFig] = saveAverageSimplexFigure( ...
        result.aggregateMetrics, "edge", runDirectory, ...
        showStandardDeviation);
    [trianglePdf, triangleFig] = saveAverageSimplexFigure( ...
        result.aggregateMetrics, "triangle", runDirectory, ...
        showStandardDeviation);
    result.figureFiles.averageEdgeMetricsPdf = edgePdf;
    result.figureFiles.averageEdgeMetricsFig = edgeFig;
    result.figureFiles.averageTriangleMetricsPdf = trianglePdf;
    result.figureFiles.averageTriangleMetricsFig = triangleFig;

    aggregateMetrics = result.aggregateMetrics;
    save(fullfile(runDirectory, 'aggregate_metrics.mat'), ...
        'aggregateMetrics');
end
save(fullfile(runDirectory, 'result.mat'), 'result', 'cfg', '-v7.3');

fprintf('Finished. Results saved in:\n%s\n', runDirectory);
end

function [pdfPath, figPath] = saveAverageSimplexFigure( ...
    aggregateMetrics, simplexLevel, runDirectory, showStandardDeviation)
figureHandle = dynsc.plot.plotAverageSimplexMetrics( ...
    aggregateMetrics, simplexLevel, false, showStandardDeviation);
baseName = "average_" + simplexLevel + "_metrics";
pdfPath = fullfile(runDirectory, baseName + ".pdf");
figPath = fullfile(runDirectory, baseName + ".fig");
dynsc.plot.exportFigurePdf(figureHandle, pdfPath);
savefig(figureHandle, figPath);
close(figureHandle);
end

function showStandardDeviation = averageStandardDeviationSetting(cfg)
% Older configurations retain the historical shaded-band behavior.
showStandardDeviation = true;
if isfield(cfg, 'plot') && isfield(cfg.plot, 'average') && ...
        isfield(cfg.plot.average, 'showStandardDeviation')
    showStandardDeviation = ...
        logical(cfg.plot.average.showStandardDeviation);
end
end

function enabled = staticBaselineEnabled(cfg)
% Older configurations predate this switch and retain the former behavior.
enabled = true;
if isfield(cfg, 'baseline') && isfield(cfg.baseline, 'static') && ...
        isfield(cfg.baseline.static, 'enabled')
    enabled = logical(cfg.baseline.static.enabled);
end
end

function enabled = identityTransitionBaselineEnabled(cfg)
% Keep the historical config field for backward compatibility. The baseline
% is the identity-mean, fixed-Q persistence prior.
enabled = isfield(cfg, 'baseline') && ...
    isfield(cfg.baseline, 'identityTransition') && ...
    isfield(cfg.baseline.identityTransition, 'enabled') && ...
    logical(cfg.baseline.identityTransition.enabled);
end

function enabled = constrainedRlsBaselineEnabled(cfg)
enabled = isfield(cfg, 'baseline') && ...
    isfield(cfg.baseline, 'constrainedRls') && ...
    isfield(cfg.baseline.constrainedRls, 'enabled') && ...
    logical(cfg.baseline.constrainedRls.enabled);
end

function enabled = previousUpdateBaselineEnabled(cfg)
enabled = isfield(cfg, 'baseline') && ...
    isfield(cfg.baseline, 'previousUpdate') && ...
    isfield(cfg.baseline.previousUpdate, 'enabled') && ...
    logical(cfg.baseline.previousUpdate.enabled) && ...
    string(cfg.update.method) ~= "constrained_map";
end

function enabled = exactHessianBaselineEnabled(cfg)
enabled = isfield(cfg, 'baseline') && ...
    isfield(cfg.baseline, 'exactHessianUpdate') && ...
    isfield(cfg.baseline.exactHessianUpdate, 'enabled') && ...
    logical(cfg.baseline.exactHessianUpdate.enabled) && ...
    string(cfg.update.method) ~= "barrier_exact_hessian";
end

function enabled = standardEkfBarrierBaselineEnabled(cfg)
enabled = isfield(cfg, 'baseline') && ...
    isfield(cfg.baseline, 'standardEkfBarrier') && ...
    isfield(cfg.baseline.standardEkfBarrier, 'enabled') && ...
    logical(cfg.baseline.standardEkfBarrier.enabled) && ...
    string(cfg.update.method) ~= "barrier_standard_ekf_projection";
end

function metrics = selectStaticMetrics(staticBaseline)
if isempty(staticBaseline)
    metrics = [];
else
    metrics = staticBaseline.metrics;
end
end

function result = selectStaticResult(staticBaseline)
if isempty(staticBaseline)
    result = [];
else
    result = staticBaseline.result;
end
end

function metrics = selectIdentityMetrics(identityBaseline)
if isempty(identityBaseline)
    metrics = [];
else
    metrics = identityBaseline.metrics;
end
end

function result = selectIdentityResult(identityBaseline)
if isempty(identityBaseline)
    result = [];
else
    result = identityBaseline.result;
end
end

function metrics = selectConstrainedRlsMetrics(constrainedRlsBaseline)
if isempty(constrainedRlsBaseline)
    metrics = [];
else
    metrics = constrainedRlsBaseline.metrics;
end
end

function result = selectConstrainedRlsResult(constrainedRlsBaseline)
if isempty(constrainedRlsBaseline)
    result = [];
else
    result = constrainedRlsBaseline.result;
end
end

function metrics = selectPreviousUpdateMetrics(previousUpdateBaseline)
if isempty(previousUpdateBaseline)
    metrics = [];
else
    metrics = previousUpdateBaseline.metrics;
end
end

function result = selectPreviousUpdateResult(previousUpdateBaseline)
if isempty(previousUpdateBaseline)
    result = [];
else
    result = previousUpdateBaseline.result;
end
end


function metrics = selectExactHessianMetrics(exactHessianBaseline)
if isempty(exactHessianBaseline)
    metrics = [];
else
    metrics = exactHessianBaseline.metrics;
end
end

function result = selectExactHessianResult(exactHessianBaseline)
if isempty(exactHessianBaseline)
    result = [];
else
    result = exactHessianBaseline.result;
end
end

function metrics = selectStandardEkfMetrics(standardEkfBaseline)
if isempty(standardEkfBaseline)
    metrics = [];
else
    metrics = standardEkfBaseline.metrics;
end
end

function result = selectStandardEkfResult(standardEkfBaseline)
if isempty(standardEkfBaseline)
    result = [];
else
    result = standardEkfBaseline.result;
end
end

function [enabled, interval] = stateSnapshotSettings(cfg)
% Older configurations inherit the requested five-iteration default.
enabled = true;
interval = 5;

if isfield(cfg, 'plot') && isfield(cfg.plot, 'stateSnapshots')
    enabled = cfg.plot.stateSnapshots.enabled;
    interval = cfg.plot.stateSnapshots.interval;
end
end

function runDirectory = createRunDirectory(cfg)
if ~isfolder(cfg.output.root)
    [created, message] = mkdir(cfg.output.root);
    assert(created, 'dynsc:OutputDirectory', ...
        'Could not create output root: %s', message);
end

safeName = matlab.lang.makeValidName(char(cfg.experiment.name));
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
baseDirectory = fullfile(cfg.output.root, ...
    sprintf('%s_%s', safeName, timestamp));
runDirectory = baseDirectory;
suffix = 1;

while isfolder(runDirectory)
    runDirectory = sprintf('%s_%02d', baseDirectory, suffix);
    suffix = suffix + 1;
end

[created, message] = mkdir(runDirectory);
assert(created, 'dynsc:OutputDirectory', ...
    'Could not create run directory: %s', message);
end
