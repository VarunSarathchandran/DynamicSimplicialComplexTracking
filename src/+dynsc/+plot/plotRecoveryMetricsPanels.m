function fig = plotRecoveryMetricsPanels(runs, varargin)
%PLOTRECOVERYMETRICSPANELS Side-by-side recovery panels with one legend.
%
%   fig = dynsc.plot.plotRecoveryMetricsPanels(runs) draws one column per
%   entry of runs, with support F1 on the top row and unthresholded
%   relative error on the bottom row, and a single legend shared by every
%   panel. The line styling, colours, markers, fonts and metric choices are
%   the same as dynsc.plot.plotAverageRecoveryMetrics, so a panel here is
%   visually interchangeable with the single-run figures already used in
%   the paper.
%
%   runs is a cell array whose entries may each be
%     - an aggregate metrics structure,
%     - a result structure from a tracking runner,
%     - a path to aggregate_metrics.mat or result.mat, or
%     - a run directory containing either of those files.
%
%   Name-value arguments:
%     'Labels'          string array, one column title per run.
%     'Phase'           "test" (default) or "calibration"; used only when a
%                       run is supplied as a runner result structure.
%     'ShowStandardDeviation'
%                       false (default). Single-realization runs have no
%                       spread, so bands are off unless asked for.
%     'ShareRelativeErrorAxis'
%                       false (default). true gives every column the same
%                       relative-error limits, which is only meaningful
%                       when the runs are comparable in scale.
%     'RelativeErrorLimit'
%                       [lo hi] forced on every column; implies sharing.
%     'FigureSize'      [width height] in inches, default [6.9, 4.05].
%     'Visible'         false (default) or true.
%     'SavePath'        PDF path to write; omitted means nothing is written.
%
%   No estimation is repeated: every curve comes from stored metrics.

parser = inputParser;
parser.addRequired('runs', @(x) iscell(x) && ~isempty(x));
parser.addParameter('Labels', strings(0), @(x) isstring(x) || iscellstr(x)); %#ok<ISCLSTR>
parser.addParameter('Phase', "test", @(x) any(string(x) == ["test", "calibration"]));
parser.addParameter('Methods', strings(0), @(x) isstring(x) || iscellstr(x)); %#ok<ISCLSTR>
parser.addParameter('Metrics', ["supportF1", "relativeError"], ...
    @(x) (isstring(x) || iscellstr(x)) && ~isempty(x)); %#ok<ISCLSTR>
% Row titles. Narrow multi-column layouts need shorter ones than the
% single-run figures use, so they are overridable.
parser.addParameter('MetricTitles', strings(0), ...
    @(x) isstring(x) || iscellstr(x)); %#ok<ISCLSTR>
parser.addParameter('ShowStandardDeviation', false, @(x) islogical(x) || isnumeric(x));
parser.addParameter('ShareRelativeErrorAxis', false, @(x) islogical(x) || isnumeric(x));
parser.addParameter('RelativeErrorLimit', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && x(2) > x(1)));
% "grid" is metrics-as-rows by runs-as-columns. "row" flattens every panel
% into a single line, which is much shorter on the page; panels of one run
% stay adjacent and the run label is drawn once above each pair.
parser.addParameter('Layout', "grid", ...
    @(x) any(string(x) == ["grid", "row"]));
% "none" drops the axis label text, for when the panel title already names
% the metric and the abscissa. Tick labels are always kept.
parser.addParameter('AxisLabels', "both", ...
    @(x) any(string(x) == ["both", "none"]));
parser.addParameter('LegendColumns', 4, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
% IEEE requires figure text to be no smaller than the body text, so the
% default matches the axis font rather than shrinking the legend.
parser.addParameter('LegendFontSize', 8, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('FigureSize', [6.9, 4.05], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
parser.addParameter('Visible', false, @(x) islogical(x) || isnumeric(x));
parser.addParameter('SavePath', "", @(x) isstring(x) || ischar(x));
parser.parse(runs, varargin{:});
opts = parser.Results;

phase = string(opts.Phase);
showBand = logical(opts.ShowStandardDeviation);
keepMethods = string(opts.Methods);
metrics = string(opts.Metrics);
for iMetric = 1:numel(metrics)
    assert(any(metrics(iMetric) == ["supportF1", "relativeError"]), ...
        'dynsc:UnknownMetric', ...
        'Metrics must be supportF1 and/or relativeError.');
end
nRow = numel(metrics);
layoutMode = string(opts.Layout);
showAxisLabels = string(opts.AxisLabels) == "both";
metricTitles = string(opts.MetricTitles);
if isempty(metricTitles)
    metricTitles = arrayfun(@defaultMetricTitle, metrics);
end
assert(numel(metricTitles) == nRow, 'dynsc:MetricTitleCountMismatch', ...
    'One metric title is required per metric row.');
nColumn = numel(runs);
labels = string(opts.Labels);
if isempty(labels)
    labels = "run" + string(1:nColumn);
end
assert(numel(labels) == nColumn, 'dynsc:LabelCountMismatch', ...
    'One label is required per run.');

aggregates = cell(1, nColumn);
for iColumn = 1:nColumn
    aggregates{iColumn} = loadAggregate(runs{iColumn}, phase);
end

% ---- styling, copied from plotAverageRecoveryMetrics --------------------
fontSize = 8;
legendFontSize = opts.LegendFontSize;
fontName = 'Times New Roman';

if logical(opts.Visible)
    visibility = 'on';
else
    visibility = 'off';
end

fig = figure('Color', 'w', 'Visible', visibility, ...
    'Units', 'inches', ...
    'Position', [1, 1, opts.FigureSize(1), opts.FigureSize(2)]);
if layoutMode == "row"
    layout = tiledlayout(fig, 1, nRow * nColumn, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    layout.TileSpacing = 'tight';
    layout.Padding = 'tight';
    runLabelBand = 0.085;
    layout.OuterPosition = [0, 0, 1, 1 - runLabelBand];
else
    layout = tiledlayout(fig, nRow, nColumn, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
end

axesGrid = gobjects(nRow, nColumn);
columnMethods = cell(1, nColumn);
for iColumn = 1:nColumn
    columnMethods{iColumn} = ...
        selectMethods(methodSpecifications(aggregates{iColumn}), keepMethods);
    assert(~isempty(columnMethods{iColumn}), 'dynsc:NoMethodsSelected', ...
        'No requested method is present in run %d.', iColumn);
end

for iRow = 1:nRow
    metricName = metrics(iRow);
    for iColumn = 1:nColumn
        aggregate = aggregates{iColumn};
        methods = columnMethods{iColumn};
        if layoutMode == "row"
            % Keep a run's metrics next to each other.
            tileIndex = (iColumn - 1) * nRow + iRow;
        else
            tileIndex = (iRow - 1) * nColumn + iColumn;
        end
        axisHandle = nexttile(layout, tileIndex);
        hold(axisHandle, 'on');
        for iMethod = 1:numel(methods)
            plotEstimatorPair(axisHandle, aggregate.time, ...
                methods(iMethod), metricName, metricLimits(metricName), ...
                showBand);
        end
        formatAxis(axisHandle);
        xlim(axisHandle, [1, max(2, numel(aggregate.time))]);
        firstOfKind = (layoutMode == "row") || (iColumn == 1);
        if metricName == "supportF1"
            ylim(axisHandle, [0, 1]);
        end
        if showAxisLabels && firstOfKind
            if metricName == "supportF1"
                ylabel(axisHandle, 'Support F1 score');
            else
                ylabel(axisHandle, 'Relative error');
            end
        elseif ~firstOfKind
            axisHandle.YTickLabel = [];
        end
        if layoutMode == "row"
            % The run label is drawn once per pair, after layout settles.
            title(axisHandle, metricTitles(iRow), 'FontWeight', 'normal');
            if showAxisLabels
                xlabel(axisHandle, '{\it k}', 'Interpreter', 'tex');
            end
        elseif iRow == 1
            title(axisHandle, labels(iColumn));
            subtitle(axisHandle, metricTitles(iRow));
            axisHandle.Title.FontWeight = 'bold';
            axisHandle.Subtitle.FontWeight = 'normal';
        else
            title(axisHandle, metricTitles(iRow), ...
                'FontWeight', 'normal');
        end
        if layoutMode ~= "row"
            if iRow == nRow && showAxisLabels
                xlabel(axisHandle, 'Time step, {\it k}', 'Interpreter', 'tex');
            elseif iRow ~= nRow
                axisHandle.XTickLabel = [];
            end
        end
        axesGrid(iRow, iColumn) = axisHandle;
    end
end

errorRow = find(metrics == "relativeError", 1);
if ~isempty(errorRow)
    errorAxes = axesGrid(errorRow, :);
else
    errorAxes = gobjects(1, 0);
end

applyRelativeErrorLimits(errorAxes, aggregates, columnMethods, opts, layoutMode);

% ---- one legend for the whole figure ------------------------------------
legendHandle = addSharedLegend(axesGrid(1, 1), ...
    unionMethods(columnMethods), opts.LegendColumns);

set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);
set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
legendHandle.FontSize = legendFontSize;
legendHandle.FontName = fontName;
if layoutMode == "row"
    annotateRunGroups(fig, axesGrid, labels, fontName, fontSize);
else
    for iColumn = 1:nColumn
        axesGrid(1, iColumn).Title.FontWeight = 'bold';
        if ~isempty(axesGrid(1, iColumn).Subtitle.String)
            axesGrid(1, iColumn).Subtitle.FontWeight = 'normal';
        end
    end
end

if strlength(string(opts.SavePath)) > 0
    savePath = char(string(opts.SavePath));
    folder = fileparts(savePath);
    if ~isempty(folder) && ~isfolder(folder)
        mkdir(folder);
    end
    dynsc.plot.exportFigurePdf(fig, savePath);
    fprintf('Wrote %s\n', savePath);
end
end

% ------------------------------------------------------------------------

function applyRelativeErrorLimits( ...
    errorAxes, aggregates, columnMethods, opts, layoutMode)
if isempty(errorAxes)
    return;
end
share = logical(opts.ShareRelativeErrorAxis) || ~isempty(opts.RelativeErrorLimit);
if ~share
    for iAxis = 1:numel(errorAxes)
        limit = columnErrorLimit(aggregates{iAxis}, columnMethods{iAxis});
        ylim(errorAxes(iAxis), limit);
    end
    return;
end
if ~isempty(opts.RelativeErrorLimit)
    limit = opts.RelativeErrorLimit;
else
    limit = [0, 0];
    for iAxis = 1:numel(errorAxes)
        candidate = columnErrorLimit(aggregates{iAxis}, columnMethods{iAxis});
        limit(2) = max(limit(2), candidate(2));
    end
end
for iAxis = 1:numel(errorAxes)
    ylim(errorAxes(iAxis), limit);
    if iAxis > 1 && layoutMode ~= "row"
        errorAxes(iAxis).YTickLabel = [];
    end
end
fprintf('Shared relative-error limit: [%g %g]\n', limit);
end

function annotateRunGroups(fig, axesGrid, labels, fontName, fontSize)
% Centre each run's label over the panels belonging to that run. Positions
% are read back from the laid-out axes so the label tracks the real tiles.
drawnow;
nColumn = size(axesGrid, 2);
for iColumn = 1:nColumn
    positions = vertcat(axesGrid(:, iColumn).Position);
    left = min(positions(:, 1));
    right = max(positions(:, 1) + positions(:, 3));
    annotation(fig, 'textbox', ...
        [left, 1 - 0.075, right - left, 0.06], ...
        'String', labels(iColumn), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'EdgeColor', 'none', ...
        'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'bold');
end
end

function limits = metricLimits(metricName)
if metricName == "supportF1"
    limits = [0, 1];
else
    limits = [0, Inf];
end
end

function name = defaultMetricTitle(metricName)
if metricName == "supportF1"
    name = "Support recovery";
else
    name = "Unthresholded Topology Error";
end
end

function methods = selectMethods(methods, keepMethods)
% Keep only the requested estimator labels, in the canonical display order.
if isempty(keepMethods)
    return;
end
methods = methods(ismember(string({methods.label}), keepMethods));
end

function limit = columnErrorLimit(aggregate, methods)
largest = 0;
for iMethod = 1:numel(methods)
    for level = ["edge", "triangle"]
        series = methods(iMethod).data.(level).relativeError.mean;
        series = series(isfinite(series));
        if ~isempty(series)
            largest = max(largest, max(series));
        end
    end
end
if ~isfinite(largest) || largest <= 0
    limit = [0, 1];
else
    limit = [0, ceil(largest * 10) / 10];
end
end

function aggregate = loadAggregate(runSpecification, phase)
if isstruct(runSpecification)
    aggregate = structToAggregate(runSpecification, phase);
    return;
end
path = char(string(runSpecification));
if isfolder(path)
    candidates = {fullfile(path, 'aggregate_metrics.mat'), ...
        fullfile(path, 'result.mat')};
    path = '';
    for iCandidate = 1:numel(candidates)
        if isfile(candidates{iCandidate})
            path = candidates{iCandidate};
            break;
        end
    end
    assert(~isempty(path), 'dynsc:MissingRunResult', ...
        'Neither aggregate_metrics.mat nor result.mat was found in %s.', ...
        char(string(runSpecification)));
end
assert(isfile(path), 'dynsc:MissingRunResult', 'No file at %s', path);
loaded = load(path);
names = fieldnames(loaded);
for iName = 1:numel(names)
    candidate = loaded.(names{iName});
    if isstruct(candidate) && isscalar(candidate) && ...
            (isfield(candidate, 'realizations') || isfield(candidate, 'time'))
        aggregate = structToAggregate(candidate, phase);
        return;
    end
end
error('dynsc:UnrecognisedMetricsFile', ...
    'No aggregate or result structure was found in %s.', path);
end

function aggregate = structToAggregate(candidate, phase)
if isfield(candidate, 'time') && isfield(candidate, 'numRealizations')
    aggregate = candidate;
    return;
end
assert(isfield(candidate, 'realizations'), 'dynsc:UnrecognisedMetrics', ...
    'A run must be an aggregate or a runner result structure.');
realization = candidate.realizations{1};
if phase == "test"
    methods = realization.testMethods;
else
    methods = realization.calibrationMethods;
end
aggregate = singleRealizationAggregate(methods);
end

function aggregate = singleRealizationAggregate(methods)
% Mirrors the private helper in runPrimarySchoolTrackingExperiment so that
% a saved run can be replotted outside the runner.
realization = struct();
assert(~isempty(methods.iekf), 'dynsc:MissingIekfForPlot', ...
    'The paper plotting pipeline requires IEKF metrics as its reference.');
realization.metrics = methods.iekf.metrics;
if isfield(methods, 'constrainedRls') && ~isempty(methods.constrainedRls)
    realization.constrainedRlsBaseline = ...
        struct('metrics', methods.constrainedRls.metrics);
end
if isfield(methods, 'constrainedMap') && ~isempty(methods.constrainedMap)
    realization.previousUpdateBaseline = ...
        struct('metrics', methods.constrainedMap.metrics);
end
if isfield(methods, 'ekf') && ~isempty(methods.ekf)
    realization.standardEkfBarrierBaseline = ...
        struct('metrics', methods.ekf.metrics);
end
if isfield(methods, 'laplace') && ~isempty(methods.laplace)
    realization.exactHessianBaseline = ...
        struct('metrics', methods.laplace.metrics);
end
aggregate = dynsc.metrics.aggregateRecoveryMetrics({realization});
end

% ---- styling helpers, kept identical to plotAverageRecoveryMetrics -----

function methods = methodSpecifications(aggregate)
methods = struct('data', {}, 'label', {}, 'color', {}, ...
    'edgeStyle', {}, 'triangleStyle', {}, 'marker', {});

methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'static', 'static'), "Static", ...
    [0.30, 0.30, 0.30], 's');
[identityEstimator, identityLabel] = resolveIdentityEstimator(aggregate);
methods = appendMethod(methods, identityEstimator, ...
    identityLabel, [0.4940, 0.1840, 0.5560], 'd');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'constrainedRls'), ...
    "Constrained RLS", [0.1500, 0.1500, 0.1500], 'x');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'constrainedMap', 'previousUpdate'), ...
    "Constrained MAP", [0.6350, 0.0780, 0.1840], 'p');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'ekf', 'standardEkfBarrier'), "EKF", ...
    [0.5490, 0.3373, 0.2941], '^');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'iekf', 'primary'), "IEKF", ...
    [0.1216, 0.4667, 0.7059], 'o');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'laplace', 'exactHessian'), "Laplace", ...
    [0.0000, 0.6000, 0.6000], 'h');
end

function methods = unionMethods(columnMethods)
% Legend entries in the canonical display order, over every method that
% appears in at least one panel.
methods = columnMethods{1};
for iColumn = 2:numel(columnMethods)
    candidates = columnMethods{iColumn};
    for iMethod = 1:numel(candidates)
        if ~any(string({methods.label}) == candidates(iMethod).label)
            methods(end + 1) = candidates(iMethod); %#ok<AGROW>
        end
    end
end
order = ["Static", "Persistence", "Constrained RLS", "Constrained MAP", ...
    "EKF", "IEKF", "Laplace"];
present = string({methods.label});
[~, rank] = ismember(present, order);
rank(rank == 0) = numel(order) + 1;
[~, sorted] = sort(rank);
methods = methods(sorted);
end

function methods = appendMethod(methods, data, label, color, marker)
if isempty(data)
    return;
end
entry.data = data;
entry.label = label;
entry.color = color;
entry.edgeStyle = '-';
entry.triangleStyle = ':';
entry.marker = marker;
methods(end + 1) = entry;
end

function [estimator, label] = resolveIdentityEstimator(aggregate)
estimator = resolveEstimator(aggregate, 'persistence', 'identity');
label = estimatorLabel(estimator, "Persistence");
if label == "Constrained RLS" && isfield(aggregate, 'constrainedRls') && ...
        ~isempty(aggregate.constrainedRls)
    estimator = [];
    label = "Persistence";
end
end

function label = estimatorLabel(estimator, fallback)
label = fallback;
if ~isempty(estimator) && isfield(estimator, 'name') && ...
        strlength(string(estimator.name)) > 0
    label = string(estimator.name);
end
end

function estimator = resolveEstimator(aggregate, varargin)
estimator = [];
for iName = 1:numel(varargin)
    fieldName = varargin{iName};
    if isfield(aggregate, fieldName) && ~isempty(aggregate.(fieldName))
        estimator = aggregate.(fieldName);
        return;
    end
end
end

function plotEstimatorPair( ...
    axisHandle, iterations, method, metricName, limits, showBand)
plotMeanAndBand(axisHandle, iterations, ...
    method.data.edge.(metricName), method.color, ...
    method.edgeStyle, method.marker, method.label + " edges", ...
    limits, showBand);
plotMeanAndBand(axisHandle, iterations, ...
    method.data.triangle.(metricName), method.color, ...
    method.triangleStyle, method.marker, method.label + " triangles", ...
    limits, showBand);
end

function plotMeanAndBand(axisHandle, iterations, summary, color, ...
    lineStyle, marker, displayName, limits, showBand)
if showBand
    lower = summary.mean - summary.std;
    upper = summary.mean + summary.std;
    lower = max(lower, limits(1));
    if isfinite(limits(2))
        upper = min(upper, limits(2));
    end
    fill(axisHandle, [iterations, fliplr(iterations)], ...
        [lower, fliplr(upper)], color, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
plot(axisHandle, iterations, summary.mean, ...
    'Color', color, 'LineStyle', lineStyle, ...
    'Marker', marker, 'MarkerIndices', markerIndices(iterations), ...
    'MarkerSize', 2.4, 'MarkerFaceColor', 'none', 'LineWidth', 0.75, ...
    'DisplayName', displayName, 'HandleVisibility', 'off');
end

function indices = markerIndices(iterations)
indices = unique([1:4:numel(iterations), numel(iterations)]);
end

function legendHandle = addSharedLegend(axisHandle, methods, numColumns)
% Colour encodes the method, line style encodes the simplex level.
numMethods = numel(methods);
legendHandles = gobjects(numMethods + 2, 1);
legendLabels = strings(numMethods + 2, 1);
for iMethod = 1:numMethods
    legendHandles(iMethod) = plot(axisHandle, nan, nan, ...
        'Color', methods(iMethod).color, 'LineStyle', '-', ...
        'LineWidth', 0.8, 'Marker', methods(iMethod).marker, ...
        'MarkerSize', 2.8, 'MarkerFaceColor', 'none');
    legendLabels(iMethod) = methods(iMethod).label;
end
legendHandles(end - 1) = plot(axisHandle, nan, nan, ...
    'Color', [0.15, 0.15, 0.15], 'LineStyle', '-', ...
    'LineWidth', 0.8, 'Marker', 'none');
legendLabels(end - 1) = "Edges";
legendHandles(end) = plot(axisHandle, nan, nan, ...
    'Color', [0.15, 0.15, 0.15], 'LineStyle', ':', ...
    'LineWidth', 0.8, 'Marker', 'none');
legendLabels(end) = "Triangles";

legendHandle = legend(axisHandle, legendHandles, legendLabels, ...
    'Orientation', 'horizontal', 'NumColumns', numColumns, ...
    'Location', 'southoutside');
legendHandle.Box = 'off';
legendHandle.ItemTokenSize = [12, 6];
legendHandle.Layout.Tile = 'south';
end

function formatAxis(axisHandle)
grid(axisHandle, 'on');
box(axisHandle, 'on');
axisHandle.FontSize = 12;
end
