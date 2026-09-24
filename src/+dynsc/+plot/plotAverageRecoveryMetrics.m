function fig = plotAverageRecoveryMetrics( ...
    aggregate, visible, showStandardDeviation)
%PLOTAVERAGERECOVERYMETRICS Plot metrics averaged over realizations.

if nargin < 2
    visible = true;
end
if nargin < 3
    showStandardDeviation = true;
end
validateattributes(showStandardDeviation, ...
    {'logical', 'numeric'}, {'scalar'});
showStandardDeviation = logical(showStandardDeviation);

if visible
    visibility = 'on';
else
    visibility = 'off';
end

validateAggregate(aggregate);
iterations = aggregate.time;
methods = methodSpecifications(aggregate);
fontSize = 8;
legendFontSize = 6.5;
fontName = 'Times New Roman';

fig = figure('Color', 'w', 'Visible', visibility, ...
    'Units', 'inches', 'Position', [1, 1, 3.45, 4.05]);
layout = tiledlayout(fig, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

f1Axis = nexttile(layout);
hold(f1Axis, 'on');
for iMethod = 1:numel(methods)
    plotEstimatorPair(f1Axis, iterations, methods(iMethod), ...
        'supportF1', [0, 1], showStandardDeviation);
end
formatAxis(f1Axis, 'Support F1 score');
ylim(f1Axis, [0, 1]);
title(f1Axis, 'Support recovery');

errorAxis = nexttile(layout);
hold(errorAxis, 'on');
for iMethod = 1:numel(methods)
    plotEstimatorPair(errorAxis, iterations, methods(iMethod), ...
        'relativeError', [0, Inf], showStandardDeviation);
end
formatAxis(errorAxis, 'Relative error');
xlabel(errorAxis, 'Time step, {\it k}', 'Interpreter', 'tex');
title(errorAxis, 'Unthresholded Topology Error');

linkaxes([f1Axis, errorAxis], 'x');
xlim(f1Axis, [1, max(2, numel(iterations))]);
legendHandle = addSharedLegend(f1Axis, methods);
set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);
set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
legendHandle.FontSize = legendFontSize;
legendHandle.FontName = fontName;
end

function legendHandle = addSharedLegend(axisHandle, methods)
% Use color for method and line style for simplex level.
numMethods = numel(methods);
legendHandles = gobjects(numMethods + 2, 1);
legendLabels = strings(numMethods + 2, 1);
for iMethod = 1:numMethods
    legendHandles(iMethod) = plot(axisHandle, nan, nan, ...
        'Color', methods(iMethod).color, 'LineStyle', '-', ...
        'LineWidth', 0.8, ...
        'Marker', methods(iMethod).marker, 'MarkerSize', 2.8, ...
        'MarkerFaceColor', 'none');
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
    'Orientation', 'horizontal', ...
    'NumColumns', 4, ...
    'Location', 'southoutside');
legendHandle.Box = 'off';
legendHandle.ItemTokenSize = [12, 6];
legendHandle.Layout.Tile = 'south';
end

function methods = methodSpecifications(aggregate)
% Display order: Static, Persistence, Constrained RLS, Constrained MAP, ...
% EKF, IEKF, Laplace.
methods = struct('data', {}, 'label', {}, ...
    'color', {}, ...
    'edgeStyle', {}, 'triangleStyle', {}, ...
    'marker', {});

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

function methods = appendMethod( ...
    methods, data, label, color, marker)
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
% In the short-lived legacy RLS format, persistence was an alias for the
% same constrainedRls summary. Avoid plotting that saved curve twice.
if label == "Constrained RLS" && ...
        isfield(aggregate, 'constrainedRls') && ...
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
    if isfield(aggregate, fieldName) && ...
            ~isempty(aggregate.(fieldName))
        estimator = aggregate.(fieldName);
        return;
    end
end
end

function plotEstimatorPair( ...
    axisHandle, iterations, method, metricName, limits, ...
    showStandardDeviation)
plotMeanAndBand(axisHandle, iterations, ...
    method.data.edge.(metricName), method.color, ...
    method.edgeStyle, method.marker, method.label + " edges", ...
    limits, showStandardDeviation);
plotMeanAndBand(axisHandle, iterations, ...
    method.data.triangle.(metricName), method.color, ...
    method.triangleStyle, method.marker, ...
    method.label + " triangles", limits, showStandardDeviation);
end

function lineHandle = plotMeanAndBand( ...
    axisHandle, iterations, summary, color, lineStyle, marker, ...
    displayName, limits, showStandardDeviation)
if showStandardDeviation
    lower = summary.mean - summary.std;
    upper = summary.mean + summary.std;
    lower = max(lower, limits(1));
    if isfinite(limits(2))
        upper = min(upper, limits(2));
    end
    fill(axisHandle, [iterations, fliplr(iterations)], ...
        [lower, fliplr(upper)], color, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
lineHandle = plot(axisHandle, iterations, summary.mean, ...
    'Color', color, 'LineStyle', lineStyle, ...
    'Marker', marker, 'MarkerIndices', markerIndices(iterations), ...
    'MarkerSize', 2.4, 'MarkerFaceColor', 'none', ...
    'LineWidth', 0.75, ...
    'DisplayName', displayName);
end

function indices = markerIndices(iterations)
indices = unique([1:4:numel(iterations), numel(iterations)]);
end

function validateAggregate(aggregate)
requiredFields = ["numRealizations", "time", "supportThreshold"];
for iField = 1:numel(requiredFields)
    assert(isfield(aggregate, requiredFields(iField)), ...
        'dynsc:InvalidAggregateMetrics', ...
        'Missing aggregate metric field: %s.', requiredFields(iField));
end
assert((isfield(aggregate, 'iekf') && ~isempty(aggregate.iekf)) || ...
    (isfield(aggregate, 'primary') && ~isempty(aggregate.primary)), ...
    'dynsc:InvalidAggregateMetrics', ...
    'The aggregate metrics must contain IEKF results.');
validateattributes(aggregate.numRealizations, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
end

function formatAxis(axisHandle, yLabel)
ylabel(axisHandle, yLabel);
grid(axisHandle, 'on');
box(axisHandle, 'on');
axisHandle.FontSize = 12;
end
