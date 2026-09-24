function fig = plotAverageSimplexMetrics( ...
    aggregate, simplexLevel, visible, showStandardDeviation)
%PLOTAVERAGESIMPLEXMETRICS Plot one simplex level across realizations.

if nargin < 3
    visible = true;
end
if nargin < 4
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

simplexLevel = lower(string(simplexLevel));
assert(isscalar(simplexLevel) && ...
    any(simplexLevel == ["edge", "triangle"]), ...
    'dynsc:InvalidSimplexPlotLevel', ...
    'simplexLevel must be "edge" or "triangle".');
assert(isfield(aggregate, 'time') && ...
    isfield(aggregate, 'numRealizations') && ...
    isfield(aggregate, 'supportThreshold'), ...
    'dynsc:InvalidAggregateMetrics', ...
    'The aggregate metric structure is incomplete.');

if simplexLevel == "edge"
    simplexField = 'edge';
    simplexLabel = "Edge";
else
    simplexField = 'triangle';
    simplexLabel = "Triangle";
end

iterations = aggregate.time;
methods = methodSpecifications(aggregate);
if simplexLevel == "edge"
    % At single-column scale, MATLAB's broken line patterns are too sparse
    % for these slowly varying curves. Color and markers identify methods.
    for iMethod = 1:numel(methods)
        methods(iMethod).lineStyle = '-';
    end
end
fontSize = 9;
legendFontSize = 7;
fontName = 'Times New Roman';

fig = figure('Color', 'w', 'Visible', visibility, ...
    'Units', 'inches', 'Position', [1, 1, 3.45, 3.8]);
layout = tiledlayout(fig, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

f1Axis = nexttile(layout);
hold(f1Axis, 'on');
legendHandles = gobjects(numel(methods), 1);
for iMethod = 1:numel(methods)
    legendHandles(iMethod) = plotEstimator(f1Axis, iterations, ...
        methods(iMethod).data.(simplexField).supportF1, ...
        methods(iMethod), [0, 1], showStandardDeviation);
end
formatAxis(f1Axis, 'Support F1 score');
ylim(f1Axis, [0, 1]);
title(f1Axis, simplexLabel + " support recovery");

errorAxis = nexttile(layout);
hold(errorAxis, 'on');
for iMethod = 1:numel(methods)
    plotEstimator(errorAxis, iterations, ...
        methods(iMethod).data.(simplexField).relativeError, ...
        methods(iMethod), [0, Inf], showStandardDeviation);
end
formatAxis(errorAxis, 'Relative error');
xlabel(errorAxis, 'Time step, {\it k}', 'Interpreter', 'tex');
title(errorAxis, simplexLabel + " Unthresholded Relative Error");

linkaxes([f1Axis, errorAxis], 'x');
xlim(f1Axis, [1, max(2, numel(iterations))]);
legendHandle = addSharedLegend(f1Axis, legendHandles, methods);
set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);
set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
legendHandle.FontSize = legendFontSize;
legendHandle.FontName = fontName;
end

function legendHandle = addSharedLegend( ...
    axisHandle, legendHandles, methods)
legendLabels = strings(numel(methods), 1);
for iMethod = 1:numel(methods)
    legendLabels(iMethod) = methods(iMethod).label;
end
legendHandle = legend(axisHandle, legendHandles, legendLabels, ...
    'Orientation', 'horizontal', ...
    'NumColumns', 3, ...
    'Location', 'southoutside');
legendHandle.Box = 'off';
legendHandle.ItemTokenSize = [12, 6];
legendHandle.Layout.Tile = 'south';
end

function methods = methodSpecifications(aggregate)
% Display order: Static, Persistence, Constrained RLS, Constrained MAP, ...
% EKF, IEKF, Laplace.
methods = struct('data', {}, 'label', {}, 'color', {}, ...
    'lineStyle', {}, 'marker', {});
colors = {
    [0.30, 0.30, 0.30]
    [0.4940, 0.1840, 0.5560]
    [0.1500, 0.1500, 0.1500]
    [0.6350, 0.0780, 0.1840]
    [0.5490, 0.3373, 0.2941]
    [0.1216, 0.4667, 0.7059]
    [0.0000, 0.6000, 0.6000]
};

methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'static', 'static'), ...
    "Static", colors{1}, ':', 's');
[identityEstimator, identityLabel] = resolveIdentityEstimator(aggregate);
methods = appendMethod(methods, identityEstimator, ...
    identityLabel, colors{2}, '--', 'd');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'constrainedRls'), ...
    "Constrained RLS", colors{3}, ':', 'x');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'constrainedMap', 'previousUpdate'), ...
    "Constrained MAP", colors{4}, '-.', 'p');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'ekf', 'standardEkfBarrier'), ...
    "EKF", colors{5}, '--', '^');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'iekf', 'primary'), ...
    "IEKF", colors{6}, '-', 'o');
methods = appendMethod(methods, ...
    resolveEstimator(aggregate, 'laplace', 'exactHessian'), ...
    "Laplace", colors{7}, '-', 'h');
end

function methods = appendMethod( ...
    methods, data, label, color, lineStyle, marker)
if isempty(data)
    return;
end
entry.data = data;
entry.label = label;
entry.color = color;
entry.lineStyle = lineStyle;
entry.marker = marker;
methods(end + 1) = entry;
end

function [estimator, label] = resolveIdentityEstimator(aggregate)
estimator = resolveEstimator(aggregate, 'persistence', 'identity');
label = estimatorLabel(estimator, "Persistence");
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

function lineHandle = plotEstimator( ...
    axisHandle, iterations, summary, method, limits, ...
    showStandardDeviation)
if showStandardDeviation
    lower = max(summary.mean - summary.std, limits(1));
    upper = summary.mean + summary.std;
    if isfinite(limits(2))
        upper = min(upper, limits(2));
    end
    fill(axisHandle, [iterations, fliplr(iterations)], ...
        [lower, fliplr(upper)], method.color, ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
lineHandle = plot(axisHandle, iterations, summary.mean, ...
    'Color', method.color, 'LineStyle', method.lineStyle, ...
    'Marker', method.marker, ...
    'MarkerIndices', markerIndices(iterations), ...
    'MarkerSize', 2.4, 'MarkerFaceColor', 'none', ...
    'LineWidth', 0.75, ...
    'DisplayName', method.label);
end

function indices = markerIndices(iterations)
indices = unique([1:4:numel(iterations), numel(iterations)]);
end

function formatAxis(axisHandle, yLabel)
ylabel(axisHandle, yLabel);
grid(axisHandle, 'on');
box(axisHandle, 'on');
axisHandle.FontSize = 12;
end
