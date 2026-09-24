function fig = plotRecoveryMetrics( ...
    metrics, baselineMetrics, visible, identityMetrics, ...
    previousUpdateMetrics, exactHessianMetrics, standardEkfMetrics, ...
    constrainedRlsMetrics)
%PLOTRECOVERYMETRICS Plot one realization using scientific method names.

if nargin < 2 || isempty(baselineMetrics)
    baselineMetrics = [];
end
% Preserve the original plotRecoveryMetrics(metrics, visible) call.
if nargin == 2 && ~isstruct(baselineMetrics)
    visible = logical(baselineMetrics);
    baselineMetrics = [];
elseif nargin < 3
    visible = true;
end
if nargin < 4 || isempty(identityMetrics)
    identityMetrics = [];
end
if nargin < 5 || isempty(previousUpdateMetrics)
    previousUpdateMetrics = [];
end
if nargin < 6 || isempty(exactHessianMetrics)
    exactHessianMetrics = [];
end
if nargin < 7 || isempty(standardEkfMetrics)
    standardEkfMetrics = [];
end
if nargin < 8 || isempty(constrainedRlsMetrics)
    constrainedRlsMetrics = [];
end

if visible
    visibility = 'on';
else
    visibility = 'off';
end

nTime = numel(metrics.time);
iterations = metrics.time;
methods = methodSpecifications( ...
    baselineMetrics, identityMetrics, previousUpdateMetrics, ...
    standardEkfMetrics, metrics, exactHessianMetrics, ...
    constrainedRlsMetrics);
for iMethod = 1:numel(methods)
    validateMetricSequence(methods(iMethod).data, metrics, nTime, ...
        methods(iMethod).label);
end
showNormalizedSHD = hasMetric(methods, 'normalizedSHD');
numRows = 3 + double(showNormalizedSHD);

fig = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [100, 100, 1000, 760 + 190 * showNormalizedSHD]);
layout = tiledlayout(fig, numRows, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

residualAxis = nexttile(layout);
hold(residualAxis, 'on');
for iMethod = 1:numel(methods)
    residual = methods(iMethod).data.priorResidualNorm;
    if isempty(residual)
        continue;
    end
    plot(residualAxis, iterations, residual, ...
        'Color', methods(iMethod).edgeColor, ...
        'LineStyle', methods(iMethod).edgeStyle, ...
        'Marker', methods(iMethod).edgeMarker, ...
        'MarkerFaceColor', 'none', 'LineWidth', 1.8, ...
        'MarkerSize', 6, 'DisplayName', methods(iMethod).label);
end
legend(residualAxis, 'Location', 'best');
formatAxis(residualAxis, 'Prior-residual norm');
title(residualAxis, 'Correction from the predicted state');

f1Axis = nexttile(layout);
hold(f1Axis, 'on');
for iMethod = 1:numel(methods)
    plotEstimatorPair(f1Axis, iterations, methods(iMethod), ...
        'supportF1');
end
formatAxis(f1Axis, 'Support F1 score');
ylim(f1Axis, [0, 1]);
legend(f1Axis, 'Location', 'best');
title(f1Axis, sprintf('Support recovery (threshold = %.3g)', ...
    metrics.supportThreshold));

linkedAxes = [residualAxis, f1Axis];
if showNormalizedSHD
    shdAxis = nexttile(layout);
    hold(shdAxis, 'on');
    for iMethod = 1:numel(methods)
        if ~hasEstimatorMetric(methods(iMethod), 'normalizedSHD')
            continue;
        end
        plotEstimatorPair(shdAxis, iterations, methods(iMethod), ...
            'normalizedSHD');
    end
    formatAxis(shdAxis, 'Normalized SHD');
    ylim(shdAxis, [0, 1]);
    legend(shdAxis, 'Location', 'best');
    title(shdAxis, 'Normalized structural Hamming distance');
    linkedAxes(end + 1) = shdAxis;
end

errorAxis = nexttile(layout);
hold(errorAxis, 'on');
for iMethod = 1:numel(methods)
    plotEstimatorPair(errorAxis, iterations, methods(iMethod), ...
        'relativeError');
end
formatAxis(errorAxis, 'Relative error');
xlabel(errorAxis, 'Filter iteration, k');
legend(errorAxis, 'Location', 'best');
title(errorAxis, 'Unthresholded Topolgy Error');

changeTimes = metrics.topologyChangeTimes;
addChangeMarkers(residualAxis, changeTimes);
addChangeMarkers(f1Axis, changeTimes);
if showNormalizedSHD
    addChangeMarkers(shdAxis, changeTimes);
end
addChangeMarkers(errorAxis, changeTimes);

linkedAxes(end + 1) = errorAxis;
linkaxes(linkedAxes, 'x');
xlim(residualAxis, [1, max(2, nTime)]);
if isempty(changeTimes)
    figureTitle = 'Recovery diagnostics over time (no topology changes)';
else
    changeText = strjoin(string(changeTimes), ', ');
    figureTitle = sprintf( ...
        'Recovery diagnostics over time | Triangle swaps at k = %s', ...
        changeText);
end
title(layout, figureTitle);
end

function available = hasMetric(methods, metricName)
available = false;
for iMethod = 1:numel(methods)
    available = available || ...
        hasEstimatorMetric(methods(iMethod), metricName);
end
end

function available = hasEstimatorMetric(method, metricName)
available = isfield(method.data.edge, metricName) && ...
    ~isempty(method.data.edge.(metricName)) && ...
    isfield(method.data.triangle, metricName) && ...
    ~isempty(method.data.triangle.(metricName));
end

function methods = methodSpecifications( ...
    staticMetrics, identityMetrics, constrainedMapMetrics, ...
    ekfMetrics, iekfMetrics, laplaceMetrics, constrainedRlsMetrics)
% Display order: Static, Persistence, Constrained RLS, Constrained MAP, ...
% EKF, IEKF, Laplace.
methods = struct('data', {}, 'label', {}, ...
    'edgeColor', {}, 'triangleColor', {}, ...
    'edgeStyle', {}, 'triangleStyle', {}, ...
    'edgeMarker', {}, 'triangleMarker', {});
methods = appendMethod(methods, staticMetrics, "Static", ...
    [0.30, 0.30, 0.30], [0.60, 0.60, 0.60], ...
    ':', '-.', 'o', 's');
identityLabel = metricLabel(identityMetrics, "Persistence");
methods = appendMethod(methods, identityMetrics, identityLabel, ...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], ...
    '--', ':', 'd', '^');
methods = appendMethod(methods, constrainedRlsMetrics, ...
    "Constrained RLS", [0.1500, 0.1500, 0.1500], ...
    [0.8500, 0.4500, 0.1000], ':', '-.', 'x', '+');
methods = appendMethod(methods, constrainedMapMetrics, ...
    "Constrained MAP", [0.6350, 0.0780, 0.1840], ...
    [0.3010, 0.7450, 0.9330], '-.', ':', 'p', 'v');
methods = appendMethod(methods, ekfMetrics, "EKF", ...
    [0.5490, 0.3373, 0.2941], [0.8900, 0.4670, 0.7610], ...
    '--', ':', 's', '^');
methods = appendMethod(methods, iekfMetrics, "IEKF", ...
    [0.1216, 0.4667, 0.7059], [1.0000, 0.4980, 0.0549], ...
    '-', '--', 'o', 's');
methods = appendMethod(methods, laplaceMetrics, "Laplace", ...
    [0.0000, 0.6000, 0.6000], [0.8500, 0.6500, 0.1000], ...
    '-', '--', 'h', '>');
end

function label = metricLabel(metrics, fallback)
label = fallback;
if ~isempty(metrics) && isfield(metrics, 'estimatorName') && ...
        strlength(string(metrics.estimatorName)) > 0
    label = string(metrics.estimatorName);
end
end

function methods = appendMethod( ...
    methods, data, label, edgeColor, triangleColor, ...
    edgeStyle, triangleStyle, edgeMarker, triangleMarker)
if isempty(data)
    return;
end
entry.data = data;
entry.label = label;
entry.edgeColor = edgeColor;
entry.triangleColor = triangleColor;
entry.edgeStyle = edgeStyle;
entry.triangleStyle = triangleStyle;
entry.edgeMarker = edgeMarker;
entry.triangleMarker = triangleMarker;
methods(end + 1) = entry; %#ok<AGROW>
end

function validateMetricSequence(candidate, reference, nTime, label)
assert(isequal(candidate.time, reference.time) && ...
    numel(candidate.edge.supportF1) == nTime && ...
    numel(candidate.triangle.supportF1) == nTime && ...
    numel(candidate.edge.relativeError) == nTime && ...
    numel(candidate.triangle.relativeError) == nTime, ...
    'dynsc:PlotMetricSizeMismatch', ...
    '%s metrics must use the common time indices.', label);
if ~isempty(candidate.priorResidualNorm)
    assert(numel(candidate.priorResidualNorm) == nTime, ...
        'dynsc:PlotMetricSizeMismatch', ...
        '%s prior residuals must use the common time indices.', label);
end
referenceHasNormalizedSHD = ...
    isfield(reference.edge, 'normalizedSHD') && ...
    isfield(reference.triangle, 'normalizedSHD');
candidateHasNormalizedSHD = ...
    isfield(candidate.edge, 'normalizedSHD') && ...
    isfield(candidate.triangle, 'normalizedSHD');
assert(candidateHasNormalizedSHD == referenceHasNormalizedSHD, ...
    'dynsc:PlotMetricSizeMismatch', ...
    '%s normalized-SHD availability must match the reference.', label);
if candidateHasNormalizedSHD
    assert(numel(candidate.edge.normalizedSHD) == nTime && ...
        numel(candidate.triangle.normalizedSHD) == nTime, ...
        'dynsc:PlotMetricSizeMismatch', ...
        '%s normalized-SHD metrics must use the common time indices.', ...
        label);
end
end

function plotEstimatorPair(axisHandle, iterations, method, metricName)
plot(axisHandle, iterations, method.data.edge.(metricName), ...
    'Color', method.edgeColor, 'LineStyle', method.edgeStyle, ...
    'Marker', method.edgeMarker, 'MarkerFaceColor', 'none', ...
    'LineWidth', 1.8, 'MarkerSize', 6, ...
    'DisplayName', method.label + " edges");
plot(axisHandle, iterations, method.data.triangle.(metricName), ...
    'Color', method.triangleColor, 'LineStyle', method.triangleStyle, ...
    'Marker', method.triangleMarker, 'MarkerFaceColor', 'none', ...
    'LineWidth', 1.8, 'MarkerSize', 6, ...
    'DisplayName', method.label + " triangles");
end

function formatAxis(axisHandle, yLabel)
ylabel(axisHandle, yLabel);
grid(axisHandle, 'on');
box(axisHandle, 'on');
axisHandle.FontSize = 12;
end

function addChangeMarkers(axisHandle, changeTimes)
for iChange = 1:numel(changeTimes)
    xline(axisHandle, changeTimes(iChange), ':', ...
        'Color', [0.45, 0.10, 0.55], ...
        'LineWidth', 1.4, 'HandleVisibility', 'off');
end
end
