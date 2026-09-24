function [figures, snapshotIndices] = plotSingleStateSnapshots( ...
    filterResult, topology, methodName, visible, snapshotInterval)
%PLOTSINGLESTATESNAPSHOTS Plot one estimate against time-matched truth.
%
% The visual style mirrors plotUpdatedStateSnapshots without changing the
% existing simulated-experiment plotting function.

if nargin < 3 || strlength(string(methodName)) == 0
    methodName = "Estimate";
end
if nargin < 4
    visible = true;
end
if nargin < 5
    snapshotInterval = 5;
end
validateattributes(snapshotInterval, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});

visibility = "off";
if visible
    visibility = "on";
end

state = filterResult.updatedState;
nEdges = size(topology.s1, 1);
nTriangles = size(topology.s2, 1);
nTime = size(state, 2);
assert(isequal(size(state), [nEdges + nTriangles, nTime]) && ...
    size(topology.s1, 2) == nTime && ...
    size(topology.s2, 2) == nTime, ...
    'dynsc:StateSnapshotSizeMismatch', ...
    'The estimated state and ground-truth topology must have matching sizes.');

snapshotIndices = unique( ...
    [1, snapshotInterval:snapshotInterval:nTime, nTime], 'stable');
figures = gobjects(1, numel(snapshotIndices));
edgeColor = [0.1216, 0.4667, 0.7059];
triangleColor = [1.0000, 0.4980, 0.0549];
triangleRows = nEdges + (1:nTriangles);

for iSnapshot = 1:numel(snapshotIndices)
    k = snapshotIndices(iSnapshot);
    figures(iSnapshot) = figure( ...
        'Color', 'w', 'Visible', visibility, ...
        'Position', [100, 100, 1200, 650]);
    layout = tiledlayout(figures(iSnapshot), 2, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    plotBlock(nexttile(layout), state(1:nEdges, k), ...
        topology.s1(:, k), edgeColor, '-', 'd', methodName, ...
        'Candidate edge index', 'Updated edge weight', 'Edge state');
    plotBlock(nexttile(layout), state(triangleRows, k), ...
        topology.s2(:, k), triangleColor, '--', 's', methodName, ...
        'Candidate triangle index', 'Updated triangle weight', ...
        'Triangle state');

    if isfield(topology, 'windowStart') && ...
            numel(topology.windowStart) >= k
        clockLabel = char(string(topology.windowStart(k), 'HH:mm'));
        title(layout, sprintf( ...
            '%s state and ground truth at k = %d (%s)', ...
            methodName, k, clockLabel));
    else
        title(layout, sprintf( ...
            '%s state and ground truth at k = %d', methodName, k));
    end
end
end

function plotBlock(axisHandle, estimate, truth, color, lineStyle, marker, ...
    methodName, xLabel, yLabel, axisTitle)
nEntries = numel(estimate);
hold(axisHandle, 'on');
stem(axisHandle, 1:nEntries, estimate, ...
    'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.6, ...
    'Marker', marker, 'MarkerSize', 6, 'MarkerFaceColor', 'none', ...
    'ShowBaseLine', 'off', 'DisplayName', methodName);

trueSupport = find(truth > 0);
if isempty(trueSupport)
    supportX = nan;
    supportY = nan;
else
    supportX = trueSupport;
    supportY = ones(size(trueSupport));
end
plot(axisHandle, supportX, supportY, 'ko', 'LineStyle', 'none', ...
    'MarkerFaceColor', 'none', 'MarkerSize', 9, 'LineWidth', 1.5, ...
    'DisplayName', 'Ground-truth support');

xlabel(axisHandle, xLabel);
ylabel(axisHandle, yLabel);
title(axisHandle, axisTitle);
xlim(axisHandle, [0.5, max(2, nEntries + 0.5)]);
ylim(axisHandle, [0, 1.08]);
grid(axisHandle, 'on');
box(axisHandle, 'on');
axisHandle.FontSize = 12;
legend(axisHandle, 'Location', 'northeastoutside');
end
