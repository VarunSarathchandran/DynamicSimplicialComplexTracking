function selection = selectTwoClassCohort( ...
    events, nodeMap, selectionCfg)
%SELECTTWOCLASSCOHORT Select a fixed reproducible two-class cohort.

classLabels = string(selectionCfg.classLabels(:));
assert(numel(classLabels) == 2 && classLabels(1) ~= classLabels(2), ...
    'dynsc:InvalidPrimarySchoolClasses', ...
    'Exactly two distinct class labels are required.');
validateattributes(selectionCfg.nodesPerClass, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(selectionCfg.seed, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'});

eventDays = dateshift(events.time, 'start', 'day');
uniqueDays = unique(eventDays, 'stable');
presentByDay = false(height(nodeMap), numel(uniqueDays));
for iDay = 1:numel(uniqueDays)
    dayEvents = eventDays == uniqueDays(iDay);
    nodesThisDay = unique([events.simplices{dayEvents}]);
    presentByDay(nodesThisDay, iDay) = true;
end

if selectionCfg.requirePresenceEveryDay
    available = all(presentByDay, 2);
else
    available = any(presentByDay, 2);
end

stream = RandStream('mt19937ar', 'Seed', selectionCfg.seed);
selectedRows = zeros(2 * selectionCfg.nodesPerClass, 1);
cursor = 0;
for iClass = 1:2
    candidateRows = find( ...
        nodeMap.Class == classLabels(iClass) & available);
    assert(numel(candidateRows) >= selectionCfg.nodesPerClass, ...
        'dynsc:InsufficientClassParticipants', ...
        ['Class %s has only %d eligible participants; %d were ' ...
         'requested.'], classLabels(iClass), numel(candidateRows), ...
        selectionCfg.nodesPerClass);
    order = randperm(stream, numel(candidateRows), ...
        selectionCfg.nodesPerClass);
    destination = cursor + (1:selectionCfg.nodesPerClass);
    selectedRows(destination) = candidateRows(order);
    cursor = cursor + selectionCfg.nodesPerClass;
end

selectedNodes = nodeMap(selectedRows, :);
selectedNodes.LocalId = (1:height(selectedNodes)).';
selectedNodes = movevars(selectedNodes, 'LocalId', 'Before', 'AhornId');

selection = struct();
selection.classLabels = classLabels.';
selection.nodes = selectedNodes;
selection.seed = selectionCfg.seed;
selection.nodesPerClass = selectionCfg.nodesPerClass;
selection.requirePresenceEveryDay = ...
    logical(selectionCfg.requirePresenceEveryDay);
selection.presentByDay = presentByDay(selectedRows, :);
end
