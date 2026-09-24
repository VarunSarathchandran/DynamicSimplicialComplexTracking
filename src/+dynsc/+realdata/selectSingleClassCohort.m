function selection = selectSingleClassCohort(events, nodeMap, selectionCfg)
%SELECTSINGLECLASSCOHORT Select one reproducible within-class cohort.

classLabel = string(selectionCfg.classLabel);
assert(isscalar(classLabel), 'dynsc:InvalidPrimarySchoolClass', ...
    'Exactly one class label is required.');
validateattributes(selectionCfg.numNodes, {'numeric'}, ...
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
candidateRows = find(nodeMap.Class == classLabel & available);
assert(numel(candidateRows) >= selectionCfg.numNodes, ...
    'dynsc:InsufficientClassParticipants', ...
    ['Class %s has only %d eligible participants; %d were ' ...
     'requested.'], classLabel, numel(candidateRows), ...
    selectionCfg.numNodes);

stream = RandStream('mt19937ar', 'Seed', selectionCfg.seed);
order = randperm(stream, numel(candidateRows), selectionCfg.numNodes);
selectedRows = candidateRows(order);
selectedNodes = nodeMap(selectedRows, :);
selectedNodes.LocalId = (1:height(selectedNodes)).';
selectedNodes = movevars(selectedNodes, 'LocalId', 'Before', 'AhornId');

selection = struct();
selection.classLabels = classLabel;
selection.nodes = selectedNodes;
selection.seed = selectionCfg.seed;
selection.numNodes = selectionCfg.numNodes;
selection.requirePresenceEveryDay = ...
    logical(selectionCfg.requirePresenceEveryDay);
selection.presentByDay = presentByDay(selectedRows, :);
end
