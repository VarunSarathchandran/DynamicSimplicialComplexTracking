function validateAhornPairTimeEvents(events, nodeMap, rawContacts)
%VALIDATEAHORNPAIRTIMEEVENTS Compare AHORN faces with official contacts.
%
% Expanding every maximal simplex into pair faces must reproduce exactly the
% official set of pair-time contact observations. Duplicates caused by an
% edge belonging to multiple maximal cliques are removed before comparison.

nNodes = height(nodeMap);
maxOriginalId = max(nodeMap.OriginalId);
originalToAhorn = zeros(maxOriginalId, 1);
originalToAhorn(nodeMap.OriginalId) = nodeMap.AhornId;

rawId1 = originalToAhorn(rawContacts.OriginalId1);
rawId2 = originalToAhorn(rawContacts.OriginalId2);
rawLow = min(rawId1, rawId2);
rawHigh = max(rawId1, rawId2);
rawPairCode = sub2ind([nNodes, nNodes], rawLow, rawHigh);
rawKeys = unique([rawContacts.TimeSeconds, rawPairCode], 'rows');

numPairFaces = sum(arrayfun( ...
    @(n) nchoosek(n, 2), events.numVertices));
eventSeconds = zeros(numPairFaces, 1);
eventPairCode = zeros(numPairFaces, 1);
cursor = 0;
secondsFromEpoch = posixtime(events.time);

for iEvent = 1:events.numEvents
    pairs = nchoosek(events.simplices{iEvent}, 2);
    nPairs = size(pairs, 1);
    destination = cursor + (1:nPairs);
    eventSeconds(destination) = secondsFromEpoch(iEvent);
    eventPairCode(destination) = sub2ind( ...
        [nNodes, nNodes], pairs(:, 1), pairs(:, 2));
    cursor = cursor + nPairs;
end

eventKeys = unique([eventSeconds, eventPairCode], 'rows');
assert(isequal(rawKeys, eventKeys), 'dynsc:AhornConversionMismatch', ...
    ['Expanded AHORN maximal simplices do not reproduce the official ' ...
     'pair-time contact data under the reconstructed node mapping.']);
end
