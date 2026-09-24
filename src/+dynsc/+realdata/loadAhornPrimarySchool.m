function events = loadAhornPrimarySchool(filename)
%LOADAHORNPRIMARYSCHOOL Read AHORN timestamped maximal simplices.

assert(isfile(filename), 'dynsc:MissingPrimarySchoolData', ...
    'AHORN data file not found: %s', filename);

lines = strip(readlines(filename));
lines = lines(strlength(lines) > 0);
assert(numel(lines) >= 2, 'dynsc:EmptyPrimarySchoolData', ...
    'AHORN data file contains no simplex records.');

header = jsondecode(char(lines(1)));
payload = lines(2:end);
nodeText = extractBefore(payload, ' {"time": ');
timeText = extractBetween(payload, '"time": "', '"}');

assert(all(strlength(nodeText) > 0) && all(strlength(timeText) > 0), ...
    'dynsc:InvalidAhornFormat', ...
    'Could not parse one or more AHORN simplex records.');

nEvents = numel(payload);
simplices = cell(nEvents, 1);
numVertices = zeros(nEvents, 1);
for iEvent = 1:nEvents
    simplex = sscanf(strrep(char(nodeText(iEvent)), ',', ' '), '%d').';
    assert(~isempty(simplex) && numel(unique(simplex)) == numel(simplex), ...
        'dynsc:InvalidAhornSimplex', ...
        'AHORN event %d contains no nodes or repeated nodes.', iEvent);
    simplices{iEvent} = simplex;
    numVertices(iEvent) = numel(simplex);
end

eventTime = datetime(timeText, ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ssXXX', 'TimeZone', 'UTC');

events = struct();
headerFields = fieldnames(header);
assert(isscalar(headerFields), 'dynsc:InvalidAhornHeader', ...
    'Expected one format-version field in the AHORN header.');
events.formatVersion = string(header.(headerFields{1}));
events.time = eventTime;
events.simplices = simplices;
events.numVertices = numVertices;
events.numEvents = nEvents;
events.nodeIds = unique([simplices{:}]);
events.rawResolutionSeconds = 20;

assert(issorted(eventTime), 'dynsc:UnsortedAhornEvents', ...
    'AHORN simplex events must be sorted by timestamp.');
end
