function dataset = buildPrimarySchoolDataset(cfg)
%BUILDPRIMARYSCHOOLDATASET Construct fixed-cohort windowed empirical SCs.

events = dynsc.realdata.loadAhornPrimarySchool(cfg.data.ahornFile);
assert(cfg.window.rawResolutionSeconds == events.rawResolutionSeconds, ...
    'dynsc:PrimarySchoolResolutionMismatch', ...
    'The configured and observed raw temporal resolutions disagree.');
[nodeMap, rawContacts] = dynsc.realdata.buildPrimarySchoolNodeMap( ...
    cfg.data.officialContactsFile, cfg.data.metadataFile);

assert(max(events.nodeIds) == height(nodeMap) && ...
    numel(events.nodeIds) == height(nodeMap), ...
    'dynsc:PrimarySchoolNodeMapMismatch', ...
    'AHORN node identifiers and reconstructed metadata mapping disagree.');

if cfg.validation.compareRawPairTimeEvents
    dynsc.realdata.validateAhornPairTimeEvents( ...
        events, nodeMap, rawContacts);
end

selection = dynsc.realdata.selectTwoClassCohort( ...
    events, nodeMap, cfg.selection);
eventDay = dateshift(events.time, 'start', 'day');
uniqueDays = unique(eventDay, 'stable');
nDays = numel(uniqueDays);

days = repmat(struct( ...
    'index', [], 'sourceDate', NaT, 'topology', []), 1, nDays);
diagnosticTables = cell(nDays, 1);

for iDay = 1:nDays
    keep = eventDay == uniqueDays(iDay);
    [topology, diagnostics] = ...
        dynsc.realdata.buildWindowedOrderTwoSC( ...
            events.time(keep), events.simplices(keep), selection, ...
            cfg.window.durationMinutes, iDay);
    days(iDay).index = iDay;
    days(iDay).sourceDate = uniqueDays(iDay);
    days(iDay).topology = topology;
    diagnosticTables{iDay} = diagnostics;
end

dataset = struct();
dataset.name = "contact-primary-school";
dataset.description = [ ...
    "Real temporal contact topology with synthetic signals to be " ...
    "generated in a later experiment stage."];
dataset.sourceRepresentation = "AHORN timestamped maximal simplices";
dataset.rawResolutionSeconds = events.rawResolutionSeconds;
dataset.windowMinutes = cfg.window.durationMinutes;
dataset.nodeMap = nodeMap;
dataset.selection = selection;
dataset.days = days;
dataset.windowDiagnostics = vertcat(diagnosticTables{:});
dataset.construction = [ ...
    "Edges are pair contacts observed anywhere in each coarse window. " ...
    "Triangles are triples co-occurring in a common raw 20-second " ...
    "maximal simplex; aggregate-graph 3-cliques are not filled."];
end
