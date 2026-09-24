function cfg = primary_school_config(projectRoot)
%PRIMARY_SCHOOL_CONFIG Configuration for empirical SC construction only.
%
% This configuration is intentionally separate from default_config. It does
% not invoke or modify any filtering, prediction, update, or metric code.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
end

dataFolder = fullfile(projectRoot, "src", "+dynsc", "+data");

cfg = struct();
cfg.experiment.name = "primary_school_sc_construction";

% AHORN/ScHoLP timestamped maximal-simplex stream.
cfg.data.ahornFile = fullfile( ...
    dataFolder, "contact-primary-school.txt");

% Official SocioPatterns source files. The contact file is used only to
% reconstruct and validate the AHORN 1,...,242 node renumbering; topology is
% constructed from the AHORN maximal simplices above.
cfg.data.officialContactsFile = fullfile( ...
    dataFolder, "primaryschool.csv.gz");
cfg.data.metadataFile = fullfile( ...
    dataFolder, "primaryschool_metadata.txt");

% Final cohort: ten students from each of the two classes in one grade.
% Change only these labels to study another same-grade class pair.
cfg.selection.classLabels = ["1A", "1B"];
cfg.selection.nodesPerClass = 10;
cfg.selection.requirePresenceEveryDay = true;
cfg.selection.seed = 17;

% Edges are unions of all contacts within a non-overlapping 10-minute
% window. A triangle is present only if its three nodes co-occurred in one
% raw maximal simplex at a common 20-second timestamp within that window.
cfg.window.durationMinutes = 10;
cfg.window.rawResolutionSeconds = 20;

% These checks verify the external data conversion and the constructed SCs.
cfg.validation.compareRawPairTimeEvents = true;

% Visual diagnostics only; these do not feed back into node selection.
cfg.plot.visible = "off";
cfg.plot.snapshotDay = 1;
cfg.plot.snapshotClockTimes = ["09:00", "10:30", "12:30", "15:30"];

cfg.output.root = fullfile( ...
    projectRoot, "results", cfg.experiment.name);
end
