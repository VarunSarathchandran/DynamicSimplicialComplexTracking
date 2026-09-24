function projectRoot = startup_dynamic_sc()
%STARTUP_DYNAMIC_SC Add the dynamic code and required generation helpers.

projectRoot = fileparts(mfilename('fullpath'));

folders = [
    string(projectRoot)
    fullfile(projectRoot, "configs")
    fullfile(projectRoot, "experiments")
    fullfile(projectRoot, "src")
    fullfile(projectRoot, "tests")
    fullfile(projectRoot, "utils", "Generation")
];

for iFolder = 1:numel(folders)
    if isfolder(folders(iFolder))
        addpath(folders(iFolder));
    end
end

addGurobiPathIfAvailable();

fprintf('Dynamic-SC project ready: %s\n', projectRoot);
fprintf('Added utils/Generation; other legacy folders were not added.\n');
end

function addGurobiPathIfAvailable()
if exist('gurobi', 'file') == 2
    return;
end

gurobiHome = getenv('GUROBI_HOME');
if ~isempty(gurobiHome)
    matlabFolder = fullfile(gurobiHome, 'matlab');
    if isfile(fullfile(matlabFolder, 'gurobi.m'))
        addpath(matlabFolder);
        return;
    end
end

installations = dir( ...
    '/Library/gurobi*/macos_universal2/matlab/gurobi.m');
if ~isempty(installations)
    addpath(installations(end).folder);
end
end
