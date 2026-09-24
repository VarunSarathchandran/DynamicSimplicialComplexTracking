function summary = run_initial_experiment(cfg)
%RUN_INITIAL_EXPERIMENT Run the minimal synthetic dynamic-SC experiment.
%
%   summary = run_initial_experiment
%   summary = run_initial_experiment(cfg)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
startup_dynamic_sc();

if nargin < 1 || isempty(cfg)
    cfg = default_config(projectRoot);
end

summary = dynsc.runExperiment(cfg);
end
