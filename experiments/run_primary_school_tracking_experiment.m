function result = run_primary_school_tracking_experiment(cfg)
%RUN_PRIMARY_SCHOOL_TRACKING_EXPERIMENT Calibrate on day 1 and test on day 2.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
startup_dynamic_sc();

if nargin < 1 || isempty(cfg)
    cfg = primary_school_tracking_config(projectRoot);
end

result = dynsc.realdata.runPrimarySchoolTrackingExperiment(cfg);
end
