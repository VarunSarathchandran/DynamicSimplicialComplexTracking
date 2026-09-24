function cfg = primary_school_tracking_config(projectRoot)
%PRIMARY_SCHOOL_TRACKING_CONFIG Semi-synthetic empirical-topology experiment.
%
% Day 1 is used only for calibration. All calibrated quantities and tuned
% hyperparameters are frozen before the same clock interval on day 2 is
% processed. The 10-minute topology construction itself is unchanged.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
end

dataCfg = primary_school_config(projectRoot);
methodCfg = default_config(projectRoot);

cfg = dataCfg;
% The raw pair-event cross-check is a data-conversion audit, not part of
% the experiment; the reported run had it disabled.
cfg.validation.compareRawPairTimeEvents = false;
cfg.experiment.name = "primary_school_tracking_1B_full25_day1counts";
cfg.experiment.seed = 2718;
cfg.experiment.numRealizations = 1;

% The cohort is the whole of class 1B: all 25 of its members are present on
% both days, so numNodes equals the eligible class size and no within-class
% subsampling takes place. The selection seed therefore only permutes local
% labels and cannot influence which topology is produced. Class 1B is the
% only class whose full-cohort day-1 and day-2 segments contain at least one
% triangle in every window; that was established from the constructed
% topology alone, with no recovery result computed beforehand.
cfg.selection.mode = "single_class";
cfg.selection.classLabel = "1B";
cfg.selection.numNodes = 25;
cfg.selection.seed = 8;

% This is an experimental segment, not a new aggregation window: it contains
% 12 non-overlapping 20-minute SC states from 09:20 through 13:20. The
% 20-minute width is the narrowest that leaves no empty-triangle window on
% either day for this cohort.
cfg.window.durationMinutes = 20;
cfg.segment.startClock = "09:20";
cfg.segment.durationHours = 4;
cfg.segment.calibrationDay = 1;
cfg.segment.testDay = 2;

% Reuse exactly the signal-generation model used by simulated experiments.
% The feature counts are raised well above the simulated defaults because a
% single joint smoothness measurement built from only three features cannot
% rank candidate edges reliably.
cfg.signals = methodCfg.signals;
cfg.signals.numNodeSignals = 50;
cfg.signals.numEdgeSignals = 50;
cfg.signals.sigma = 0;
cfg.signals.verbose = false;

% Estimate beta/rho as unsmoothed maximum-likelihood ratios, successes over
% trials, by setting both pseudo-count parameters to zero. Every transition
% category in this segment supplies at least one trial, so no smoothing is
% required; a category with no trials is rejected rather than smoothed.
cfg.calibration.betaPriorAlpha = 0;
cfg.calibration.betaPriorBeta = 0;
cfg.calibration.varianceFloor = 1e-6;

% Hyperparameter handling:
%   "fixed" uses the absolute values below without a day-1 grid search.
%   "calibrated_screen" retains the original day-1 multiplier screening.
% Transition probabilities and initialization are still estimated on day 1.
% lambdaInc is set near the value at which the inclusion-barrier information
% matches the linear-measurement information, mean diag(Ainc'*diag(1./d.^2)*
% Ainc) against mean diag(H'*R^-1*H), so that inclusion actually reaches the
% posterior covariance. lambdaBox is set far below it because the box Hessian
% carries 1/x^2 + 1/(1-x)^2 and most candidate simplices are estimated near
% zero; at equal weights the box term exceeded the measurement information by
% five to six orders of magnitude and dominated every covariance.
cfg.tuning.mode = "fixed";
cfg.tuning.fixed.smoothnessVariance = 0.001;
cfg.tuning.fixed.lambdaInc = 0.03;
cfg.tuning.fixed.lambdaBox = 1e-12;

% These fields are used only when mode is "calibrated_screen". The central
% smoothness variance is estimated from day-1 residuals and multiplied by
% the candidates below. Set enabled=false to use only the first profile.
cfg.tuning.enabled = false;
cfg.tuning.estimator = "laplace";
cfg.tuning.smoothnessVarianceMultiplier = [1, 1/3, 3, 1, 1];
cfg.tuning.lambdaInc = [0.01, 0.01, 0.01, 0.003, 0.03];
cfg.tuning.lambdaBox = [0.01, 0.01, 0.01, 0.003, 0.03];

cfg.prediction = methodCfg.prediction;
cfg.prediction.initialization.method = "manual";
cfg.prediction.initialState = [];
cfg.prediction.initialCovariance = [];
cfg.prediction.skipFirstPrediction = true;
cfg.prediction.Q = "bernoulli_second_moment";

cfg.measurement = methodCfg.measurement;
cfg.measurement.C1 = [];
cfg.measurement.C2 = [];
cfg.measurement.noiseVariances = [];
% "calibrated" uses fixed day-1 means and empirical variances on day 2.
% "ground_truth" is an explicitly oracle sanity check using the exact
% time-varying counts and the small variances below.
% "calibration_matched" takes the time-matched counts observed on day 1 and
% freezes them, exactly as the transition probabilities are frozen, so no
% test-day quantity enters the measurement. The variances below are fixed
% pseudo-measurement weights rather than calibrated noise variances: they
% make the counts near-hard constraints, which the triangle block requires.
cfg.measurement.countSource = "calibration_matched";
cfg.measurement.matchedCountVariances = [0.001; 0.001];
cfg.measurement.oracleCountVariances = [0.001; 0.001];

cfg.update = methodCfg.update;
cfg.metrics = methodCfg.metrics;
cfg.metrics.supportThreshold = 0.2;

cfg.methods.constrainedRls = true;
cfg.methods.constrainedMap = true;
cfg.methods.ekf = true;
cfg.methods.iekf = true;
cfg.methods.laplace = true;

% Optionally save test-day state snapshots for one enabled method. Valid
% method names are constrainedRls, constrainedMap, ekf, iekf, and laplace.
cfg.testSnapshots.enabled = true;
cfg.testSnapshots.method = "laplace";
cfg.testSnapshots.interval = 5;

cfg.plot.visible = "off";
cfg.output.root = fullfile(projectRoot, "results", cfg.experiment.name);
cfg.output.saveDetailedFilters = false;
end
