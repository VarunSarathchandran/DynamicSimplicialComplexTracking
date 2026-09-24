function cfg = default_config(projectRoot)
%DEFAULT_CONFIG Configuration for closure-aware dynamic-SC experiments.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
end

cfg = struct();

cfg.experiment.name = "log_barrier_measurements";
cfg.experiment.seed = 7;
cfg.experiment.numRealizations = 10;

cfg.topology.generator = ...
    @dynsc.topology.generateClosureAwareMarkov;
cfg.topology.numNodes = 15;
cfg.topology.edgeProbability = 0.35;
cfg.topology.triangleFraction = 0.50;
cfg.topology.changeInterval = 10;

% Birth-survival probabilities used by the hard ground-truth generator.
cfg.transition.beta1 = 0.05;
cfg.transition.rho1 = 0.98;
cfg.transition.beta2 = 0.05;
cfg.transition.rho2 = 0.98;


cfg.signals.generator = @dynsc.signals.generateFromTopology;
cfg.signals.numTimeSteps = 30;
cfg.signals.numNodeSignals = 1;
cfg.signals.numEdgeSignals = 1;
cfg.signals.smoothnessType = "low_curl";
cfg.signals.sigma = 0;
cfg.signals.normalizeNoise = true;
cfg.signals.verbose = false;

% Prediction model. F is retained for the linear identity baseline; the
% closure-aware model uses f(s) for its mean and its state-dependent Jacobian
% for covariance propagation.
cfg.state_size = nchoosek(cfg.topology.numNodes, 2) + ...
    nchoosek(cfg.topology.numNodes, 3);
cfg.prediction.model = "closure_aware";
% Birth-survival probabilities assumed by all dynamics-informed methods.
% These match the generator by default; change them to study model mismatch.
cfg.prediction.transition.beta1 = 0.05;
cfg.prediction.transition.rho1 = 0.98;
cfg.prediction.transition.beta2 = 0.05;
cfg.prediction.transition.rho2 = 0.98;
cfg.prediction.F = eye(cfg.state_size);
% The first observed topology is the ER draw itself, so assimilate its first
% measurement before applying the Markov transition, which begins at k = 2.
cfg.prediction.skipFirstPrediction = true;
% Recompute the Gaussian moment-matched covariance of the hard Bernoulli
% transition at every prediction. To recover the previous fixed model, set Q
% to "block_diagonal" and use qEdge/qTriangle below, or supply a numeric Q.
cfg.prediction.Q = "bernoulli_second_moment";
cfg.prediction.processNoise.qEdge = 0.1;
cfg.prediction.processNoise.qTriangle = 0.1;

% Initialize every edge with the Erdos-Renyi probability p and every triangle
% with the approximate marginal probability triangleFraction * p^3. Under
% the independent-edge Erdos-Renyi prior, p^3 is the probability that all
% three boundary edges are present. The corresponding Bernoulli variances
% form the diagonal of P0. The "ground_truth" oracle option replaces only s0
% with the first true topology and retains this same prior covariance. Use
% method = "manual" to supply initialState and initialCovariance directly.
 cfg.prediction.initialization.method = "bernoulli";
 cfg.prediction.initialization.varianceFloor = 1e-6;
 cfg.prediction.initialState = [];
 cfg.prediction.initialCovariance = [];
% cfg.prediction.initialization.method = "ground_truth";
%cfg.prediction.initialization.method = "manual";
%cfg.prediction.initialState = zeros(cfg.state_size,1);
%cfg.prediction.initialCovariance = eye(cfg.state_size);
% Three virtual measurements: smoothness, total edge weight, and total
% triangle weight. C1 and C2 are resolved from the ground-truth realization.
cfg.measurement.builder = @dynsc.measurement.buildModel;
cfg.measurement.smoothnessType = "low_curl";
cfg.measurement.C1 = "ground_truth";
cfg.measurement.C2 = "ground_truth";

% [Var(eta_k); Var(mu_k^1); Var(mu_k^2)]. Supply these before enabling the
% measurement model; no values are assumed here.
cfg.measurement.noiseVariances = [
    0.1
    0.01
    0.01
];

% Correction model. "constrained_map" uses the existing Gurobi mean update
% and ordinary linear-measurement covariance approximation.
% "barrier_virtual_measurement" adds the inclusion and box log barriers to
% the MAP mean objective and uses their final-estimate Jacobians in the
% IEKF/Gauss-Newton covariance update. "barrier_exact_hessian" keeps that
% same MAP objective but uses its exact log-barrier Hessian for covariance.
% "barrier_standard_ekf_projection" performs a one-shot EKF update for the
% two square-root barriers and projects its mean into the SC feasible set.
%cfg.update.method = "constrained_map";
cfg.update.method = "barrier_virtual_measurement";

% If h_inc(s)^2 = -sum(log(A_inc*s)) and
% h_box(s)^2 = -sum(log(s) + log(1-s)), then lambda = 1/(2*sigma^2)
% is the corresponding Gaussian virtual-measurement weight. Thus the
% implied variances are 1/(2*lambda_inc) and 1/(2*lambda_box).
cfg.update.barrier.lambda_inc = 0.01;
cfg.update.barrier.lambda_box = 0.01;
cfg.update.barrier.interiorMargin = 1e-8;
cfg.update.barrier.maxIterations = 500;
cfg.update.barrier.maxFunctionEvaluations = 5000;
cfg.update.barrier.optimalityTolerance = 1e-9;
cfg.update.barrier.stepTolerance = 1e-10;

cfg.metrics.evaluator = @dynsc.metrics.evaluateRecoveryOverTime;
cfg.metrics.supportThreshold = 0.2;

% The measurement-only static problem is disabled by default because it can
% be severely underdetermined. Set enabled = true only when it is a relevant
% benchmark for a particular experiment.
cfg.baseline.static.enabled = false;
cfg.baseline.solver = @dynsc.baseline.runStaticMeasurementOnly;
% Persistence-prior baseline: identity state transition and Jacobian with a
% fixed block-diagonal process covariance. It shares the initialization and
% measurements with the primary filter, but deliberately uses the historical
% constrained-MAP mean and ordinary covariance correction. qEdge and
% qTriangle are inherited from cfg.prediction.processNoise.
cfg.baseline.identityTransition.enabled = true;

% Constrained-RLS baseline: F = I and Q = 0. It uses the same resolved
% initialization and measurements as the other filters, together with the
% constrained-MAP mean and ordinary linear-measurement covariance update.
cfg.baseline.constrainedRls.enabled = true;

% Historical-correction baseline. When the primary filter uses the barrier
% virtual-measurement update, this baseline retains the same initialization,
% prediction f, Jacobian F_k, and process covariance Q_k, but uses the former
% constrained-MAP mean and LambdaNormal^{-1} covariance correction. It is
% skipped automatically when the primary already uses constrained_map.
cfg.baseline.previousUpdate.enabled = true;

% Exact-Hessian covariance ablation. This uses the same closure-aware
% prediction and exactly the same fixed log-barrier MAP objective as the
% barrier_virtual_measurement method. Only the corrected covariance changes:
% the exact barrier Hessians replace the IEKF/Gauss-Newton outer products.
cfg.baseline.exactHessianUpdate.enabled = true;

% Standard-EKF nonlinear-measurement baseline. The square-root log-barrier
% measurements are linearized once at the predicted mean. After
% the ordinary augmented Kalman update, its mean is projected into the box
% and inclusion feasible set in the PPost^{-1} metric. The EKF covariance is
% retained unchanged by that mean-only projection.
cfg.baseline.standardEkfBarrier.enabled = true;

cfg.plot.stateSnapshots.enabled = true;
cfg.plot.stateSnapshots.interval = 2;
% Show mean +/- one standard deviation in plots aggregated over realizations.
% Set this to false to plot only the realization-averaged curves.
cfg.plot.average.showStandardDeviation = false;

cfg.output.root = fullfile(projectRoot, "results");
cfg.output.saveRealizations = true;
end
