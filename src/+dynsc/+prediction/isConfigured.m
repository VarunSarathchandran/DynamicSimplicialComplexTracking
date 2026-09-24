function tf = isConfigured(cfg)
%ISCONFIGURED True when the transition model and initialization are set.

model = dynsc.prediction.predictionModel(cfg);
if isempty(cfg.Q) && isempty(cfg.F)
    tf = false;
    return;
end

assert(~isempty(cfg.Q), 'dynsc:IncompletePredictionConfiguration', ...
    'Prediction requires Q.');
if model == "linear"
    assert(~isempty(cfg.F), 'dynsc:IncompletePredictionConfiguration', ...
        'Linear prediction requires F and Q.');
end

method = initializationMethod(cfg);
if method == "manual"
    manualProvided = [
        ~isempty(cfg.initialState)
        ~isempty(cfg.initialCovariance)
    ];
    assert(all(manualProvided), ...
        'dynsc:IncompletePredictionConfiguration', ...
        ['Manual initialization requires both initialState and ' ...
         'initialCovariance.']);
elseif ~any(method == ["bernoulli", "ground_truth"])
    error('dynsc:InvalidInitializationMethod', ...
        ['cfg.prediction.initialization.method must be "bernoulli" ' ...
         '"ground_truth", or "manual".']);
end

tf = true;
end

function method = initializationMethod(cfg)
% Configurations created before initialization modes remain manual.
if isfield(cfg, 'initialization') && ...
        isfield(cfg.initialization, 'method')
    method = string(cfg.initialization.method);
else
    method = "manual";
end

assert(isscalar(method), 'dynsc:InvalidInitializationMethod', ...
    'The prediction initialization method must be scalar text.');
end
