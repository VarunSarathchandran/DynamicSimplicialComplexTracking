function model = predictionModel(cfg)
%PREDICTIONMODEL Resolve the configured prediction model with legacy default.

if isfield(cfg, 'model') && ~isempty(cfg.model)
    model = string(cfg.model);
else
    model = "linear";
end

assert(isscalar(model) && any(model == ["linear", "closure_aware"]), ...
    'dynsc:InvalidPredictionModel', ...
    'cfg.prediction.model must be "linear" or "closure_aware".');
end
