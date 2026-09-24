function parameters = resolveParameters(parameters, nEdges, nTriangles)
%RESOLVEPARAMETERS Validate and expand closure-aware transition parameters.
%
% Scalar inputs are used by the current experiments. Accepting a scalar or
% one value per simplex keeps the transition utilities ready for later
% edge- or triangle-specific parameters without changing their interfaces.

parameters.beta1 = expandProbability( ...
    parameters, 'beta1', nEdges, 'edge');
parameters.rho1 = expandProbability( ...
    parameters, 'rho1', nEdges, 'edge');
parameters.beta2 = expandProbability( ...
    parameters, 'beta2', nTriangles, 'triangle');
parameters.rho2 = expandProbability( ...
    parameters, 'rho2', nTriangles, 'triangle');
parameters.phi1 = parameters.rho1 - parameters.beta1;
parameters.phi2 = parameters.rho2 - parameters.beta2;
end

function value = expandProbability(parameters, name, nEntries, level)
assert(isfield(parameters, name), 'dynsc:MissingTransitionParameter', ...
    'Missing transition parameter: %s.', name);

value = parameters.(name);
validateattributes(value, {'numeric'}, ...
    {'vector', 'real', 'finite', '>=', 0, '<=', 1});

if isscalar(value)
    value = repmat(double(value), nEntries, 1);
else
    assert(numel(value) == nEntries, ...
        'dynsc:TransitionParameterSizeMismatch', ...
        '%s must be scalar or contain one value per %s.', name, level);
    value = double(value(:));
end
end
