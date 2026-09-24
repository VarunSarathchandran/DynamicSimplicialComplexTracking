function measurement = buildModel(topology, signals, cfg)
%BUILDMODEL Construct the three stacked virtual measurements at every time.

assert(dynsc.measurement.isConfigured(cfg), ...
    'dynsc:MeasurementNotConfigured', ...
    'Supply the three measurement-noise variances before building the model.');

noiseVariances = cfg.noiseVariances(:);
validateattributes(noiseVariances, {'numeric'}, ...
    {'real', 'finite', 'nonnegative', 'numel', 3});

nEdges = size(topology.B1_full, 2);
nTriangles = size(topology.B2_full, 2);
nState = nEdges + nTriangles;
nTime = numel(signals.time);

assert(size(topology.s1, 2) == nTime && ...
    size(topology.s2, 2) == nTime, ...
    'dynsc:MeasurementTopologyTimeMismatch', ...
    ['The topology support and signal sequence must contain the same ' ...
     'number of time steps.']);

A1 = [speye(nEdges), sparse(nEdges, nTriangles)];
A2 = [sparse(nTriangles, nEdges), speye(nTriangles)];
edgeTotalRow = ones(1, nEdges) * A1;
triangleTotalRow = ones(1, nTriangles) * A2;

C1 = resolveTotal(cfg.C1, topology.s1, nTime, 'C1');
C2 = resolveTotal(cfg.C2, topology.s2, nTime, 'C2');

h = zeros(nState, nTime);
h1 = zeros(nEdges, nTime);
h2 = zeros(nTriangles, nTime);
H = zeros(3, nState, nTime);
y = [zeros(1, nTime); C1; C2];

for iTime = 1:nTime
    [h(:, iTime), h1(:, iTime), h2(:, iTime)] = ...
        dynsc.measurement.buildSmoothnessVector( ...
            topology, signals.Y_node(:, :, iTime), ...
            signals.Y_edge(:, :, iTime), cfg.smoothnessType);

    H(:, :, iTime) = [
        h(:, iTime)'
        full(edgeTotalRow)
        full(triangleTotalRow)
    ];
end

measurement = struct();
measurement.smoothnessType = string(cfg.smoothnessType);
measurement.C1 = C1;
measurement.C2 = C2;
measurement.A1 = A1;
measurement.A2 = A2;
measurement.h = h;
measurement.h1 = h1;
measurement.h2 = h2;
measurement.H = H;
measurement.y = y;
measurement.noiseMean = zeros(3, 1);
measurement.noiseCovariance = diag(noiseVariances);
end

function value = resolveTotal( ...
    configuredValue, groundTruthState, nTime, name)
if (isstring(configuredValue) || ischar(configuredValue)) && ...
        string(configuredValue) == "ground_truth"
    value = sum(double(groundTruthState), 1);
    return;
end

validateattributes(configuredValue, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, name);

if isscalar(configuredValue)
    value = repmat(double(configuredValue), 1, nTime);
else
    assert(numel(configuredValue) == nTime, ...
        'dynsc:MeasurementTotalTimeMismatch', ...
        '%s must be scalar or contain one value per time step.', name);
    value = reshape(double(configuredValue), 1, []);
end
end
