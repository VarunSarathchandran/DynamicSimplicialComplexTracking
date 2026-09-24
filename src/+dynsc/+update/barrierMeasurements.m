function measurement = barrierMeasurements(state, inclusionSlackMatrix)
%BARRIERMEASUREMENTS Evaluate scalar log-barrier measurements and Jacobians.
%
%   h_inc(s) = sqrt(-sum(log(A_inc*s)))
%   h_box(s) = sqrt(-sum(log(s) + log(1-s)))
%
% Jacobians are returned as columns so their information contributions are
% (1/sigma^2) * j * j'.

state = state(:);
validateattributes(state, {'numeric'}, ...
    {'column', 'real', 'finite'});
assert(size(inclusionSlackMatrix, 2) == numel(state), ...
    'dynsc:BarrierInclusionSizeMismatch', ...
    'The inclusion matrix must have one column per state component.');

inclusionSlack = inclusionSlackMatrix * state;
assert(all(inclusionSlack > 0), ...
    'dynsc:NonPositiveInclusionSlack', ...
    'All edge-minus-triangle inclusion slacks must be strictly positive.');
assert(all(state > 0 & state < 1), ...
    'dynsc:StateOutsideOpenBox', ...
    'Every state component must lie strictly between zero and one.');

inclusionBarrier = -sum(log(inclusionSlack));
boxBarrier = -sum(log(state) + log1p(-state));
assert(inclusionBarrier > 0 && isfinite(inclusionBarrier), ...
    'dynsc:InvalidInclusionBarrier', ...
    'The aggregate inclusion barrier must be finite and positive.');
assert(boxBarrier > 0 && isfinite(boxBarrier), ...
    'dynsc:InvalidBoxBarrier', ...
    'The aggregate box barrier must be finite and positive.');

inclusionValue = sqrt(inclusionBarrier);
boxValue = sqrt(boxBarrier);
inclusionBarrierGradient = -inclusionSlackMatrix' * ...
    (1 ./ inclusionSlack);
boxBarrierGradient = -1 ./ state + 1 ./ (1 - state);

measurement = struct();
measurement.inclusionSlack = inclusionSlack;
measurement.inclusionBarrier = inclusionBarrier;
measurement.boxBarrier = boxBarrier;
measurement.inclusionValue = inclusionValue;
measurement.boxValue = boxValue;
measurement.inclusionJacobian = ...
    inclusionBarrierGradient / (2 * inclusionValue);
measurement.boxJacobian = boxBarrierGradient / (2 * boxValue);
end
