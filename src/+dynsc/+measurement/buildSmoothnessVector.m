function [h, h1, h2] = buildSmoothnessVector( ...
    topology, nodeSignals, edgeSignals, smoothnessType)
%BUILDSMOOTHNESSVECTOR Reproduce the h-vector used by the MILP solvers.

fullB1 = topology.B1_full;
fullB2 = topology.B2_full;

h1 = diag(fullB1' * nodeSignals * nodeSignals' * fullB1);
h1 = h1(:);

switch string(smoothnessType)
    case "low_curl"
        Z = edgeSignals' * fullB2;
        h2 = sum(Z.^2, 1)';

    case {"new", "new_smoothness"}
        C = abs(fullB2);
        Z = edgeSignals;
        u = sum(Z.^2, 2);
        q = C' * u;
        M = Z' * C;
        m2 = sum(M.^2, 1)';
        h2 = 3 * q - m2;

    otherwise
        error('dynsc:UnknownMeasurementSmoothnessType', ...
            'Unknown measurement smoothness type: %s', smoothnessType);
end

h1 = h1 ./ size(nodeSignals, 2);
h2 = h2 ./ size(edgeSignals, 2);
h  = [h1; h2];
end
