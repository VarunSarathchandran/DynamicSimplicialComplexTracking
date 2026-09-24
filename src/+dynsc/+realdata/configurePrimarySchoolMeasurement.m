function cfg = configurePrimarySchoolMeasurement( ...
    cfg, calibration, smoothnessVariance)
%CONFIGUREPRIMARYSCHOOLMEASUREMENT Resolve calibrated or oracle counts.

validateattributes(smoothnessVariance, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
countSource = string(cfg.countSource);
assert(isscalar(countSource), 'dynsc:InvalidCountSource', ...
    'The measurement count source must be scalar text.');

switch countSource
    case "calibrated"
        cfg.C1 = calibration.C1;
        cfg.C2 = calibration.C2;
        countVariances = calibration.countVariances;
    case "ground_truth"
        cfg.C1 = "ground_truth";
        cfg.C2 = "ground_truth";
        countVariances = cfg.oracleCountVariances(:);
        validateattributes(countVariances, {'numeric'}, ...
            {'real', 'finite', 'positive', 'numel', 2});
    case "calibration_matched"
        % Time-matched counts observed on the calibration day, frozen and
        % applied to the test segment covering the same clock interval.
        % This uses no test-day quantity and is therefore not an oracle.
        cfg.C1 = calibration.C1Sequence;
        cfg.C2 = calibration.C2Sequence;
        countVariances = cfg.matchedCountVariances(:);
        validateattributes(countVariances, {'numeric'}, ...
            {'real', 'finite', 'positive', 'numel', 2});
    otherwise
        error('dynsc:InvalidCountSource', ...
            ['countSource must be "calibrated", "calibration_matched" ' ...
             'or "ground_truth".']);
end

cfg.countSource = countSource;
cfg.noiseVariances = [smoothnessVariance; countVariances];
end
