function [predictedState, predictedCovariance] = predictState( ...
    previousState, previousCovariance, F, Q)
%PREDICTSTATE Perform one linear-Gaussian prediction step.
%
%   x_{k|k-1} = F x_{k-1|k-1}
%   P_{k|k-1} = F P_{k-1|k-1} F' + Q

predictedState = F * previousState;
predictedCovariance = F * previousCovariance * F' + Q;
predictedCovariance = 0.5 * ...
    (predictedCovariance + predictedCovariance');
end
