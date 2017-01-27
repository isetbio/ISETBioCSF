function params = thresholdParamsGenerate(varargin)
% params = thresholdParamsGenerate(varargin)
%
% Generate parameters that tell us how to find thresholds.
%
%   method - basic classifier method
%     'svm' - Support vector machine
%     'mlpt' - Maximum likelihood, signal/noise known exactly
%     'mlpe' - Maximum likelihood, signal/noise learned from examples
%   signalSource - What response signals to use
%     'isomerizations'
%     'photocurrents'
%   nIntervals - Simulate 1 or 2 interval task
%   kFold - Number of cross-validation folds
%   PCAComponents - Number of PCA components, 0 for no PCA.
%   criterionFraction - Threshold is at this fraction correct.
%   evidenceIntegrationTime - How many milliseconds of the response to use
%   in the classified. An empty matrix results in using all of it.
%
% Parameter struct type
params.type = 'threshold';

% Some reasonable parameters
params.method = 'svm';
params.signalSource = 'isomerizations';
params.V1filterBank = [];
params.nIntervals = 2;
params.kFold = 5;
params.STANDARDIZE = true;
params.PCAComponents = 60;
params.criterionFraction = 0.75;
params.evidenceIntegrationTime = [];
end
