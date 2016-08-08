function params = thresholdParamsGenerate(varargin)
% params = thresholdParamsGenerate(varargin)
%
% Generate parameters that tell us how to find thresholds.
%
%   method - basic classifier method
%     'svm'
%   signalSource - What response signals to use
%     'isomerizations'
%     'photocurrents'
%   nIntervals - Simulate 1 or 2 interval task
%   kFold - Number of cross-validation folds
%   PCAComponents - Number of PCA components, 0 for no PCA.

% Parameter struct type
params.type = 'threshold';

% Some reasonable parameters
params.method = 'svm';
params.signalSource = 'isomerizations';
params.nIntervals = 2;
params.kFold = 5;
params.PCAComponents = 60;