function params = LMPlaneInstanceParamsGenerate(varargin)
% params = LMPlaneInstanceParamsGenerate(varargin)
%
% Generate parameters to describe color directions, contrasts, etc for
% studying performance in the LM plane.
%
%   trialsNum - Number of noisy response instances of each type to generate
%   deltaAngle - Sampling angle in LM contrast plane
%   baseStimulusLength - Vector length of base stimuli in each color dir
%   scaleIntoMonitorGamut - Scale base contrast into monitor gamut?
%                           If this is true, baseStimulusLength has no
%                           effect.
%   nContrastsPerDirection - Number of contrasts per color direction.
%   lowContrast - Low constrast for contrast sampling
%   highContrast - High contrast for contrast sampling
%   contrastScale - 'log' or 'linear', how contrasts are sampled

% Parameter struct type
params.type = 'LMPlaneInstance';

% Define how many noisy data instances to generate
params.trialsNum = 100;

% Delta angle sampling in LM plane (samples between 0 and 180 degrees)
params.deltaAngle = 45;
params.baseStimulusLength = 1;
params.scaleIntoMonitorGamut = true;

% Number of contrasts to run in each color direction
params.nContrastsPerDirection = 5; 
params.lowContrast = 0.001;
params.highContrast = 0.2;
params.contrastScale = 'log';    % choose between 'linear' and 'log'
