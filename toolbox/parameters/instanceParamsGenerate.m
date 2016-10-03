function params = instanceParamsGenerate(varargin)
% params = instanceParamsGenerate(varargin)
%
% Generate parameters to describe color directions, contrasts, etc for
% studying performance.
%
% Key/value pairs
%   'instanceType' - String (default 'LMPlane') Type of spatial pattern
%     'LMPlane' - Gabor pattern
%     'contrasts' - Spatial spot on background
%
% If instanceType is 'LMPlane', these are the parameters:
%   trialsNum - Number of noisy response instances of each type to generate
%   deltaAngle - Sampling angle in LM contrast plane
%   baseStimulusLength - Vector length of base stimuli in each color dir
%   scaleIntoMonitorGamut - Scale base contrast into monitor gamut?
%     If this is true, baseStimulusLength has no effect.
%   nContrastsPerDirection - Number of contrasts per color direction.
%   lowContrast - Low constrast for contrast sampling
%   highContrast - High contrast for contrast sampling
%   contrastScale - 'log' or 'linear', how contrasts are sampled
%
% If instance type is 'contrasts', these are the parameters:
%   trialsNum - Number of noisy response instances of each type to generate
%   nContrastsPerDirection - Number of contrasts per color direction.
%   lowContrast - Low constrast for contrast sampling
%   highContrast - High contrast for contrast sampling
%   contrastScale - 'log' or 'linear', how contrasts are sampled
%
% See also
%   testConeContastsFromTestDirectionParams

% Parse input
p = inputParser; p.KeepUnmatched = true;
p.addParameter('instanceType','LMPlane',@ischar);
p.parse(varargin{:});
            
% Parameter struct type
params.type = 'Instance';
params.instanceType = p.Results.instanceType;

switch (params.instanceType)
    case 'LMPlane'
        % Define how many noisy data instances to generate
        params.trialsNum = 100;
        
        % Delta angle sampling in LM plane (samples between 0 and 180 degrees)
        params.startAngle = 0;
        params.deltaAngle = 90;
        params.nAngles = 2;
        params.baseStimulusLength = 1;
        params.scaleIntoMonitorGamut = true;
        
        % Number of contrasts to run in each color direction
        % Choose between 'linear' and 'log'
        params.nContrastsPerDirection = 5;
        params.lowContrast = 0.001;
        params.highContrast = 0.15;
        params.contrastScale = 'log';    
    case 'contrasts'
        % Define how many noisy data instances to generate
        params.trialsNum = 100;
        
        % This must be false, but is here to allow code to run
        params.scaleIntoMonitorGamut = false;

        % Number of contrasts to run in each color direction
        % Choose between 'linear' and 'log'
        params.nContrastsPerDirection = 5;
        params.lowContrast = 0.001;
        params.highContrast = 0.15;
        params.contrastScale = 'log';  
        
        if (params.scaleIntoMonitorGamut)
            error('Cannot do scaling for contrasts instance type');
        end
    otherwise
        error('Unknown instance type passed.');
end


