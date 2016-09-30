function colorModulationParams = colorModulationParamsGenerate(varargin)
% colorModulationParams = colorModulationParamsGenerate(varargin)
%
% Generate parameters for a color modulation.

% Key/value pairs
%   'modulationType' - String (default 'monitor') Type of modulation
%     'monitor' - Specify properties of a modulatoin on a monitor
%     'AO' - Specify parameters of adaptive optics rig stimulus
%
% If modulationType is 'monitor', these are the parameters:
%   contrast - Contrast specfied relative to coneContrasts.
%              Can be a vector of contrasts.
%   coneContrasts - Color direction of grating in cone contrast space
%                   Can be a 3 by N matrix of contrast directions.
%
% If modulationType is 'AO', these are the parameters
%   spotWavelengthNm - Vector of wavelengths of light superimposed
%     in the background.
%   spotCornealIrradianceUW - Vector of corneal irradiance for
%     spot at full power, less any light leakage, units of UW/cm2.

% Parse input
p = inputParser;
p.addParameter('modulationType','monitor',@isstring);
p.parse(varargin{:});

colorModulationParams.type = 'ColorModulation';
colorModulationParams.modulationType = p.Results.modulationType;

switch (colorModulationParams.modulationType )
    case 'monitor'
        colorModulationParams.contrast = 1;
        colorModulationParams.coneContrasts = [0.05 -0.05 0]';
    case 'AO'
        colorModulationParams.spotWavelengthNm = 680;
        colorModulationParams.spotCornealIrradianceUW = 20;
    otherwise
        error('Unknown modulation type specified');
end


