function backgroundParams = backgroundParamsGenerate(varargin)
% backgroundParams = backgroundParamsGenerate(varargin)
%
% Generate background spectral parameters
%
% Key/value pairs
%   'backgroundType' - String (default 'monitor') Type of background
%     'monitor' - Specify properties of a background on a monitor
%     'AO' - Specify parameters of adaptive optics rig background
%
% If backgroundType is 'monitor', these are the parameters:
%   backgroundxYY - Colorimetric specification of background, in CIE xyY (cd/m2)
%   monitorFile - Isetbio display description of monitor on which grating is shown.
%   leakageLum - Luminance when monitor input is zero.
%   lumFactor - Multiply background luminance by this factor.
%
% If backgroundType is 'AO', these are the parameters
%   backgroundWavelengthsNm - Vector of wavelengths of light superimposed
%     in the background.
%   backgroundCornealPowerUW - Vector of corneal irradiance for each
%     of the monochromatic lights in the background, units of UW/cm2.

% Parse input
p = inputParser; p.KeepUnmatched = true;
p.addParameter('backgroundType','monitor',@ischar);
p.parse(varargin{:});

backgroundParams.type = 'Background';
backgroundParams.backgroundType = p.Results.backgroundType;

switch (backgroundParams.backgroundType)
    case 'monitor'
        backgroundParams.backgroundxyY = [0.27 0.30 49.8]';
        backgroundParams.monitorFile = 'CRT-MODEL';
        backgroundParams.leakageLum = 1.0;
        backgroundParams.lumFactor = 1;
    case 'AO'
        backgroundParams.backgroundWavelengthsNm = [830 790 680];
        backgroundParams.backgroundCornealPowerUW = [20 10 0.02];
    otherwise
        error('Unknown background type specified');
end