function spatialParams = spatialParamsGenerate(varargin)
% spatialParams = spatialParamsGenerate(varargin)
%
% Define parameters of a spatial pattern.
%
% Key/value pairs
%   'spatialType' - String (default 'Gabor') Type of spatial pattern
%     'Gabor' - Gabor pattern
%     'spot' - Spatial spot on background
%     'pedestalDisk' - Disk on background
%   'windowType' - String (default 'Gaussian') Type of spatial pattern.
%     May not have an effect for all spatial types.
%     'Gaussian' - Gaussian spatial window (true Gabor)
%     'halfcos' - Half-cosine spatial window
%
% If spatialType is 'Gabor', these are the parameters:
%   windowType - What type of spatial window
%     'Gaussian' - Gaussian window, this is a real Gabor
%     'halfcos' - Half-cosine window instead.
%   fieldOfViewDegs - Field of view in degrees, linear horizontal dimension.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%                      For half-cosine, this is one half the full width.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   viewingDistance - Viewing distance of observer from monitor in meters.
%
% If spatialType is 'spot', these are the parameters
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   spotSizeDegs - Diameter of the spot, in degrees.
%   backgroundSizeDegs - Linear size of the background square, in degrees.
%   fieldOfViewDegs - Linear size of the full image, in degrees
%   row - Row dimension of image in pixels.
%   col - Col dimension of image in pixels
%   Row and colomn dimensions should be equal.

% Parse input
p = inputParser; p.KeepUnmatched = true;
p.addParameter('spatialType','Gabor',@ischar);
p.addParameter('windowType','Gaussian',@ischar);
p.parse(varargin{:});
 
% Type setup
spatialParams.type = 'Spatial';
spatialParams.spatialType = p.Results.spatialType;

switch (spatialParams.spatialType)
    case 'Gabor'
        spatialParams.windowType = p.Results.windowType;
        spatialParams.fieldOfViewDegs = 4;
        spatialParams.cyclesPerDegree = 2;
        spatialParams.gaussianFWHMDegs = 1.5;
        spatialParams.row = 128;
        spatialParams.col = 128;
        spatialParams.ang = 0;
        spatialParams.ph = 0;
        spatialParams.viewingDistance = 0.75;
        
    case 'spot'
        spatialParams.spotSizeDegs = 0.05;
        spatialParams.backgroundSizeDegs = 2;
        spatialParams.fieldOfViewDegs = 3;
        spatialParams.row = 128;
        spatialParams.col = 128;
        if (spatialParams.row ~= spatialParams.col)
            error('Image should be square');
        end
        spatialParams.viewingDistance = 0.75;
        
    case 'pedestalDisk'
        spatialParams.windowType = p.Results.windowType;
        spatialParams.windowType = 'softCircle';
        spatialParams.fieldOfViewDegs = 4;
        spatialParams.pedestalDiameterDegs = 2;
        spatialParams.row = 128;
        spatialParams.col = 128;
        spatialParams.viewingDistance = 0.75;
        
    otherwise
        error('Unknown spatial type passed: ''%s''.', spatialParams.spatialType);
end
