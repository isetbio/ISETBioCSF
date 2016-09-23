function spotParams = spotParamsGenerate(varargin)
% spotParams = spotParamsGenerate(varargin)
%
% Define parameters of a spatial gabor pattern.  Actually, the spatial
% window can also be a half-cosine.
%
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%                      For half-cosine, this is one half the full width.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   viewingDistance - Viewing distance of observer from monitor in meters.

spotParams.type = 'Spot';

% Whole field of view
spotParams.viewingDistance = 0.75;

% The spot
spotParams.spotSizeDegs = 1.5;
spotParams.spotWavelengthNm = 680;
spotParams.spotCornealIrradianceUW = 20;

% Background parameters
spotParams.backgroundSizeDegs = 2;
spotParams.backgroundWavelengthsNm = [830 790 680];
spotParams.backgroundCornealIrradianceUW = [20 10 20];

% Pixel resolution of scene image
spotParams.row = 128;
spotParams.col = 128;
