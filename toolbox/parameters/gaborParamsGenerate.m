function spatialParams = spatialParamsGenerate(varargin)
% spatialParams = spatialParamsGenerate(varargin)
%
% Define parameters of a spatial gabor pattern.  Actually, the spatial
% window can also be a half-cosine.
%
%   windowType - What type of spatial window
%     'Gaussian' - Gaussian window, this is a real Gabor
%     'halfcos' - Half-cosine window instead.
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%                      For half-cosine, this is one half the full width.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   viewingDistance - Viewing distance of observer from monitor in meters.

spatialParams.type = 'Gabor';

spatialParams.windowType = 'Gaussian';
spatialParams.fieldOfViewDegs = 4;
spatialParams.cyclesPerDegree = 2;
spatialParams.gaussianFWHMDegs = 1.5;
spatialParams.row = 128;
spatialParams.col = 128;
spatialParams.ang = 0;
spatialParams.ph = 0;
spatialParams.viewingDistance = 0.75;