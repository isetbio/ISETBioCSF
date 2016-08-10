function gaborParams = gaborParamsGenerate(varargin)
% gaborParams = gaborParamsGenerate(varargin)
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

gaborParams.type = 'Gabor';

gaborParams.windowType = 'Gaussian';
gaborParams.fieldOfViewDegs = 4;
gaborParams.cyclesPerDegree = 2;
gaborParams.gaussianFWHMDegs = 1.5;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.viewingDistance = 0.75;