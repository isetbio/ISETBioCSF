function gaborParams = gaborParamsGenerate(varargin)
% gaborParams = gaborParamsGenerate(varargin)
%
% Define parameters of a spatial gabor pattern
%
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   viewingDistance - Viewing distance of observer from monitor in meters.

gaborParams.type = 'Gabor';

gaborParams.fieldOfViewDegs = 4;
gaborParams.cyclesPerDegree = 2;
gaborParams.gaussianFWHMDegs = 1.5;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.viewingDistance = 0.75;