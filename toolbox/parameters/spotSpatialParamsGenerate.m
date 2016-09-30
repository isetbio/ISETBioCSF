function spotSpatialParams = spotSpatialParamsGenerate(varargin)
% spotSpatialParams = spotSpatialParamsGenerate(varargin)
%
% Define spatial parameters of a simple circular spot on a background.
% Area outside of background is assumed to be black.
%
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   spotSizeDegs - Diameter of the spot, in degrees.
%   backgroundSizeDegs - Linear size of the background square, in degrees.
%   imageSizeDegs - Linear size of the image, in degrees
%   row - Row dimension of image in pixels.
%   col - Col dimension of image in pixels
%
% Row and coloum dimensions should be equal.

spotSpatialParams.type = 'SpotSpatial';

% Whole field of view
spotSpatialParams.viewingDistance = 0.75;

% The spot
spotSpatialParams.spotSizeDegs = 1;

% Background parameters
spotSpatialParams.backgroundSizeDegs = 2;

% Image size 
spotSpatialParams.imageSizeDegs

% Pixel resolution of scene image
spotSpatialParams.row = 128;
spotSpatialParams.col = 128;
if (spotSpatialParams.row ~= spotSpatialParams.col)
    error('Image should be square');
end