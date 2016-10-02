function [spotPattern] = drawSpot(spotParams)
% [spotPattern] = drawSpot(spotParams)
% 
% This function will draw a grayscale spot (increment) on top of a
% background for simulating the spatial summation experiments of Davila and
% Geisler (1991).
%
% Background is set to zero, spot to 1 in the returend image.
%
% The input should be a spotParams structure.  See spotParamsGenerate.

% 9/23/16    wst wrote it

% Get spot parameters
spotRadiusDeg = spotParams.spotSizeDegs/2;
spotRadiusPixels = (spotParams.row/spotParams.backgroundSizeDegs).*spotRadiusDeg;

% Create canvas in which to place the spot
spotPattern = zeros(spotParams.row, spotParams.col);

% Create and place spot in canvas
croppedSpot = Circle(spotRadiusPixels);
centerRow = round(spotParams.row/2);
centerCol = round(spotParams.col/2);

if isodd(size(croppedSpot,1)) == 0; %spot diameter is even;
    minusHW = size(croppedSpot,1)/2;
    plusHW = (size(croppedSpot,1)/2)-1;
else %spot diameter is odd;
    minusHW = (size(croppedSpot,1)-1)/2;
    plusHW = minusHW;
end

spotPattern(centerRow-minusHW:centerRow+plusHW, centerCol-minusHW:centerCol+plusHW) = croppedSpot;
end