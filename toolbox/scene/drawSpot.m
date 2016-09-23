function [spotPattern] = drawSpot(rParams)
%this function will draw a grayscale spot (increment) on top of a
%background for simulating the spatial summation experiments of Davila and
%Geisler (1991)

%9/23/16    wst wrote it

%get spot parameters
spotRadiusDeg = rParams.spotParams.spotSizeDegs/2;
spotRadiusPixels = (rParams.spotParams.row/rParams.spotParams.fieldOfViewDegs).*spotRadiusDeg;

%create canvas in which to place the spot
spotPattern = zeros(rParams.spotParams.row, rParams.spotParams.col);

%create and place spot in canvas
croppedSpot = Circle(spotRadiusPixels);
centerRow = round(rParams.spotParams.row/2);
centerCol = round(rParams.spotParams.col/2);

if isodd(size(croppedSpot,1)) == 0; %spot diameter is even;
    minusHW = size(croppedSpot,1)/2;
    plusHW = (size(croppedSpot,1)/2)-1;
else %spot diameter is odd;
    minusHW = (size(croppedSpot,1)-1)/2;
    plusHW = minusHW;
end

spotPattern(centerRow-minusHW:centerRow+plusHW, centerCol-minusHW:centerCol+plusHW) = croppedSpot;
end