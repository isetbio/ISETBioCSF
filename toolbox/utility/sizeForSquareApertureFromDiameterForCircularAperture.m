function squareApertureSize = sizeForSquareApertureFromDiameterForCircularAperture(circularApertureDiameter)
% squareApertureDiameter = sizeForSquareApertureFromDiameterForCircularAperture(circularApertureDiameter)
% 
% Given the diameter of a circular aperature, find the linear size of a
% square aperture that produces the same area. We use this because in
% isetbio, we specify photorceptor diameter by height and width of a square
% collecting aperature that are then multiplied to give us the area.

    area = pi * (circularApertureDiameter/2)^2;
    squareApertureSize = sqrt(area);
end

