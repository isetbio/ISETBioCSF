function squareApertureDiameter = diameterForSquareApertureFromDiameterForCircularAperture(circularApertureDiameter)
    area = pi * (circularApertureDiameter/2)^2;
    squareApertureDiameter = sqrt(area);
end

