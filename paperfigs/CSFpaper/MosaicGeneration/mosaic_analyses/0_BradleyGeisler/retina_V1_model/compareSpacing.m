function compareSpacing

    eccDegs = 2;
    angle = 180;
    micronsPerDegree = 300;
    
    eccMeters = eccDegs * micronsPerDegree * 1e-6;
    
    xDegs = eccDegs*cos(angle);
    yDegs = eccDegs*sin(angle);
    
    ganglionSpacing = spacing_fn(xDegs, yDegs);
    coneSpacing = coneSizeReadData('eccentricity', eccMeters, 'angle', angle)*1e6/micronsPerDegree;
    
    [coneSpacing ganglionSpacing]
    
    
end

