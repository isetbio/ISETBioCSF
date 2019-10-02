function compareSpacing
% compare the ganglion cell spacing.
%
% Syntax:
%   compareSpacing
%
% Description:
%    Compare the ganglion cell spacing using the provided data of:
%       eccentricity: 2 degrees
%       angle: 180 degrees
%       microns per degree: 300.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    eccDegs = 2;
    angle = 180;
    micronsPerDegree = 300;
    
    eccMeters = eccDegs * micronsPerDegree * 1e-6;
    
    xDegs = eccDegs * cos(angle);
    yDegs = eccDegs * sin(angle);
    
    ganglionSpacing = spacing_fn(xDegs, yDegs);
    coneSpacing = coneSizeReadData('eccentricity', eccMeters, ...
        'angle', angle) * 1e6 / micronsPerDegree;

    [coneSpacing ganglionSpacing]

end
