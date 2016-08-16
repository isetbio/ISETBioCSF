function isomerizations = IsomerizationsFromLuminanceGeisler(luminanceCdM2,durationInSeconds,pupilDiameterMm)
% isomerizations = IsomerizationsFromLuminanceGeisler(luminanceCdM2,durationInSeconds,pupilDiameterMm)
%
% Geisler 1984 gives a formula for isomerizations from luminance (his
% equation 2), which he obtained from Wyszecki and Stiles.
%
% This function computes this with the 1984 parameters, ignoring
% convolution by the optical point spread function but including lens
% transmittance.
%
% A number of parameters from Geisler's paper are hard coded into this
% routine, but one could add key value pairs to allow them to be set.

% Geisler's choice of occular transimttance parameter
occularTransmittance = 0.68;

% Geisler's choice of quantum efficiency at 555 nm parameter
quantumEfficiency555 = 0.5;

receptorApertureMinutes = 0.6;
receptorAreaMinutes2 = pi*((receptorApertureMinutes/2)^2);

% Pupil area
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);

% This is a magic unit constant in the formula
magicUnitConstant = 347.8;

%% Apply the formula
isomerizations = receptorAreaMinutes2*durationInSeconds*pupilAreaMm2*occularTransmittance*quantumEfficiency555*magicUnitConstant*luminanceCdM2;