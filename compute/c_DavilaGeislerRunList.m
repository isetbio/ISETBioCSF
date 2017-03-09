% c_DavilaGeislerRunList
%
% Run the Davila-Geisler summation computations with various set parameters.

% Clear and close
clear; close all;

% Common parameters
params.computeResponses = false;
params.findPerformance = false;
params.fitPsychometric = false;

params.useScratchTopLevelDirName = false;
params.thresholdMethod = 'mlpt';

% Background for production should be 2.1, but can use smaller to make
% tests run faster.
params.spotDiametersMinutes = [0.25 0.5 0.75 1 1.5 2 3 5 7.5 10 15 20 40 60];
params.backgroundSizeDegs = 1.5;
params.luminances = 10;
 
params.mosaicRotationDegs = 30;
params.coneSpacingMicrons = 3;
params.innerSegmentSizeMicrons = sizeForSquareApertureFromDiameterForCircularAperture(params.coneSpacingMicrons);
params.coneDarkNoiseRate = [0 0 0];
params.LMSRatio = [0.67 0.33 0];
params.pupilDiamMm = 3;
params.thresholdCriterionFraction = 0.701;
params.freezeNoise = true;
params.generatePlots = true;
params.plotSpatialSummation = true;
params.visualizeResponses = false;

params.imagePixels = 1000;

% Check that pixels per minute isn't insane, given smallest spot size
% The D-G code sets the image to 1.1 times the background size at present
imageSizeMinutes = 60*1.1*params.backgroundSizeDegs;
minutesPerPixel = imageSizeMinutes/params.imagePixels;
fprintf('Smallest spot diameter %g minutes, minutesPerPixel %g, pixels per spot %0.1f\n',min(params.spotDiametersMinutes),minutesPerPixel,min(params.spotDiametersMinutes)/minutesPerPixel);

% With blur, no dark noise, with aperture blur and Davila-Geisler optics
params.blur = true;
params.apertureBlur = true;
params.coneDarkNoiseRate = [0 0 0];
params.opticsModel = 'DavilaGeisler';
c_DavilaGeislerReplicate(params);
params = rmfield(params,'opticsModel');

% No blur, no dark noise, with aperture blur
params.blur = false;
params.apertureBlur = true;
params.coneDarkNoiseRate = [0 0 0];
params.opticsModel = 'DavilaGeisler';
c_DavilaGeislerReplicate(params);
params = rmfield(params,'opticsModel');

% No blur, no dark noise, with aperture blur
params.blur = true;
params.apertureBlur = true;
params.coneDarkNoiseRate = [0 0 0];
params.opticsModel = 'DavilaGeislerLsfAsPsf';
c_DavilaGeislerReplicate(params);
params = rmfield(params,'opticsModel');

