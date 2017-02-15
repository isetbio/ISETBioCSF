% c_DavilaGeislerRunList
%
% Run the Davila-Geisler summation computations with various set parameters.

% Clear and close
clear; close all;

% Common parameters
params.computeResponses = true;
params.findPerformance = true;
params.fitPsychometric = true;

params.useScratchTopLevelDirName = false;
params.thresholdMethod = 'mlpt';

% Background for production should be 2.1, but can use smaller to make
% tests run faster.
params.spotDiametersMinutes = [0.25 0.5 0.75 1 2 5 7.5 10 15 20 40];
params.backgroundSizeDegs = 1;
params.luminances = 10;

params.nTrainingSamples = 500;
params.conePacking = 'hexReg';
params.mosaicRotationDegs = 30;
params.coneSpacingMicrons = 3;
params.innerSegmentSizeMicrons = sizeForSquareApertureFromDiameterForCircularAperture(params.coneSpacingMicrons);
params.coneDarkNoiseRate = [0 0 0];
params.LMSRatio = [0.67 0.33 0];
params.pupilDiamMm = 3;
params.thresholdCriterionFraction = 0.701;
params.freezeNoise = true;
params.useTrialBlocks = true;
params.generatePlots = true;
params.plotSpatialSummation = true;
params.visualizeResponses = false;

% With blur, no dark noise, with aperture blur and Davila-Geisler optics
params.blur = true;
params.apertureBlur = true;
params.coneDarkNoiseRate = [0 0 0];
params.opticsModel = 'DavilaGeisler';
c_DavilaGeislerReplicate(params);
params = rmfield(params,'opticsModel');

