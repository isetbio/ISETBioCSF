% c_BanksEtAlRunList
%
% Run the Banks et al. computations with various set parameters.

% Common parameters
params.useScratchTopLevelDirName = false;
params.computeResponses = true;
params.findPerformance = true;
params.fitPsychometric = true;
params.nTrainingSamples = 1000;
params.conePacking = 'hexReg';
params.mosaicRotationDegs = 30;
params.coneSpacingMicrons = 3;
params.innerSegmentSizeMicrons = sizeForSquareApertureFromDiameterForCircularAperture(params.coneSpacingMicron);
params.cyclesPerDegree = [3 5 10 20 40 50 60];
params.luminances = [3.4 34 340];
params.pupilDiamMm = 2;
params.thresholdCriterionFraction = 0.75;
params.freezeNoise = true;
params.useTrialBlocks = true;
params.generatePlots = true;
params.visualizeResponses = false;

% No blur, no dark noise, no aperture blur
params.blur = false;
params.apertureBlur = false;
params.darkNoise = [0 0 0];
c_BanksEtAlReplicate(params);

% No blur, no dark noise, with aperture blur
params.blur = false;
params.apertureBlur = true;
params.darkNoise = [0 0 0];
c_BanksEtAlReplicate(params);

% With blur, no dark noise, with aperture blur
params.blur = true;
params.apertureBlur = true;
params.darkNoise = [0 0 0];
c_BanksEtAlReplicate(params);
