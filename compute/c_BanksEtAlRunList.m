% c_BanksEtAlRunList
%
% Run the Banks et al. computations with various set parameters.

%% General parameters
params.computeResponses = true;
params.findPerformance = true;
params.fitPsychometric = true;
params.nTrainingSamples = 10;
params.thresholdCriterionFraction = 0.75;
params.freezeNoise = true;
params.useTrialBlocks = true;
params.generatePlots = true;

params.blur = false;
params.apertureBlur = false;
params.visualizeResponses = true;
c_BanksEtAlReplicate(params, ...
    'useScratchTopLevelDirName',true, ...
    'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ...
    'cyclesPerDegree',50,'luminances',340,'pupilDiamMm',2,...
    'nContrastsPerDirection',2,'lowContrast',0.0001,'highContrast',0.1,'contrastScale','linear');