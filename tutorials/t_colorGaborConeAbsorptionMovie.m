% t_colorGaborConeAbsorptionsMovie
%
% Create scene sequence for a Gaussian windowed color gabor and then from
% it generate an optical image sequence and finally a cone reponse
% movie.
%
% The scene sequence generation logic illustrated here is encapsulated in a
% fancier manner in colorGaborSceneSequenceCreate.
%
% See also t_colorGaborScene, colorGaborSceneSequenceCreate
%
% 7/8/16  dhb  Wrote it.

%% Initialize
ieInit; clear; close all;

% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'../toolbox')); 

%% Define parameters of a gabor pattern
%
% Parameters in degrees.  The field of view is the horizontal dimension.
gaborParams.fieldOfViewDegs = 1.5;
gaborParams.gaussianFWHMDegs = 0.7;
gaborParams.cyclesPerDegree = 2;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.contrast = 1;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.coneContrasts = [0.06 -0.06 0]';
gaborParams.backgroundxyY = [0.27 0.30 49.8]';
gaborParams.monitorFile = 'CRT-MODEL';
gaborParams.viewingDistance = 0.75;

% Temporal stimulus parameters
%
% 50 msec sampling intervals are probably too long for real work, but OK
% for this tutorial.
frameRate = 60;
temporalParams.windowTauInSeconds = 0.165;
temporalParams.stimulusDurationInSeconds = 5*temporalParams.windowTauInSeconds;
temporalParams.stimulusSamplingIntervalInSeconds = 1/frameRate;
[sampleTimes,gaussianTemporalWindow] = gaussianTemporalWindowCreate(temporalParams);
nSampleTimes = length(sampleTimes);

% Plot the temporal window, just to make sure it looks right
figure(1); clf;
plot(sampleTimes,gaussianTemporalWindow,'r');
xlabel('Time (seconds)');
ylabel('Window Amplitude');
title('Stimulus Temporal Window');

%% Loop over time and build a cell array of scenes
gaborScene = cell(nSampleTimes,1);
for ii = 1:nSampleTimes
    % Make the sence for this time
    gaborParams.contrast = gaussianTemporalWindow(ii);
    fprintf('Computing scene %d of %d, time %0.3f, windowVal %0.3f\n',ii,nSampleTimes,sampleTimes(ii),gaussianTemporalWindow(ii));
    gaborScene{ii} = colorGaborSceneCreate(gaborParams);
end

%% Create the OI object we'll use to compute the retinal images from the scenes
%
% Then oop over scenes and compute the optical image for each one 
oiParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
oiParams.offAxis = false;
oiParams.blur = false;
oiParams.lens = true;
theBaseOI = colorDetectOpticalImageConstruct(oiParams);

theOI = cell(nSampleTimes,1);
for ii = 1:nSampleTimes
    % Compute retinal image
    fprintf('Computing optical image %d of %d, time %0.3f\n',ii,nSampleTimes,sampleTimes(ii));
    theOI{ii} = oiCompute(theBaseOI,gaborScene{ii});
end

%% Create the coneMosaic object we'll use to compute cone respones
mosaicParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
mosaicParams.macular = true;
mosaicParams.LMSRatio = [0.6 0.3 0.1];
mosaicParams.timeStepInSeconds = temporalParams.stimulusSamplingIntervalInSeconds;
mosaicParams.integrationTimeInSeconds = mosaicParams.timeStepInSeconds;
mosaicParams.photonNoise = false;
mosaicParams.osModel = 'Linear';
theMosaic = colorDetectConeMosaicConstruct(mosaicParams);

for ii = 1:nSampleTimes      
    % Compute mosaic response for each stimulus frame
    % For a real calculation, we would save these so that we could use them
    % to do something.  But here we just (see next line) compute the
    % contrast seen by each class of cone in the mosaic, just to show we
    % can do something.
    fprintf('Computing absorptions %d of %d, time %0.3f\n',ii,nSampleTimes,sampleTimes(ii));
    gaborConeAbsorptions(:,:,ii) = theMosaic.compute(theOI{ii},'currentFlag',false);
end

%% Make a movie of the stimulus sequence
conditionDir = paramsToDirName(gaborParams,temporalParams,oiParams,mosaicParams,[]);
showLuminanceMap = false;
visualizeSceneOrOpticalImageSequence(conditionDir, 'scene', gaborScene, sampleTimes, showLuminanceMap, 'gaborStimulusMovie');

%% Make a movie of the stimulus sequence
showLuminanceMap = false;
visualizeSceneOrOpticalImageSequence(conditionDir,'optical image', theOI, sampleTimes, showLuminanceMap, 'gaborOpticalImageMovie');

%% Make a movie of the isomerizations
eyeMovementSequence = [];
visualizeMosaicResponseSequence(conditionDir, 'isomerizations (R*/cone)', gaborConeAbsorptions, eyeMovementSequence, theMosaic.pattern, sampleTimes, [theMosaic.width theMosaic.height], theMosaic.fov, mosaicParams.integrationTimeInSeconds, 'gaborIsomerizations');

%% Plot cone contrasts as a function of time, as a check
%
% We are not quite sure why they have the scalloped look that they do.
% Maybe monitor quantization?
for ii = 1:nSampleTimes      
    LMSContrasts(:,ii) = mosaicUnsignedConeContrasts(gaborConeAbsorptions(:,:,ii),theMosaic);
end
vcNewGraphWin; hold on;
plot(sampleTimes,LMSContrasts(1,:)','r');
plot(sampleTimes,LMSContrasts(2,:)','g');
plot(sampleTimes,LMSContrasts(3,:)','b');
xlabel('Time (seconds)');
ylabel('Contrast');
title('LMS Cone Contrasts');

