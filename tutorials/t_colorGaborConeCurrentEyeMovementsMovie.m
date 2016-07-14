%% t_coneGaborConeCurrentEyeMovementsMovie
%
% Show how to generate a movie with the cone absoprtions and photocurrent
% to a stimulus, with eye movements and optional CRT raster effects.
%
% The ideas in this tutorial are encapsulated in expanded form
% in the routine
%   colorDetectResponseInstanceArrayConstruct
% whose use is demonstrated in the tutorial
%   t_colorGaborConeCurrentEyeMovementsResponseInstances.
%
%  7/9/16  npc Wrote it.

%% Initialize
ieInit; clear; close all;

% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'../toolbox')); 

%% Define parameters of simulation  
%
% The time step at which to compute eyeMovements and osResponses
simulationTimeStep = 5/1000;

% Stimulus (gabor) params
% 
% See t_colorGaborScene for what these means.  One additional
% parameter here is leakageLum, which is the luminance of the 
% display device when there is zero input.
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
gaborParams.leakageLum = 1.0;
gaborParams.monitorFile = 'CRT-MODEL';
gaborParams.viewingDistance = 0.75;
theBaseGaborParams = gaborParams;

% Temporal modulation and stimulus sampling parameters
frameRate = 60;
temporalParams.windowTauInSeconds = 0.165;
temporalParams.stimulusDurationInSeconds = 5*temporalParams.windowTauInSeconds;
temporalParams.stimulusSamplingIntervalInSeconds = 1/frameRate;

% Optionally, have zero amplitude eye movements
temporalParams.eyesDoNotMove = false; 

% Optional CRT raster effects.
% 
% The underlying routine that generates temporal samples 
% can simulate the fact that CRTs produce an impulse during
% each frame, although this simulation works on a frame basis
% not on a pixel-by-pixel basis.  
% 
% The parameer rasterSamples is the number
% of raster samples generated per CRT refresh
% interval.
temporalParams.addCRTrasterEffect = false;
temporalParams.rasterSamples = 5; 
if (temporalParams.addCRTrasterEffect)
    simulationTimeStep = simulationTimeStep/temporalParams.rasterSamples;
end

% Optical image parameters
oiParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
oiParams.offAxis = false;
oiParams.blur = false;
oiParams.lens = true;

% Cone mosaic parameters
mosaicParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
mosaicParams.macular = true;
mosaicParams.LMSRatio = [0.6 0.3 0.1];
mosaicParams.timeStepInSeconds = simulationTimeStep;
mosaicParams.integrationTimeInSeconds = mosaicParams.timeStepInSeconds;
mosaicParams.photonNoise = true;
mosaicParams.osNoise = true;
mosaicParams.osModel = 'Linear';

%% Create the optics
theOI = colorDetectOpticalImageConstruct(oiParams);

%% Create the cone mosaic
theMosaic = colorDetectConeMosaicConstruct(mosaicParams);

%% Create stimulus temporal window
[stimulusSampleTimes, gaussianTemporalWindow, rasterModulation] = gaussianTemporalWindowCreate(temporalParams);
if (temporalParams.addCRTrasterEffect)
    temporalParams.stimulusSamplingIntervalInSeconds = stimulusSampleTimes(2)-stimulusSampleTimes(1);
end
stimulusFramesNum = length(stimulusSampleTimes);

%% Generate eye movements for the entire stimulus duration
eyeMovementsPerStimFrame = temporalParams.stimulusSamplingIntervalInSeconds/simulationTimeStep;
eyeMovementsTotalNum = round(eyeMovementsPerStimFrame*stimulusFramesNum);
eyeMovementSequence = theMosaic.emGenSequence(eyeMovementsTotalNum);
if (isfield(temporalParams,'eyesDoNotMove') && (temporalParams.eyesDoNotMove))
    eyeMovementSequence = eyeMovementSequence * 0;
end
        
%% Loop over our stimulus frames
for stimFrameIndex = 1:stimulusFramesNum
    fprintf('Computing isomerizations for frame %d of %d\n', stimFrameIndex, stimulusFramesNum);
    
    % Modulate stimulus contrast
    gaborParams.contrast = theBaseGaborParams.contrast * gaussianTemporalWindow(stimFrameIndex);
    
    % Apply CRT raster modulation
    if (~isempty(rasterModulation))
        gaborParams.contrast = theBaseGaborParams.contrast * gaussianTemporalWindow(stimFrameIndex) * rasterModulation(stimFrameIndex);
        gaborParams.backgroundxyY(3) = gaborParams.leakageLum + theBaseGaborParams.backgroundxyY(3)*rasterModulation(stimFrameIndex);
    end
    
    % Create a scene for the current frame
    theScene = colorGaborSceneCreate(gaborParams);
    
    % Compute the optical image
    theOI = oiCompute(theOI, theScene);
    
    % Apply current frame eye movements to the mosaic
    eyeMovementIndices = (round((stimFrameIndex-1)*eyeMovementsPerStimFrame)+1 : round(stimFrameIndex*eyeMovementsPerStimFrame));
    theMosaic.emPositions = eyeMovementSequence(eyeMovementIndices,:);
    
    % Compute isomerizations for the current frame
    frameIsomerizationSequence = theMosaic.compute(theOI,'currentFlag',false);
    
    if (stimFrameIndex==1)
        coneIsomerizationSequence = frameIsomerizationSequence;
    else
        coneIsomerizationSequence = cat(3, coneIsomerizationSequence, frameIsomerizationSequence);
    end
end 

%% Compute photocurrent sequence
fprintf('Computing photocurrent sequence ...\n');
coneIsomerizationRate = coneIsomerizationSequence/theMosaic.integrationTime;
photocurrentSequence = theMosaic.os.compute(coneIsomerizationRate,theMosaic.pattern);
timeAxis = (1:size(photocurrentSequence,3))*mosaicParams.timeStepInSeconds;
timeAxis = timeAxis - (timeAxis(end)-timeAxis(1))/2;

% Visualize and render video of the isomerizations
conditionDir = paramsToDirName(gaborParams,temporalParams,oiParams,mosaicParams,[]);
visualizeMosaicResponseSequence(conditionDir, 'isomerizations (R*/cone)', coneIsomerizationSequence, eyeMovementSequence, ...
                                theMosaic.pattern, timeAxis, [theMosaic.width theMosaic.height], ...
                                theMosaic.fov, mosaicParams.integrationTimeInSeconds, ...
                                'gaborIsomerizationsWithEyeMovements');

% Visualize and render video of the photocurrents
visualizeMosaicResponseSequence(conditionDir, 'photocurrent (pAmps)', photocurrentSequence, eyeMovementSequence, ...
                                theMosaic.pattern, timeAxis, [theMosaic.width theMosaic.height], ...
                                theMosaic.fov, mosaicParams.integrationTimeInSeconds, ...
                                'gaborPhotocurrentsWithEyeMovements');
