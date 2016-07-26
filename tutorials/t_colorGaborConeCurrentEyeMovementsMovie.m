function validationData = t_coneGaborConeCurrentEyeMovementsMovie(rParams)
% validationData = t_coneGaborConeCurrentEyeMovementsMovie(rParams)
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
% The returned validation structure allows this routine to be called from a
% validation script driven by the UnitTest toolbox.
%
% The tutorial produces output according to a scheme controlled by the
% specified IBIOColorDetect rwObject.
%
% See also:
%   t_colorGaborRespnseGenerationParams
%   t_colorGaborScene
%	t_colorGaborConeIsomerizationsMovie
%   t_colorGaborConeCurrentEyeMovementsResponseInstances
%   colorGaborSceneCreate 
%   colorDetectOpticalImageConstruct
%   colorDetectConeMosaicConstruct
%   colorDetectResponseInstanceArrayConstruct
%   colorDetectResponseInstanceArrayFastConstruct
%
%  7/9/16  npc Wrote it.

%% Initialize
ieInit; clear; close all;

%% Fix random number generator so we can validate output exactly
rng(1);

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (nargin < 1 | isempty(rParams))
    rParams = t_colorGaborResponseGenerationParams;
end

% Override some of the defaults
rParams.mosaicParams.isomerizationNoise = true;
rParams.mosaicParams.osNoise = true;
rParams.mosaicParams.osModel = 'Linear';

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
rwObject.readProgram = '';
rwObject.writeProgram = mfilename;
rwObject.parentParamsList = {};
rwObject.currentParamsList = {rParams, rParams.colorModulationParams};

%% Create the optics
theOI = colorDetectOpticalImageConstruct(rParams.oiParams);

%% Create the cone mosaic
theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);

%% Create stimulus temporal window
[stimulusSampleTimes, gaussianTemporalWindow, rasterModulation] = gaussianTemporalWindowCreate(rParams.temporalParams);
if (rParams.temporalParams.addCRTrasterEffect)
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = stimulusSampleTimes(2)-stimulusSampleTimes(1);
end
stimulusFramesNum = length(stimulusSampleTimes);

%% Generate eye movements for the entire stimulus duration
eyeMovementsPerStimFrame = rParams.temporalParams.stimulusSamplingIntervalInSeconds/rParams.temporalParams.simulationTimeStepSecs;
eyeMovementsTotalNum = round(eyeMovementsPerStimFrame*stimulusFramesNum);
eyeMovementSequence = theMosaic.emGenSequence(eyeMovementsTotalNum);
if (isfield(rParams.temporalParams,'eyesDoNotMove') && (rParams.temporalParams.eyesDoNotMove))
    eyeMovementSequence = eyeMovementSequence * 0;
end
        
%% Loop over our stimulus frames
baseColorModulationParams = rParams.colorModulationParams;
for stimFrameIndex = 1:stimulusFramesNum
    fprintf('Computing isomerizations for frame %d of %d\n', stimFrameIndex, stimulusFramesNum);
    
    % Modulate stimulus contrast
    colorModulationParamsTemp = baseColorModulationParams;
    colorModulationParamsTemp.contrast = baseColorModulationParams.contrast * gaussianTemporalWindow(stimFrameIndex);
    
    % Apply CRT raster modulation
    if (~isempty(rasterModulation))
        colorModulationParamsTemp.contrast = theBasecolorModulationParams.contrast * gaussianTemporalWindow(stimFrameIndex) * rasterModulation(stimFrameIndex);
        colorModulationParamsTemp.backgroundxyY(3) = gaborParams.leakageLum + theBasecolorModulationParams.backgroundxyY(3)*rasterModulation(stimFrameIndex);
    end
    
    % Create a scene for the current frame
    theScene = colorGaborSceneCreate(rParams.gaborParams,colorModulationParamsTemp);
    
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
timeAxis = (1:size(photocurrentSequence,3))*rParams.mosaicParams.timeStepInSeconds;
timeAxis = timeAxis - (timeAxis(end)-timeAxis(1))/2;

%% Visualize and render video of the isomerizations
visualizeMosaicResponseSequence(rwObject, 'isomerizations (R*/cone)', coneIsomerizationSequence, eyeMovementSequence, ...
                                theMosaic.pattern, timeAxis, [theMosaic.width theMosaic.height], ...
                                theMosaic.fov, rParams.mosaicParams.integrationTimeInSeconds, ...
                                'gaborIsomerizationsWithEyeMovements');

%% Visualize and render video of the photocurrents
visualizeMosaicResponseSequence(rwObject, 'photocurrent (pAmps)', photocurrentSequence, eyeMovementSequence, ...
                                theMosaic.pattern, timeAxis, [theMosaic.width theMosaic.height], ...
                                theMosaic.fov, rParams.mosaicParams.integrationTimeInSeconds, ...
                                'gaborPhotocurrentsWithEyeMovements');

%% Return validation if desired                           
if (nargout > 0)
    validationData = [];
end
