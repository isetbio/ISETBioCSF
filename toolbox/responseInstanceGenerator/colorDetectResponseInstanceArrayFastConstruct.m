function [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic, varargin)
% [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
% 
% Construct an array of nTrials response instances given the
% simulationTimeStep, spatialParams, temporalParams, theOI, theMosaic.
%
% The noise free isomerizations response is returned for the first frame
% in the temporal sequence.  It is for debugging and probably not of general interest.
%
% This is a sped up version of colorDetectResponseInstanceArrayConstruct.
%
% 7/10/16  npc Wrote it.

p = inputParser;
p.addParameter('seed',1, @isnumeric);   
p.parse(varargin{:});
currentSeed = p.Results.seed;

% Start computation time measurement
tic

% Get pupil size out of OI, which is sometimes needed by colorSceneCreate
oiParamsTemp.pupilDiamMm = 1000*opticsGet(oiGet(theOI,'optics'),'aperture diameter');

% Create the background scene
theBaseColorModulationParams = colorModulationParams;
theBaseColorModulationParams.coneContrasts = [0 0 0]';
theBaseColorModulationParams.contrast = 0;
backgroundScene = colorSceneCreate(spatialParams, backgroundParams, ...
    theBaseColorModulationParams, oiParamsTemp);

% Create the modulated scene
modulatedScene = colorSceneCreate(spatialParams, backgroundParams, ...
    colorModulationParams, oiParamsTemp);

% Compute the background OI
oiBackground = theOI;
oiBackground = oiCompute(oiBackground, backgroundScene);

% Compute the modulated OI
oiModulated = theOI;
oiModulated = oiCompute(oiModulated, modulatedScene);

% Generate the stimulus modulation function
[stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(temporalParams);

% Compute the oiSequence
theOIsequence = oiSequence(oiBackground, oiModulated, stimulusTimeAxis, ...
                                stimulusModulationFunction, 'composition', 'blend');
                            
% Generate eye movement paths for all instances
eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theMosaic.integrationTime);
theEMpaths = colorDetectMultiTrialEMPathGenerate(theMosaic, nTrials, eyeMovementsNum, temporalParams.emPathType, 'seed', currentSeed);

[isomerizations, photocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'seed', currentSeed, ...
                    'emPaths', theEMpaths, ...
                    'currentFlag', true);
isomerizationsTimeAxis = theMosaic.timeAxis + theOIsequence.timeAxis(1);
photoCurrentTimeAxis = isomerizationsTimeAxis;
 
% Remove unwanted portions of the responses
isomerizationsTimeIndicesToKeep = find(abs(isomerizationsTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
photocurrentsTimeIndicesToKeep = find(abs(photoCurrentTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);

if (isa(theMosaic , 'coneMosaicHex'))
     isomerizations = isomerizations(:,:,isomerizationsTimeIndicesToKeep);
     photocurrents = photocurrents(:,:,photocurrentsTimeIndicesToKeep);
else
     isomerizations = isomerizations(:,:,:,isomerizationsTimeIndicesToKeep);
     photocurrents = photocurrents(:,:,:,photocurrentsTimeIndicesToKeep);
end

% Form the responseInstanceArray struct
responseInstanceArray.theMosaicIsomerizations = single(isomerizations);
responseInstanceArray.theMosaicPhotoCurrents = single(photocurrents);
responseInstanceArray.theMosaicEyeMovements = theEMpaths(:,isomerizationsTimeIndicesToKeep,:);
responseInstanceArray.timeAxis = isomerizationsTimeAxis(isomerizationsTimeIndicesToKeep);
responseInstanceArray.photocurrentTimeAxis = photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep);

% Compute the noise-free isomerizations & photocurrents using the first emPath
theEMpaths = theEMpaths(1,:,:);

% Noise-free responses
% Save original noise flags
originalIsomerizationNoiseFlag = theMosaic.noiseFlag;
originalPootocurrentNoiseFlag = theMosaic.os.noiseFlag;

% Set noiseFlags to none
theMosaic.noiseFlag = 'none';
theMosaic.os.noiseFlag = 'none';

% Compute the noise-free responses
[noiseFreeIsomerizations, noiseFreePhotocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths', theEMpaths, ...
                    'currentFlag', true);
% Restore noiseFlags to none
theMosaic.noiseFlag = originalIsomerizationNoiseFlag;
theMosaic.os.noiseFlag = originalPootocurrentNoiseFlag;

noiseFreeIsomerizations = squeeze(noiseFreeIsomerizations);
if (~isempty(noiseFreePhotocurrents))
    noiseFreePhotocurrents = squeeze(noiseFreePhotocurrents);
end

fprintf('Response instance array generation (%d instances) took %2.3f minutes to compute.\n', nTrials, toc/60);

end
