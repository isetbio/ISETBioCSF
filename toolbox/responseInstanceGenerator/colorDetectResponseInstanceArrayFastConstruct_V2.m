function [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct_V2(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
% [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
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
theEMpaths = zeros(nTrials, eyeMovementsNum, 2);     

[isomerizations, photocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths', theEMpaths, ...
                    'currentFlag', true, ...
                    'newNoise', true ...
                    );
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

% Bring the data to old format
for iTrial = 1:nTrials
    trialIsomerizations = single(squeeze(isomerizations(iTrial,:,:,:)));
    if (isempty(photocurrents))
        trialPhotocurrents = [];
    else
        trialPhotocurrents = single(squeeze(photocurrents(iTrial,:,:,:)));
    end
    if (iTrial == 1)
        theFirstInstance = struct(...
            'theMosaicIsomerizations', trialIsomerizations, ...
            'theMosaicPhotoCurrents', trialPhotocurrents, ...
            'theMosaicEyeMovements', squeeze(theEMpaths(iTrial,isomerizationsTimeIndicesToKeep,:)), ...
            'timeAxis', isomerizationsTimeAxis(isomerizationsTimeIndicesToKeep), ...
            'photocurrentTimeAxis', photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep) ...
        );
        responseInstanceArray = repmat(theFirstInstance, nTrials, 1);
        responseInstanceArray(1) = theFirstInstance;
    else
        responseInstanceArray(iTrial) = struct(...
            'theMosaicIsomerizations', trialIsomerizations, ...
            'theMosaicPhotoCurrents', trialPhotocurrents, ...
            'theMosaicEyeMovements', squeeze(theEMpaths(iTrial,isomerizationsTimeIndicesToKeep,:)), ...
            'timeAxis', isomerizationsTimeAxis(isomerizationsTimeIndicesToKeep), ...
            'photocurrentTimeAxis', photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep) ...
        );
    end
end % iTrial

% Compute the noise-free isomerizations indirectly, from the mean over all trials
% Relax, for now only ...
noiseFreeIsomerizations = squeeze(mean(isomerizations,1));

fprintf('Response instance array generation (%d instances) took %2.3f minutes to compute.\n', nTrials, toc/60);

end
