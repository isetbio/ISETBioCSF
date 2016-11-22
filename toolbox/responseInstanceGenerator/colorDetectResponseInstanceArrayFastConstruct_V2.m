function [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct_V2(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
% [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
% 
% Construct an array of nTrials response instances given the
% simulationTimeStep, spatialParams, temporalParams, theOI, theMosaic.
%
% The noise free isomerizations response is returned for the first frame
% in the temporal sequence.  It is for debugging and probably not of
% general interest.
%
% This is a sped up version of colorDetectResponseInstanceArrayConstruct.
%
%  7/10/16  npc Wrote it.

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

% Compyte the modulated OI
oiModulated = theOI;
oiModulated = oiCompute(oiModulated, modulatedScene);

% Generate the stimulus modulation function
[stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(temporalParams);

% Compute the oiSequence
theOIsequence = oiSequence(oiBackground, oiModulated, stimulusTimeAxis, ...
                                stimulusModulationFunction, 'composition', 'blend');
                            
% Generate eye movement paths for all instances        
eyeMovementsNum = computeEyeMovementsNum(theMosaic.integrationTime, theOIsequence);
theEMpaths = zeros(nTrials, eyeMovementsNum, 2);     

[isomerizations, isomerizationsTimeAxis, ...
  photocurrents, photoCurrentTimeAxis] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths', theEMpaths, ...
                    'currentFlag', true, ...
                    'newNoise', true ...
                    );
      
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


function eyeMovementsNum = computeEyeMovementsNum(integrationTime, theOIsequence)
    % Generate eye movement sequence for all oi's
    if (numel(theOIsequence.oiTimeAxis) == 1)
        stimulusSamplingInterval  = integrationTime;
    else
        stimulusSamplingInterval = oiTimeAxis(2)-oiTimeAxis(1);
    end
    
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/integrationTime;
    eyeMovementsNum = round(eyeMovementsNumPerOpticalImage*theOIsequence.length);
    
    if (eyeMovementsNum < 1)
        error('Less than 1 eye movement!!! \nStimulus sampling interval:%g ms Cone mosaic integration time: %g ms\n', 1000*stimulusSamplingInterval, 1000*theConeMosaic.integrationTime);
    else 
        %fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', eyeMovementsNum, eyeMovementsNumPerOpticalImage);
    end 
end