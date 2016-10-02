function [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, simulationTimeStep, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
% [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, simulationTimeStep, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
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

% Inform user regarding the computation progress
% progressHandle = generateProgressBar('Starting computation ...');

% Start computation time measurement
tic

% Save base color modulation params
theBaseColorModulationParams = colorModulationParams;

% Create stimulus temporal window
[stimulusSampleTimes, gaussianTemporalWindow, rasterModulation] = gaussianTemporalWindowCreate(temporalParams);
if (temporalParams.addCRTrasterEffect)
    temporalParams.stimulusSamplingIntervalInSeconds = stimulusSampleTimes(2)-stimulusSampleTimes(1);
end
stimulusFramesNum = length(stimulusSampleTimes);

% Generate new, larger mosaic to deal with eye movements
padRows = 64;
padCols = 64;
theLargerMosaic = theMosaic.copy();
theLargerMosaic.pattern = zeros(theMosaic.rows+2*padRows, theMosaic.cols+2*padCols);
theLargerMosaic.noiseFlag = false;

% Loop over our stimulus frames
for stimFrameIndex = 1:stimulusFramesNum
    % waitbar(0.5*stimFrameIndex/stimulusFramesNum, progressHandle, sprintf('stimulus label: %s\ncomputing optical image for frame #%d/%d', stimulusLabel, stimFrameIndex, stimulusFramesNum));
    
    % Modulate stimulus contrast
    colorModulationParams.contrast = theBaseColorModulationParams.contrast * gaussianTemporalWindow(stimFrameIndex);
    
    % % Apply CRT raster modulation
    % if (~isempty(rasterModulation))
    %     colorModulationParams.contrast = theBaseColorModulationParams.contrast * gaussianTemporalWindow(stimFrameIndex) * rasterModulation(stimFrameIndex);
    %     backgroundParams.backgroundxyY(3) = backgroundParams.leakageLum + backgroundParams.backgroundxyY(3)*rasterModulation(stimFrameIndex);
    % end
    
    % Create a scene for the current frame
    theScene = colorSceneCreate(spatialParams,backgroundParams,colorModulationParams);
    
    % Compute the optical image for the current frame
    theOI = oiCompute(theOI, theScene);
    
    % Compute noise-free isomerizations at each cone location for the current frame
    theFrameFullMosaicIsomerizatios{stimFrameIndex} = theLargerMosaic.computeSingleFrame(theOI,'FullLMS',true);
end

% Larger mosaic not needed anymore
clearvars('theLargerMosaic');

% For each trial compute new eye movement path and obtain new response
for iTrial = 1:nTrials
    if (mod(iTrial-1,50) == 0)
        % waitbar(0.5+0.5*iTrial/nTrials, progressHandle, sprintf('stimulus label: %s\ncomputing responses for trial %d/%d', stimulusLabel, iTrial, nTrials));
    end
    
    % Generate eye movements for the entire stimulus duration of this trial
    eyeMovementsPerStimFrame = temporalParams.stimulusSamplingIntervalInSeconds/simulationTimeStep;
    eyeMovementsTotalNum = round(eyeMovementsPerStimFrame*stimulusFramesNum);
    
    % Temporary fix, don't call emGenSequence for just one eye position
    % case.  And, worry about why emGenSequence changes integration time
    % before putting it back.
    if (eyeMovementsTotalNum == 1)
        eyeMovementSequence = zeros(1,2);
    else
        eyeMovementSequence = theMosaic.emGenSequence(eyeMovementsTotalNum);
        if (isfield(temporalParams,'eyesDoNotMove') && (temporalParams.eyesDoNotMove))
            eyeMovementSequence = eyeMovementSequence * 0;
        end
    end
    
    % Loop over our stimulus frames
    for stimFrameIndex = 1:stimulusFramesNum
        % Apply current frame eye movements to the mosaic
        eyeMovementIndices = (round((stimFrameIndex-1)*eyeMovementsPerStimFrame)+1 : round(stimFrameIndex*eyeMovementsPerStimFrame));
        theMosaic.emPositions = eyeMovementSequence(eyeMovementIndices,:);
        
        % Compute noise-free isomerizations for the current frame by applying eye movements during this stimulus frame
        theFrameEyeMovementPathIsomerizations = ...
            theMosaic.applyEMPath(theFrameFullMosaicIsomerizatios{stimFrameIndex}, ...
            'padRows',padRows,'padCols',padCols);
        
        % Stash noise free frame on first sequence, for debugging return
        if (stimFrameIndex == 1)
            noiseFreeIsomerizations = theFrameEyeMovementPathIsomerizations;
        end
        
        % Add noise
        if (theMosaic.noiseFlag)
            theFrameEyeMovementPathIsomerizations = theMosaic.photonNoise(theFrameEyeMovementPathIsomerizations,'newNoise', true);
        end
        
        % Accumulate isomerizations by adding current frame isomerizations in the current eye movement path
        if (stimFrameIndex==1)
            coneIsomerizationSequence = theFrameEyeMovementPathIsomerizations;
        else
            coneIsomerizationSequence = cat(3, coneIsomerizationSequence, theFrameEyeMovementPathIsomerizations);
        end
    end
    
    % Compute photocurrent sequence
    % 
    % Whether noise is produced here is determined by by the
    % mosaicParams.osNoise flag, which is then used to set the os noise
    % flag at mosaic initialization time.
    coneIsomerizationRate = coneIsomerizationSequence/theMosaic.integrationTime;
    photocurrentSequence = theMosaic.os.compute(coneIsomerizationRate,theMosaic.pattern);
   
    % Accumulate data in cell array of structs.
    if (iTrial == 1)
        % Only include the central response
        timeAxis = (0:size(photocurrentSequence,3)-1)*theMosaic.integrationTime;
        timeAxis = timeAxis - timeAxis(end)/2;
        timeIndicesToKeep = find(abs(timeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);

        theFirstInstance = struct(...
            'theMosaicIsomerizations', single(coneIsomerizationSequence(:,:,timeIndicesToKeep)), ...
            'theMosaicPhotoCurrents', single(photocurrentSequence(:,:,timeIndicesToKeep)), ...
            'theMosaicEyeMovements', eyeMovementSequence(timeIndicesToKeep,:), ...
            'timeAxis', timeAxis(timeIndicesToKeep) ...
        );
        responseInstanceArray = repmat(theFirstInstance, nTrials, 1);
        responseInstanceArray(1) = theFirstInstance;
    else
        responseInstanceArray(iTrial) = struct(...
            'theMosaicIsomerizations', single(coneIsomerizationSequence(:,:,timeIndicesToKeep)), ...
            'theMosaicPhotoCurrents', single(photocurrentSequence(:,:,timeIndicesToKeep)), ...
            'theMosaicEyeMovements', eyeMovementSequence(timeIndicesToKeep,:), ...
            'timeAxis', timeAxis(timeIndicesToKeep) ...
            );
    end
    
end 

fprintf('Response instance array generation (%d instances) took %2.3f minutes to compute.\n', nTrials, toc/60);

% Close progress bar
% close(progressHandle);

end

function progressHandle = generateProgressBar(initialMessage)

% progressHandle = waitbar(0, '');
% titleHandle = get(findobj(progressHandle,'Type','axes'),'Title');
% set(titleHandle,'FontSize',12, 'FontName', 'Menlo');
% waitbar(0, progressHandle, initialMessage);
% pause(0.2);
    
end
