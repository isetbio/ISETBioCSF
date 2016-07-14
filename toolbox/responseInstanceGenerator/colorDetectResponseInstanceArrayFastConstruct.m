function responseInstanceArray = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, simulationTimeStep, gaborParams, temporalParams, theOI, theMosaic)
% responseInstanceArray = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, simulationTimeStep, gaborParams, temporalParams, theOI, theMosaic)
% 
% Construct an array of nTrials response instances given the simulationTimeStep, gaborParams, temporalParams, theOI, theMosaic
%
%  7/10/16  npc Wrote it.

    % Inform user regarding the computation progress
    progressHandle = generateProgressBar('Starting computation ...');
    
    % Start computation time measurement
    tic
    
    % Save base gabor params
    theBaseGaborParams = gaborParams;
    
    % Create stimulus temporal window
    [stimulusSampleTimes, gaussianTemporalWindow, rasterModulation] = gaussianTemporalWindowCreate(temporalParams);
    if (temporalParams.addCRTrasterEffect)
        temporalParams.stimulusSamplingIntervalInSeconds = stimulusSampleTimes(2)-stimulusSampleTimes(1);
    end
    stimulusFramesNum = length(stimulusSampleTimes);
    
    % Generate new, larger mosaic todeal with eye movements
    padRows = 64;
    padCols = 64;
    theLargerMosaic = theMosaic.copy();
    theLargerMosaic.pattern = zeros(theMosaic.rows+2*padRows, theMosaic.cols+2*padCols);
    theLargerMosaic.noiseFlag = false;        
            
    % Loop over our stimulus frames
    for stimFrameIndex = 1:stimulusFramesNum 
        waitbar(0.5*stimFrameIndex/stimulusFramesNum, progressHandle, sprintf('stimulus label: %s\ncomputing optical image for frame #%d/%d', stimulusLabel, stimFrameIndex, stimulusFramesNum));
        
        % Modulate stimulus contrast
        gaborParams.contrast = theBaseGaborParams.contrast * gaussianTemporalWindow(stimFrameIndex);
        
        % Apply CRT raster modulation
        if (~isempty(rasterModulation))
            gaborParams.contrast = theBaseGaborParams.contrast * gaussianTemporalWindow(stimFrameIndex) * rasterModulation(stimFrameIndex);
            gaborParams.backgroundxyY(3) = gaborParams.leakageLum + theBaseGaborParams.backgroundxyY(3)*rasterModulation(stimFrameIndex);
        end
        
        % Create a scene for the current frame
        theScene = colorGaborSceneCreate(gaborParams);   
        
        % Compute the optical image for the current frame
        theOI = oiCompute(theOI, theScene);
        
        % Compute noise-free isomerizations at each cone location for the current frame
        theFrameFullMosaicIsomerizatios{stimFrameIndex} = theLargerMosaic.computeSingleFrame(theOI,'FullLMS',true);
    end % stimFrameIndex
    
    % Larger mosaic not needed anymore
    clearvars('theLargerMosaic');
    
    % For each trial compute new eye movement path and obtain new response
    for iTrial = 1: nTrials
        waitbar(0.5+0.5*iTrial/nTrials, progressHandle, sprintf('stimulus label: %s\ncomputing responses for trial %d/%d', stimulusLabel, iTrial, nTrials));
        
        % Generate eye movements for the entire stimulus duration of this trial
        eyeMovementsPerStimFrame = temporalParams.stimulusSamplingIntervalInSeconds/simulationTimeStep;
        eyeMovementsTotalNum = round(eyeMovementsPerStimFrame*stimulusFramesNum);
        eyeMovementSequence = theMosaic.emGenSequence(eyeMovementsTotalNum);    

        if (isfield(temporalParams,'eyesDoNotMove') && (temporalParams.eyesDoNotMove))
            eyeMovementSequence = eyeMovementSequence * 0;
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
        end % stimFrameIndex
        
        % Compute photocurrent sequence
        coneIsomerizationRate = coneIsomerizationSequence/theMosaic.integrationTime;
        photocurrentSequence = theMosaic.os.compute(coneIsomerizationRate,theMosaic.pattern);
        timeAxis = (0:size(photocurrentSequence,3)-1)*theMosaic.sampleTime;
        timeAxis = timeAxis - timeAxis(end)/2;
        
        % Only include the central response
        timeIndicesToKeep = find(abs(timeAxis*1000-temporalParams.millisecondsToIncludeOffset) <= temporalParams.millisecondsToInclude/2);
        
        % Accumulate data in cell array of structs. 
        responseInstanceArray(iTrial) = struct(...
            'theMosaicIsomerizations', single(coneIsomerizationSequence(:,:,timeIndicesToKeep)), ...
             'theMosaicPhotoCurrents', single(photocurrentSequence(:,:,timeIndicesToKeep)), ...
              'theMosaicEyeMovements', eyeMovementSequence(timeIndicesToKeep,:), ...
                           'timeAxis', timeAxis(timeIndicesToKeep) ...
        );
    end % iTrial
    fprintf('Response instance array generation (%d instances) took %2.1f minutes to compute.\n', nTrials, toc/60);
    % Close progress bar
    close(progressHandle);
end

function progressHandle = generateProgressBar(initialMessage)
    progressHandle = waitbar(0, '');
    titleHandle = get(findobj(progressHandle,'Type','axes'),'Title');
    set(titleHandle,'FontSize',12, 'FontName', 'Menlo');
    waitbar(0, progressHandle, initialMessage);
    pause(0.2);
end
