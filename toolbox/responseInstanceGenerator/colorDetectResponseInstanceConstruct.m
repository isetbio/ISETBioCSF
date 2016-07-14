function responseInstance = colorDetectResponseInstanceConstruct(simulationTimeStep, gaborParams, temporalParams, oiParams, mosaicParams, theOI, theMosaic)
% responseInstance = colorDetectResponseInstanceConstruct(simulationTimeStep, gaborParams, temporalParams, mosaicParams, theOI, theMosaic)
% 
% Construct a response instance given the simulationTimeStep, gaborParams, temporalParams, mosaicParams, theOI, theMosaic
%
%
%  7/9/16  npc Wrote it.
%

    % Inform user regarding the computation progress
    progressHandle = waitbar(0,'Starting computation ...');
        
    % Save base gabor params
    theBaseGaborParams = gaborParams;
    
    % Create stimulus temporal window
    [stimulusSampleTimes, gaussianTemporalWindow, rasterModulation] = gaussianTemporalWindowCreate(temporalParams);
    if (temporalParams.addCRTrasterEffect)
        temporalParams.stimulusSamplingIntervalInSeconds = stimulusSampleTimes(2)-stimulusSampleTimes(1);
    end
    stimulusFramesNum = length(stimulusSampleTimes);
    
    % Generate eye movements for the entire stimulus duration
    eyeMovementsPerStimFrame = temporalParams.stimulusSamplingIntervalInSeconds/simulationTimeStep;
    eyeMovementsTotalNum = round(eyeMovementsPerStimFrame*stimulusFramesNum);
    eyeMovementSequence = theMosaic.emGenSequence(eyeMovementsTotalNum);
    
    if (isfield(temporalParams,'eyesDoNotMove') && (temporalParams.eyesDoNotMove))
        eyeMovementSequence = eyeMovementSequence * 0;
    end
        
    % Loop over our stimulus frames
    for stimFrameIndex = 1:stimulusFramesNum
        
        waitbar(0.9*stimFrameIndex/stimulusFramesNum, progressHandle, sprintf('Computing isomerizations for frame %d', stimFrameIndex));
        
        % Modulate stimulus contrast
        gaborParams.contrast = theBaseGaborParams.contrast * gaussianTemporalWindow(stimFrameIndex);
        
        % Apply CRT raster modulation
        if (~isempty(rasterModulation))
            gaborParams.contrast = theBaseGaborParams.contrast * gaussianTemporalWindow(stimFrameIndex) * rasterModulation(stimFrameIndex);
            gaborParams.backgroundxyY(3) = gaborParams.leakageLum + theBaseGaborParams.backgroundxyY(3) * rasterModulation(stimFrameIndex);
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
    end % for stimFrameIndex

    %% Compute photocurrent sequence
    waitbar(0.95, progressHandle, sprintf('Computing photocurrent sequence'));
    coneIsomerizationRate = coneIsomerizationSequence/theMosaic.integrationTime;
    photocurrentSequence = theMosaic.os.compute(coneIsomerizationRate,theMosaic.pattern);

    % Return a struct with the computed response and all simulation params  
    responseInstance = struct(...
        'theMosaic', theMosaic, ...
        'theMosaicIsomerizations', coneIsomerizationSequence, ...
        'theMosaicPhotoCurrents', photocurrentSequence, ...
        'theMosaicEyeMovements', eyeMovementSequence, ...
        'gaborParams', gaborParams, ...
        'temporalParams', temporalParams, ...
        'oiParams', oiParams, ...
        'mosaicParams', mosaicParams, ...
        'simulationTimeStep', simulationTimeStep...
        );
    
    % Close progress bar
    close(progressHandle);
end
