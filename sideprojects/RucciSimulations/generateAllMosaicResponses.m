function generateAllMosaicResponses(theMosaic, nullSceneOI, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, ...
                contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir)

    % Set mosaic integration time and fixation duration
    theMosaic.integrationTime = mosaicIntegrationTimeSeconds;
    eyeMovementsNum = ceil(fixationDurationSeconds/theMosaic.integrationTime);


    % Compute fixational eye movements for desired number of trials
    [emPaths, fixEMOBJ] = theMosaic.emGenSequence(eyeMovementsNum, ...
        'nTrials', nTrials, 'centerPaths', ~true);
    timeAxis = fixEMOBJ.timeAxis;
    emPathsDegs = fixEMOBJ.emPosMicrons / theMosaic.micronsPerDegree;
    
    visualizeConeMosaicAndEMPath(theMosaic, fixEMOBJ);


    % Compute noise-free mosaic responses for the nullOI
    % Save original noise flags
    originalIsomerizationNoiseFlag = theMosaic.noiseFlag;
    originalPhotocurrentNoiseFlag = theMosaic.os.noiseFlag;

    % Set noiseFlags to none
    theMosaic.noiseFlag = 'none';
    theMosaic.os.noiseFlag = 'none';

    fname = fullfile(resourcesDir, sprintf('zeroContrast_nTrials_%d.mat',  nTrials));
    theContrastLevel = 0;
    [coneExcitations, photoCurrents, eyeMovementPaths] = ...
         computeResponses(theMosaic, 0*emPaths, nullSceneOI, 1, warmupTimeSeconds, parforWorkers);
    fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
    save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
    
    % Reset noiseFlags 
    theMosaic.noiseFlag = originalIsomerizationNoiseFlag;
    theMosaic.os.noiseFlag = originalPhotocurrentNoiseFlag;
    
    
    % Compute mosaic responses for all other OIs and all contrast levels
    nContrasts = numel(contrastLevels);
    theInstance = analyzedNoiseInstance;
    
    for theContrastLevel = 1:nContrasts
        fname = fullfile(resourcesDir, sprintf('highFrequency_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, highFrequencyOIs{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds, parforWorkers);
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
        save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');

        fname = fullfile(resourcesDir, sprintf('highFrequencyOrtho_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, highFrequencyOIsOrtho{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds, parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
        save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');

        fname = fullfile(resourcesDir, sprintf('lowFrequency_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, lowFrequencyOIs{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds,parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
        save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel','-v7.3');

        fname = fullfile(resourcesDir, sprintf('lowFrequencyOrtho_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, lowFrequencyOIsOrtho{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds,parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
        save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
    end
        
end


function [theConeExcitations, thePhotoCurrents, theEyeMovementsPaths] = computeResponses(theMosaic, emPaths, theOI, nTrials, warmupTimeSeconds, parforWorkers )
            
    % Find the non-null cone indices
    nonNullConeIndices = find(theMosaic.pattern > 1);
    allConesNum = numel(theMosaic.pattern);
    
    % Add 0.5 secosome time points 
    nullEMPaths = zeros(size(emPaths,1),round(warmupTimeSeconds/theMosaic.integrationTime),2);
    timePointsNum = size(nullEMPaths,2);
    
    % Concatenate null with stimulus
    stimulusEMPaths = emPaths;
    stimulusTimePointsNum = size(emPaths,2);
    
    emPaths = cat(2, nullEMPaths, stimulusEMPaths);
    allTimePointsNum = size(emPaths,2);
    
    % Preallocate memory for the returns
    theConeExcitations = zeros(nTrials, numel(nonNullConeIndices), stimulusTimePointsNum);
    thePhotoCurrents = theConeExcitations;
    theEyeMovementsPaths = zeros(nTrials, stimulusTimePointsNum,2);
    
    parfor (trialIndex = 1:nTrials, parforWorkers)
    %for trialIndex = 1:nTrials
        fprintf('Computing trial %d of %d\n', trialIndex, nTrials);
            
        [coneExcitations, photocurrents] = ...
            theMosaic.compute(theOI, ...
                'emPath', squeeze(emPaths(trialIndex,:,:)), ...
                'currentFlag', true);
        
        % Reshape full 3D hex activation map (coneRows x coneCols x time]
        % to 2D map (non-null cones x time)
        coneExcitations = theMosaic.reshapeHex3DmapToHex2Dmap(coneExcitations);
        photocurrents = theMosaic.reshapeHex3DmapToHex2Dmap(photocurrents);
        
        % store only causal time bins
        theConeExcitations(trialIndex,:,:) = coneExcitations(:,timePointsNum+(1:stimulusTimePointsNum));
        thePhotoCurrents(trialIndex,:,:) = photocurrents(:,timePointsNum+(1:stimulusTimePointsNum));
        theEyeMovementsPaths(trialIndex,:,:) = emPaths(trialIndex,timePointsNum+(1:stimulusTimePointsNum),:);
    end
end

