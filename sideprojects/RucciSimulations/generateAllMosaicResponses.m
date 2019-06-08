function generateAllMosaicResponses(theMosaic, nullSceneOI, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, ...
                contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir)

    % Set mosaic integration time and fixation duration
    theMosaic.integrationTime = mosaicIntegrationTimeSeconds;
    eyeMovementsNum = ceil(fixationDurationSeconds/theMosaic.integrationTime);


    % Compute fixational eye movements for desired number of trials
    [emPaths, fixEMOBJ] = theMosaic.emGenSequence(eyeMovementsNum, ...
        'nTrials', nTrials, 'centerPaths', ~true);
    timeAxis = fixEMOBJ.timeAxis;
    emPathsDegs = fixEMOBJ.emPosMicrons / theMosaic.micronsPerDegree;
    emPathsMicrons = fixEMOBJ.emPosMicrons;
    clear 'fixEMOBJ';
    
    if (strcmp(eyePosition, 'stabilized'))
        emPaths = 0*emPaths;
        emPathsDegs = 0*emPathsDegs;
        emPathsMicrons = 0*emPathsMicrons;
    end
    
    visualizeConeMosaicWithEMPath(theMosaic,  emPathsMicrons);


    % Compute noise-free mosaic responses for the nullOI
    % Save original noise flags
    originalIsomerizationNoiseFlag = theMosaic.noiseFlag;
    originalPhotocurrentNoiseFlag = theMosaic.os.noiseFlag;

    % Set noiseFlags to none
    theMosaic.noiseFlag = 'none';
    theMosaic.os.noiseFlag = 'none';

    fName = fullfile(resourcesDir, sprintf('zeroContrast_nTrials_%d.mat',  nTrials));
    theContrastLevel = 0;
    [coneExcitations, photoCurrents, eyeMovementPaths] = ...
         computeResponses(theMosaic, 0*emPaths, nullSceneOI, 1, warmupTimeSeconds, parforWorkers);
    fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
    save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
    
    % Reset noiseFlags 
    theMosaic.noiseFlag = originalIsomerizationNoiseFlag;
    theMosaic.os.noiseFlag = originalPhotocurrentNoiseFlag;
    
    
    % Compute mosaic responses for all other OIs and all contrast levels
    nContrasts = numel(contrastLevels);
    theInstance = analyzedNoiseInstance;
    
    for theContrastLevel = 1:nContrasts
        stimDescriptor = 'highFrequency';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, highFrequencyOIs{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds, parforWorkers);
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');

        stimDescriptor = 'highFrequencyOrtho';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, highFrequencyOIsOrtho{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds, parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');

        stimDescriptor = 'lowFrequency';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, lowFrequencyOIs{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds,parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel','-v7.3');

        stimDescriptor = 'lowFrequencyOrtho';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths] = ...
            computeResponses(theMosaic, emPaths, lowFrequencyOIsOrtho{theContrastLevel, theInstance}, nTrials, warmupTimeSeconds,parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
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

function visualizeConeMosaicWithEMPath(theConeMosaic, emPosMicrons)

    figure(); clf;
    ax = subplot('Position', [0.1 0.05 0.35 0.95]);
    theConeMosaic.visualizeGrid(...
        'axesHandle', ax, ...
        'labelConeTypes', false, ...
        'backgroundColor', [1 1 1], ...
        'ticksInVisualDegs', true);
    set(ax, 'FontSize', 18);
    xlabel(ax,'\it space (degs)')
    ylabel(ax,'\it space (degs)');
    
    % Get the emPath in meters so we can plot it on the same scale as the cone
    % mosaic.
    emPathsMeters = emPosMicrons * 1e-6;
    hold on;
    trialNo = 1;
    plot(emPathsMeters(trialNo,:,1),emPathsMeters(trialNo,:,2), 'c-', 'LineWidth', 3.0);
    plot(emPathsMeters(trialNo,:,1),emPathsMeters(trialNo,:,2), 'b-', 'LineWidth', 1.5);

    
    ax = subplot('Position', [0.60 0.05 0.35 0.95]);
    nTimePoints = size(emPathsMeters,2);
    timeAxis = theConeMosaic.integrationTime*(1:nTimePoints);
    emPathDegs = emPathsMeters*1e6/theConeMosaic.micronsPerDegree;
    plot(timeAxis, emPathDegs(trialNo,:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(timeAxis, emPathDegs(trialNo,:,2), 'b-', 'LineWidth', 1.5);
    legend({'x-pos', 'y-pos'});
    set(gca, 'FontSize', 18);
    xlabel('\it time (seconds)');
    ylabel('\it space (degrees)');
    set(gca, 'YLim', 0.5*max(theConeMosaic.fov)*[-1 1]);
    axis 'square'
    set(gca, 'XTick', [0:0.1:1], 'YTick', [-0.3:0.1:0.3]);
    box on; grid on;
    
end

