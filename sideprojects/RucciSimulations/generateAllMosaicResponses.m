function generateAllMosaicResponses(theMosaic, nullSceneOI, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, ...
                contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir)

    % Set mosaic integration time and fixation duration
    theMosaic.integrationTime = mosaicIntegrationTimeSeconds;
    eyeMovementsNum = ceil((fixationDurationSeconds+warmupTimeSeconds)/theMosaic.integrationTime);


    % Compute fixational eye movements for desired number of trials
    [emPaths, fixEMOBJ] = theMosaic.emGenSequence(eyeMovementsNum, ...
        'nTrials', nTrials, 'centerPaths', ~true);
    timeAxis = fixEMOBJ.timeAxis;
     
    fullTimeAxis = timeAxis - warmupTimeSeconds;
    [~, timeBinOfStimOnset] = min(abs(fullTimeAxis));
    emPathsMicrons = fixEMOBJ.emPosMicrons;
    clear 'fixEMOBJ';
    
    % Center the emPaths to (0,0) at the stimulus onset
    emPaths = bsxfun(@minus, emPaths, emPaths(:,timeBinOfStimOnset,:));
    emPathsMicrons = bsxfun(@minus, emPathsMicrons, emPathsMicrons(:,timeBinOfStimOnset,:));
    emPathsDegs = emPathsMicrons / theMosaic.micronsPerDegree;
    
    if (strcmp(eyePosition, 'stabilized'))
        emPaths = 0*emPaths;
        emPathsDegs = 0*emPathsDegs;
        emPathsMicrons = 0*emPathsMicrons;
    end
    
    % Visualize the mosaic and a single fixational eye movement path
    visualizeConeMosaicWithEMPath(theMosaic,  emPathsMicrons, fullTimeAxis);

    
    % Compute one instance of noisy cone excitations to the nullOI
    % These are used in the computeResponses() method to replace cone
    % excitations at non-causal time bins, to simulate the transition from
    % a null scene to the stimulus scene and to also let the photocurrent
    % model warmup in the presence of the null stimulus
    [coneExcitationsToNullStimulus, ~, LMSimpulseResponses, LMScurrent]  = ...
            theMosaic.compute(nullSceneOI, ...
                'emPath', 0*squeeze(emPaths(1,:,:)), ...
                'currentFlag', true);
    nonCausalTimeBins = find(fullTimeAxis<0);
    nonCausalConeExcitationsToNullStimulus = coneExcitationsToNullStimulus(:,:,nonCausalTimeBins);

    % Save the causal emPaths only
    causalTimeBins = find(fullTimeAxis>=0);
    emPathsDegs = emPathsDegs(:,causalTimeBins,:);
    
    % Compute mosaic responses for all other OIs and all contrast levels
    nContrasts = numel(contrastLevels);
    theInstance = analyzedNoiseInstance;
    
    for theContrastLevel = 1:nContrasts
        stimDescriptor = 'highFrequency';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths, timeAxis] = ...
            computeResponses(theMosaic, nonCausalConeExcitationsToNullStimulus, LMSimpulseResponses, LMScurrent, emPaths, highFrequencyOIs{theContrastLevel, theInstance}, nTrials, fullTimeAxis, parforWorkers);
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');

        stimDescriptor = 'highFrequencyOrtho';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths, timeAxis] = ...
            computeResponses(theMosaic, nonCausalConeExcitationsToNullStimulus, LMSimpulseResponses, LMScurrent, emPaths, highFrequencyOIsOrtho{theContrastLevel, theInstance}, nTrials, fullTimeAxis, parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');

        stimDescriptor = 'lowFrequency';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths, timeAxis] = ...
            computeResponses(theMosaic, nonCausalConeExcitationsToNullStimulus, LMSimpulseResponses, LMScurrent, emPaths, lowFrequencyOIs{theContrastLevel, theInstance}, nTrials, fullTimeAxis,parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel','-v7.3');

        stimDescriptor = 'lowFrequencyOrtho';
        fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevels(theContrastLevel), theInstance, nTrials, eyePosition, resourcesDir);
        [coneExcitations, photoCurrents, eyeMovementPaths, timeAxis] = ...
            computeResponses(theMosaic, nonCausalConeExcitationsToNullStimulus, LMSimpulseResponses, LMScurrent, emPaths, lowFrequencyOIsOrtho{theContrastLevel, theInstance}, nTrials, fullTimeAxis,parforWorkers );
        fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
        save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
    end
        
    % Finally, compute noise-free mosaic responses for the nullOI. These are used to
    % compute the stimulus modulated differential responses used during
    % spatial pooling via energy mechanisms.
    
    % Save original noise flags
    originalIsomerizationNoiseFlag = theMosaic.noiseFlag;
    originalPhotocurrentNoiseFlag = theMosaic.os.noiseFlag;

    % Set the mosaic's noiseFlags to none
    theMosaic.noiseFlag = 'none';
    theMosaic.os.noiseFlag = 'none';

    % Save the responses
    fName = fullfile(resourcesDir, sprintf('zeroContrast_nTrials_%d.mat',  nTrials));
    theContrastLevel = 0;
    [coneExcitations, photoCurrents, eyeMovementPaths, timeAxis] = ...
         computeResponses(theMosaic, nonCausalConeExcitationsToNullStimulus, LMSimpulseResponses, LMScurrent, 0*emPaths, nullSceneOI, 1, fullTimeAxis, parforWorkers);
    fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fName);
    save(fName, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
    
    % Restore the mosaic's noiseFlags to their original values
    theMosaic.noiseFlag = originalIsomerizationNoiseFlag;
    theMosaic.os.noiseFlag = originalPhotocurrentNoiseFlag;
end


function [theConeExcitations, thePhotoCurrents, theEyeMovementsPaths, timeAxis] = ...
    computeResponses(theMosaic, nonCausalConeExcitationsToNullStimulus, LMSimpulseResponses, LMSmeanPhotocurrents, emPaths, theOI, nTrials, fullTimeAxis, parforWorkers)
            
    nonCausalTimeBins = find(fullTimeAxis < 0);
    causalTimeBins = find(fullTimeAxis>=0);
    stimulusTimePointsNum = numel(causalTimeBins);
    
    % Preallocate memory for the returns (only the non-causal time points, and only the non-null cones)
    nonNullConeIndices = find(theMosaic.pattern > 1);
    theConeExcitations = zeros(nTrials, numel(nonNullConeIndices), stimulusTimePointsNum);
    thePhotoCurrents = theConeExcitations;
    theEyeMovementsPaths = zeros(nTrials, stimulusTimePointsNum,2);
        
    
    parfor (trialIndex = 1:nTrials, parforWorkers)
    %for trialIndex = 1:nTrials
        fprintf('Computing trial %d of %d\n', trialIndex, nTrials);
          
        % OLD COMPUTE
%         [coneExcitations, photocurrents] = ...
%             theMosaic.compute(theOI, ...
%                 'emPath', squeeze(emPaths(trialIndex,:,:)), ...
%                 'currentFlag', true);

        % Compute excitations only, we compute photocurrents directly next
        % after we replace the non-causal cone excitations with those
        % to the null stimulus.
        coneExcitations = theMosaic.compute(theOI, ...
                'emPath', squeeze(emPaths(trialIndex,:,:)), ...
                'currentFlag', ~true);
            
        % Replace non-causal cone excitations with those to the null stimulus
        coneExcitations(:,:,nonCausalTimeBins) = nonCausalConeExcitationsToNullStimulus;
     
        % Compute photocurrents directly using the @os object
        % First, reshape to [nCones x mTimeBins]
        [theConeExcitationsXW, r, c] = RGB2XWFormat(coneExcitations);
        
        % Next, compute the photocurrents using the @os object
        thePhotocurrentsXW = theMosaic.os.osCompute(theMosaic, ...
            'absorptionsInXWFormat', theConeExcitationsXW, ...
            'interpFilters', LMSimpulseResponses, ...
            'meanCur', LMSmeanPhotocurrents);
        
        % Reshape photocurrents back to [nRows x mCols x mTimeBins]
        photocurrents = XW2RGBFormat(thePhotocurrentsXW , r, c);
        
        % Reshape full 3D hex activation maps (coneRows x coneCols x time]
        % to 2D map (non-null cones x time)
        coneExcitations = theMosaic.reshapeHex3DmapToHex2Dmap(coneExcitations);
        photocurrents = theMosaic.reshapeHex3DmapToHex2Dmap(photocurrents);
            
        % Store only causal time bins
        theConeExcitations(trialIndex,:,:) = coneExcitations(:,causalTimeBins);
        thePhotoCurrents(trialIndex,:,:) = photocurrents(:,causalTimeBins);
        theEyeMovementsPaths(trialIndex,:,:) = emPaths(trialIndex,causalTimeBins,:);
    end
    
   % Return the causal time axis
   timeAxis = fullTimeAxis(causalTimeBins);
end

function visualizeConeMosaicWithEMPath(theConeMosaic, emPosMicrons, timeAxis)
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


    subplot('Position', [0.60 0.05 0.35 0.95]);
    emPathDegs = emPathsMeters*1e6/theConeMosaic.micronsPerDegree;
    plot(timeAxis, emPathDegs(trialNo,:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(timeAxis, emPathDegs(trialNo,:,2), 'b-', 'LineWidth', 1.5);
    legend({'x-pos', 'y-pos'});
    set(gca, 'FontSize', 18);
    xlabel('\it time (seconds)');
    ylabel('\it space (degrees)');
    axis 'square'
    set(gca, 'XTick', [0:0.1:1], 'YTick', [-0.3:0.1:0.3]);
    set(gca, 'YLim', 0.5*max(theConeMosaic.fov)*[-1 1], 'XLim', [timeAxis(1) timeAxis(end)]);
    box on; grid on;
end

