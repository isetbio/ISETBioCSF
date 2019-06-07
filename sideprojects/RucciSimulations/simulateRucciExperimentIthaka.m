function simulateRucciExperiment

    [rootDir,~] = fileparts(which(mfilename));
    
    s = GetComputerInfo;
    localHostName = lower(s.networkName);
    if (contains(localHostName, 'manta'))
        dropboxRoot = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)';
        dropboxRoot = fullfile(dropboxRoot, 'IBIO_analysis', 'IBIOColorDetect', 'SideProjects');
        parforWorkers = 2;
    elseif (contains(localHostName, 'leviathan'))
        parforWorkers = 12;
        dropboxRoot = '/media/dropbox_disk/Dropbox (Aguirre-Brainard Lab)';
        dropboxRoot = fullfile(dropboxRoot, 'IBIO_analysis', 'IBIOColorDetect', 'SideProjects');
    else
        parforWorkers = 1;
        fprintf(2,'Unknown dropbox location for computer named: ''%s''.\n', s.localHostName);
        dropboxRoot ='/Volumes/SamsungT3/MATLAB/projects/IBIOColorDetect/sideprojects';
    end
    
    resourcesDir = fullfile(dropboxRoot, 'RucciSimulations');
    fprintf('Resources dir: %s\n', resourcesDir);
    
    generateScenes = ~true;
    generateOpticalImages = ~true;
    generateMosaicResponses = true;
    visualizeMosaicResponses = true;
    classifyMosaicResponses = true;
    
    if (generateMosaicResponses || visualizeMosaicResponses || classifyMosaicResponses)
        % Load the mosaic
        load(fullfile(resourcesDir, 'ConeMosaic_1.0Degs_Iterations_2000_Tolerance_0.000250.mat'), 'theMosaic');
    end
    
    % Contrast levels (exploring ....)
    contrastLevels = [1.0 0.3 0.1 0.03];
    nTrials = 512;
    fixationDurationSeconds = 0.8;
    warmupTimeSeconds = 0.4;
    mosaicIntegrationTimeSeconds = 2.5/1000;
    meanLuminanceCdPerM2 = 21*4;  % match Rucc 2007 paper, which said 21 cd/m2
    % Only compute responses for the first instance of noise stimulus    
    analyzedNoiseInstance = 1;
    
    if (~contains(localHostName, 'leviathan'))
        % ----- ONLY FOR TESTING  -----
        contrastLevels = [1 0.5];
        nTrials = 2;
        fixationDurationSeconds = 0.2;
        % ----- ONLY FOR TESTING  -----
    end
    
    
    if (generateScenes)
        noiseInstances = 2;         % only computing responses for 1 though
        stimulusSizeDegs = 1.0;     % small enough to allow faster computations
        generateAllScenes(noiseInstances, stimulusSizeDegs, meanLuminanceCdPerM2, contrastLevels, resourcesDir);  
    end
    
    if (generateOpticalImages)
        %Load previously computed scenes
        scenesFile = fullfile(resourcesDir, sprintf('scenes_luminance_%2.1f.mat', meanLuminanceCdPerM2));
        load(scenesFile, 'lowFrequencyScenes', 'highFrequencyScenes', ...
             'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho', 'contrastLevels', 'noiseInstances');
        % Display scene profiles
        displayLuminanceProfiles(lowFrequencyScenes, highFrequencyScenes, contrastLevels, noiseInstances);
        % Compute ois
        generateAllOpticalImages(lowFrequencyScenes, highFrequencyScenes, lowFrequencyScenesOrtho, highFrequencyScenesOrtho, contrastLevels, noiseInstances, resourcesDir);
    end
    
    if (generateMosaicResponses)
        % Load previously computed optical images
        oisFile = fullfile(resourcesDir, 'ois.mat');
        load(oisFile, 'lowFrequencyOIs', 'highFrequencyOIs', ...
         'lowFrequencyOIsOrtho', 'highFrequencyOIsOrtho', 'contrastLevels');
     
        nTrialsPerBlock = 1;
        
        generateAllMosaicResponses(theMosaic, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, ...
                mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, contrastLevels, analyzedNoiseInstance, ...
                nTrials, nTrialsPerBlock, parforWorkers, resourcesDir);
    end
    
    if (visualizeMosaicResponses)
        trialNo = 2;
        contrastLevel = 1.0;  
        
        figNo = 1000;
        visualizeTheResponses('highFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNo, resourcesDir, figNo);
        
        figNo = 1001;
        visualizeTheResponses('highFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNo, resourcesDir, figNo);
        
        
        figNo = 2000;
        visualizeTheResponses('lowFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNo, resourcesDir, figNo);
        
        figNo = 2001;
        visualizeTheResponses('lowFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNo, resourcesDir, figNo);
    end
    
end

function visualizeTheResponses(stimDescriptor, theMosaic, contrastLevel, theInstance, nTrials, trialNo, resourcesDir, figNo)
    
    fname = fullfile(resourcesDir, sprintf('%s_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', stimDescriptor, contrastLevel, theInstance, nTrials));
    load(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths',  'emPathsDegs', 'timeAxis')

    visualizeDynamicResponse(theMosaic, coneExcitations, photoCurrents, eyeMovementPaths, emPathsDegs, timeAxis, stimDescriptor, trialNo, figNo);    

end

function visualizeDynamicResponse(theMosaic, coneExcitations, photoCurrents, eyeMovementPaths, emPathsDegs, timeAxis, stimDescriptor, trialNo, figNo)
    
    [idxOfConesAlongHorizMeridian, idxOfConesAlongVertMeridian, ...
            eccDegsOfConesAlongHorizMeridian, ...
            eccDegsOfConesAlongVertMeridian, ...
            idxOfConesAlongHorizMeridianInSerializedList, idxOfConesAlongVertMeridianInSerializedList, ] = theMosaic.indicesForConesAlongMeridians();
        
    singleTrialConeExcitationResponse = squeeze(coneExcitations(trialNo,:,:));
    coneExcitationResponseRange = prctile(singleTrialConeExcitationResponse(:), [5 95]);
    
    singleTrialPhotocurrentResponse = squeeze(photoCurrents(trialNo,:,:));
    photocurrentResponseRange = prctile(singleTrialPhotocurrentResponse(:), [5 95]);
    
    singleTrialEyeMovement = squeeze(eyeMovementPaths(trialNo,:,:));
    singleTrialEyeMovementDegs = squeeze(emPathsDegs(trialNo,:,:));
    
    hFig = figure(figNo);
    
    for timeBin = 1:size(singleTrialConeExcitationResponse,2)
%         axHandle = subplot(1,2,1);
%         theMosaic.renderActivationMap(axHandle, squeeze(singleTrialResponse(:,timeBin)), ...
%                 'mapType', 'modulated disks', ...
%                 'signalRange', responseRange, ...
%                 'showColorBar', true, ...
%                 'labelColorBarTicks', true, ...
%                 'titleForColorBar', 'R*/cone/tau');

        % The eye movement trajectory
        subplot(2,5,1);
        plot(singleTrialEyeMovementDegs(1:timeBin,1),  singleTrialEyeMovementDegs(1:timeBin,2), 'k-'); hold on;
        plot(singleTrialEyeMovementDegs(timeBin,1)+ [-0.1 0.1],  singleTrialEyeMovementDegs(timeBin,2)*[1 1], 'b-', 'LineWidth', 1.5);
        plot(singleTrialEyeMovementDegs(timeBin,1)*[1 1],  singleTrialEyeMovementDegs(timeBin,2)+ [-0.1 0.1], 'b-', 'LineWidth', 1.5); 
        hold off;
        axis 'square';
        set(gca, 'YLim', [-0.2 0.2], 'XLim', [-0.2 0.2]);
        xlabel('horizontal position (degs)');
        ylabel('vertical position (degs)');
        
        % The cone excitations for cones along the vertical meridian
        subplot(2,5,2);
        responseVector = squeeze(singleTrialConeExcitationResponse(idxOfConesAlongVertMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongVertMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', coneExcitationResponseRange);
        xlabel('vertical position (degs)');
        ylabel('cone excitations (R*/cone/tau)');
        title('cone excitation responses (vert meridian)');
        
         % The cone excitations for cones along the horizontal meridian
        subplot(2,5,3);
        responseVector = squeeze(singleTrialConeExcitationResponse(idxOfConesAlongHorizMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongHorizMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', coneExcitationResponseRange);
        xlabel('horizontal position (degs)');
        ylabel('cone excitations (R*/cone/tau)');
        title('cone excitation responses (horiz meridian)');
        
        % XT plot of cone excitations along the horizontal meridian
        subplot(2,5,4);
        imagesc(eccDegsOfConesAlongHorizMeridian, timeAxis*1000, (singleTrialConeExcitationResponse(idxOfConesAlongHorizMeridianInSerializedList,:))');
        hold on;
         % superimpose x-eye movement trajectory
        plot(singleTrialEyeMovementDegs(:,1), timeAxis*1000, 'r-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('horizontal position (degs)');
        ylabel('time (msec)');
        colormap(gray);
        title('cone excitations (horiz meridian)');
        
        % XT plot of cone excitations along the vertical meridian
        subplot(2,5,5);
        imagesc(eccDegsOfConesAlongVertMeridian, timeAxis*1000, (singleTrialConeExcitationResponse(idxOfConesAlongVertMeridianInSerializedList,:))');
        hold on;
        % superimpose y-eye movement trajectory
        plot(-singleTrialEyeMovementDegs(:,2), timeAxis*1000, 'g-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('vertical position (degs)');
        ylabel('time (msec)');
        colormap(gray);
        title('cone excitations (vert meridian)');
        
        
        % The photocurrents for cones along the vertical meridian
        subplot(2,5,7);
        responseVector = squeeze(singleTrialPhotocurrentResponse(idxOfConesAlongVertMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongVertMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', photocurrentResponseRange);
        xlabel('vertical position (degs)');
        ylabel('photocurrent (pAmps)');
        title('photocurrent responses (vert meridian)');
        
        % The photocurrents for cones along the horizontal meridian
        subplot(2,5,8);
        responseVector = squeeze(singleTrialPhotocurrentResponse(idxOfConesAlongHorizMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongHorizMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', photocurrentResponseRange);
        xlabel('horizontal position (degs)');
        ylabel('photocurrent (pAmps)');
        title('photocurrent responses (horiz meridian)');
        
        % XT plot of photocurrents along the horizontal meridian
        subplot(2,5,9);
        imagesc(eccDegsOfConesAlongHorizMeridian, timeAxis*1000, (singleTrialPhotocurrentResponse(idxOfConesAlongHorizMeridianInSerializedList,:))');
        hold on;
        % superimpose x-eye movement trajectory
        plot(singleTrialEyeMovementDegs(:,1), timeAxis*1000, 'r-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('horizontal position (degs)');
        ylabel('time (msec)');
        title('photocurrent (horiz meridian)');
        
        % XT plot of photocurrents along the vertical meridian
        subplot(2,5,10);
        imagesc(eccDegsOfConesAlongVertMeridian, timeAxis*1000, (singleTrialPhotocurrentResponse(idxOfConesAlongVertMeridianInSerializedList,:))');
        hold on;
        % superimpose y-eye movement trajectory
        plot(-singleTrialEyeMovementDegs(:,2), timeAxis*1000, 'g-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('vertical position (degs)');
        ylabel('time (msec)');
        title('photocurrent (vert meridian)');
        
        colormap(gray);
        
        drawnow;
    end
end

function findPerformance()
    
end

function generateAllMosaicResponses(theMosaic, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, ...
                contrastLevels, analyzedNoiseInstances, nTrials, nTrialsPerBlock, parforWorkers, resourcesDir)

    % Set mosaic integration time and fixation duration
    theMosaic.integrationTime = mosaicIntegrationTimeSeconds;
    eyeMovementsNum = ceil(fixationDurationSeconds/theMosaic.integrationTime);


    % Compute fixational eye movements for desired number of trials
    [emPaths, fixEMOBJ] = theMosaic.emGenSequence(eyeMovementsNum, ...
        'nTrials', nTrials, 'centerPaths', ~true);
    timeAxis = fixEMOBJ.timeAxis;
    emPathsDegs = fixEMOBJ.emPosMicrons / theMosaic.micronsPerDegree;
    
    visualizeConeMosaicAndEMPath(theMosaic, fixEMOBJ);


    % Split in blocks to fit in memory
    nBlocks = round(nTrials/nTrialsPerBlock);

    % Compute mosaic responses for all ois
    nContrasts = numel(contrastLevels);
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:analyzedNoiseInstances

            fname = fullfile(resourcesDir, sprintf('highFrequency_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
            [coneExcitations, photoCurrents, eyeMovementPaths] = ...
                computeResponses(theMosaic, emPaths, highFrequencyOIs{theContrastLevel, theInstance}, nBlocks, nTrialsPerBlock, warmupTimeSeconds, parforWorkers);
            fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
            save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs',  'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
            
            fname = fullfile(resourcesDir, sprintf('highFrequencyOrtho_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
            [coneExcitations, photoCurrents, eyeMovementPaths] = ...
                computeResponses(theMosaic, emPaths, highFrequencyOIsOrtho{theContrastLevel, theInstance}, nBlocks, nTrialsPerBlock, warmupTimeSeconds, parforWorkers );
            fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
            save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
            
            fname = fullfile(resourcesDir, sprintf('lowFrequency_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
            [coneExcitations, photoCurrents, eyeMovementPaths] = ...
                computeResponses(theMosaic, emPaths, lowFrequencyOIs{theContrastLevel, theInstance}, nBlocks, nTrialsPerBlock, warmupTimeSeconds,parforWorkers );
            fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
            save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel','-v7.3');
            
            fname = fullfile(resourcesDir, sprintf('lowFrequencyOrtho_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', contrastLevels(theContrastLevel), theInstance, nTrials));
            [coneExcitations, photoCurrents, eyeMovementPaths] = ...
                computeResponses(theMosaic, emPaths, lowFrequencyOIsOrtho{theContrastLevel, theInstance}, nBlocks, nTrialsPerBlock, warmupTimeSeconds,parforWorkers );
            fprintf('Saving mosaic responses from %d trials to %s\n', size(coneExcitations,1), fname);
            save(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths', 'emPathsDegs', 'timeAxis', 'contrastLevels', 'theContrastLevel', '-v7.3');
        end
    end
        
end

function [theConeExcitations, thePhotoCurrents, theEyeMovementsPaths] = ...
    computeResponses(theMosaic, emPaths, theOI, nBlocks, nTrialsPerBlock, warmupTimeSeconds, parforWorkers )
            
    % Find the non-null cone indices
    nonNullConeIndices = find(theMosaic.pattern > 1);

    nTrials = nBlocks*nTrialsPerBlock;
    coneExcitations = cell(1,nTrials);
    photoCurrents = cell(1,nTrials);
    eyeMovementPaths = cell(1, nTrials);
     
    % Add warm up time points 
    nullEMPaths = zeros(size(emPaths,1),round(warmupTimeSeconds/theMosaic.integrationTime),2);
    timePointsNum = size(nullEMPaths,2);
    stimulusEMPaths = emPaths;
    emPaths = cat(2, nullEMPaths, stimulusEMPaths);
    
    
    %parfor (trialIndex = 1:nTrials, parforWorkers)
    for trialIndex = 1:nTrials
        % compute responses for this block's trials
        fprintf('Computing trial %d of %d\n', trialIndex, nTrials);
            
%         [theConeExcitations, thePhotocurrents] = ...
%             theMosaic.compute(theOI, ...
%                 'emPath', emPaths(trialIndex,:,:), ...
%                 'currentFlag', true);

        % Compute excitations only
        theConeExcitations = ...
            theMosaic.compute(theOI, ...
                'emPath', emPaths(trialIndex,:,:), ...
                'currentFlag', ~true);
            
        % Reshape to [nCones x mTimeBins]
        [theConeExcitationsXW, r, c] = RGB2XWFormat(squeeze(theConeExcitations(1,:,:,:)));
        
        % Compute photocurrents 
        thePhotocurrents = theMosaic.os.osCompute(theMosaic, ...
            'absorptionsInXWFormat', theConeExcitationsXW/coneMosaic.integrat);
        
        % Reshape to [1 x nRows x mCols x mTimeBins]
        thePhotocurrents = XW2RGBFormat(thePhotocurrents , r, c);
        thePhotocurrents = reshape(thePhotocurrents, [1 size(thePhotocurrents,1) size(thePhotocurrents,2) size(thePhotocurrents,3)]);
        
        % Only get the part of the response that is related to the stimulus
        % (with some nonCausalTimeBins)
        nonCausalTimeBins = round(0.1/theMosaic.integrationTime);
        theConeExcitations = theConeExcitations(:,:,:,timePointsNum+1-nonCausalTimeBins:end);
        thePhotocurrents = thePhotocurrents(:,:,:,timePointsNum+1-nonCausalTimeBins:end);
        
        % store to cell array
        coneExcitations{trialIndex} = reformatAllTrialsMatrix(theConeExcitations, nonNullConeIndices);
        photoCurrents{trialIndex} = reformatAllTrialsMatrix(thePhotocurrents, nonNullConeIndices);
        eyeMovementPaths{trialIndex} = emPaths(trialIndex,timePointsNum+1-nonCausalTimeBins:end,:);
    end
    
    % From cell array to array
    conesNum = numel(nonNullConeIndices);
    timeBinsNum = size(stimulusEMPaths,2);
    
    % Preallocate memory
    theConeExcitations = zeros(nTrials, conesNum, timeBinsNum);
    thePhotoCurrents = theConeExcitations;
    theEyeMovementsPaths = stimulusEMPaths*0;
    for trialIndex = 1:nTrials
        theConeExcitations(trialIndex,:,:) = coneExcitations{trialIndex}; coneExcitations{trialIndex} = [];
        thePhotoCurrents(trialIndex,:,:) = photoCurrents{trialIndex}; photoCurrents{trialIndex} = [];
        theEyeMovementsPaths(trialIndex,:,:) = eyeMovementPaths{trialIndex}; eyeMovementPaths{trialIndex} = [];
    end

end

function allTrialsMatrix = reformatAllTrialsMatrix(allTrialsMatrix, nonNullConesIndices)

    [trialsNum, coneRows, coneCols, timePointsNum] = size(allTrialsMatrix);

    allTrialsMatrix = permute(allTrialsMatrix, [2 3 1 4]);
    allTrialsMatrix = reshape(allTrialsMatrix, ...
        [coneRows * coneCols, trialsNum, timePointsNum]);

    % Only get the absorptions for the non-null cones
    allTrialsMatrix = allTrialsMatrix(nonNullConesIndices, :, :);

    % Reshape to [instances x cones x timePoints]
    allTrialsMatrix = permute(allTrialsMatrix, [2 1 3]);
end
    
function generateAllOpticalImages(lowFrequencyScenes, highFrequencyScenes, lowFrequencyScenesOrtho, highFrequencyScenesOrtho, contrastLevels, noiseInstances, resourcesDir)
    nContrasts = numel(contrastLevels);
    theOI = oiCreate('wvf human');
    
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:noiseInstances       
            lowFrequencyOIs{theContrastLevel, theInstance} = oiCompute(theOI, ...
            lowFrequencyScenes{theContrastLevel, theInstance});
        
            lowFrequencyOIsOrtho{theContrastLevel, theInstance} = oiCompute(theOI, ...
            lowFrequencyScenesOrtho{theContrastLevel, theInstance});
        
            highFrequencyOIs{theContrastLevel, theInstance} = oiCompute(theOI, ...
            highFrequencyScenes{theContrastLevel, theInstance});
        
            highFrequencyOIsOrtho{theContrastLevel, theInstance} = oiCompute(theOI, ...
            highFrequencyScenesOrtho{theContrastLevel, theInstance});
        end
    end
    
    fName = fullfile(resourcesDir, 'ois.mat');
    save(fName, 'lowFrequencyOIs', 'highFrequencyOIs', ...
         'lowFrequencyOIsOrtho', 'highFrequencyOIsOrtho', 'contrastLevels', 'noiseInstances', '-v7.3');        
end

function generateAllScenes(noiseInstances, stimulusSizeDegs, meanLuminanceCdPerM2, contrastLevels, resourcesDir)
   
    viewingDistance = 75/100;
    
    % Generate stimulus spatial modulations
    noiseNorm = nan;
    
    
    nContrasts = numel(contrastLevels);
    oriDegs = 0;
    
    
    [lowFrequencySpatialModulations, lowFrequencySpatialModulationsOrtho, spatialSupportDegs, noiseNorm] = ...
        generateStimulusSpatialModulation(stimulusSizeDegs, noiseNorm, 'low frequency', oriDegs, contrastLevels, noiseInstances);
    [highFrequencySpatialModulations, highFrequencySpatialModulationsOrtho, spatialSupportDegs, noiseNorm] = ...
        generateStimulusSpatialModulation(stimulusSizeDegs, noiseNorm, 'high frequency', oriDegs, contrastLevels, noiseInstances);
 
    % Display spatial modulations
    displaySpatialModulations(lowFrequencySpatialModulations, lowFrequencySpatialModulationsOrtho, ...
                    highFrequencySpatialModulations, highFrequencySpatialModulationsOrtho);
                
    
    % Generate ISETBio display
    presentationDisplay = generateDisplay(viewingDistance, 4*meanLuminanceCdPerM2);
    
    % Generate ISETBio scenes from spatial modulations
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:noiseInstances 
            
            lowFrequencyScenes{theContrastLevel, theInstance} = generateScene(...
                squeeze(lowFrequencySpatialModulations(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
            
            lowFrequencyScenesOrtho{theContrastLevel, theInstance} = generateScene(...
                squeeze(lowFrequencySpatialModulationsOrtho(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
             
            highFrequencyScenes{theContrastLevel, theInstance} = generateScene(...
                squeeze(highFrequencySpatialModulations(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
 
            highFrequencyScenesOrtho{theContrastLevel, theInstance} = generateScene(...
                squeeze(highFrequencySpatialModulationsOrtho(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
        end
    end
   
    fName = fullfile(resourcesDir, sprintf('scenes_luminance_%2.1f.mat', meanLuminanceCdPerM2));
    save(fName, 'lowFrequencyScenes', 'highFrequencyScenes', ...
         'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho', 'contrastLevels', 'noiseInstances', '-v7.3');
end


function realizedScene = generateScene(theStimulusSpatialModulation,presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2)
    realizedScene = sceneFromFile(theStimulusSpatialModulation, 'rgb', ...
                meanLuminanceCdPerM2, presentationDisplay);
    realizedScene = sceneSet(realizedScene, 'fov', stimulusSizeDegs);
end

function displaySpatialModulations(lowFrequencySpatialModulations, lowFrequencySpatialModulationsOrtho, ...
                    highFrequencySpatialModulations, highFrequencySpatialModulationsOrtho)

    nContrasts = size(lowFrequencySpatialModulations,2);
    
    showAllStimuli = ~true;
    if (showAllStimuli)
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', 4, ...
           'rowsNum', nContrasts, ...
           'heightMargin',  0.03, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.03, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.04);
   
        for theNoiseInstance = 1:noiseInstances
            hFig = figure(1000 + theNoiseInstance); clf;
            set(hFig, 'Position', [10 10 1000 760]);

            for iContrast = 1:nContrasts
            subplot('Position', subplotPosVectors(iContrast,1).v);
            imshow(squeeze(lowFrequencySpatialModulations(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'

            subplot('Position', subplotPosVectors(iContrast,2).v);
            imshow(squeeze(lowFrequencySpatialModulationsOrtho(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'

            subplot('Position', subplotPosVectors(iContrast,3).v);
            imshow(squeeze(highFrequencySpatialModulations(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'

            subplot('Position', subplotPosVectors(iContrast,4).v);
            imshow(squeeze(highFrequencySpatialModulationsOrtho(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'
            end
        end

        colormap(gray(1024));
    end
end

function displayLuminanceProfiles(lowFrequencyScenes, highFrequencyScenes, contrastLevels, noiseInstances)

    nContrasts = numel(contrastLevels);
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:noiseInstances 
            realizedScene = lowFrequencyScenes{theContrastLevel, theInstance};
            lumMap1(theInstance, theContrastLevel,:,:) = sceneGet(realizedScene, 'luminance');

            realizedScene = highFrequencyScenes{theContrastLevel, theInstance};
            lumMap2(theInstance, theContrastLevel,:,:) = sceneGet(realizedScene, 'luminance');
        end
    end


    fovDegs = sceneGet(realizedScene, 'horizontal fov');
    cols = sceneGet(realizedScene, 'cols');

    sampleSpacing = fovDegs / cols;
    spatialSupportDegs = (1:cols)*sampleSpacing;
    spatialSupportDegs = spatialSupportDegs - mean(spatialSupportDegs);
        
    nContrasts = size(lumMap1,2);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', nContrasts, ...
       'colsNum', 2, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.04);
   
    
    N = round(size(lumMap1,3)/2);
    CLim = [0 max([max(lumMap1(:)) max(lumMap2(:))])];
    

    figure(555);
    colormap(gray);
    for theContrastLevel = 1:nContrasts 
        subplot('Position', subplotPosVectors(theContrastLevel,1).v);
        imagesc(spatialSupportDegs, spatialSupportDegs, squeeze(lumMap1(1,theContrastLevel,:,:)));
        set(gca, 'CLim', CLim);
        axis 'image';
        
        
        subplot('Position', subplotPosVectors(theContrastLevel,2).v);
        meanSpatialProfile = mean(squeeze(lumMap1(:, theContrastLevel, N, :)),1);
        
        plot(spatialSupportDegs, meanSpatialProfile, 'rs-');
        axis 'square'
        set(gca, 'YLim', CLim);
        drawnow
    end
    
    figure(666);
    colormap(gray);
    for theContrastLevel = 1:nContrasts 
        subplot('Position', subplotPosVectors(theContrastLevel,1).v);
        imagesc(spatialSupportDegs, spatialSupportDegs, squeeze(lumMap2(1,theContrastLevel,:,:)));
        set(gca, 'CLim', CLim);
        axis 'image';
        
        subplot('Position', subplotPosVectors(theContrastLevel,2).v);
        meanSpatialProfile = mean(squeeze(lumMap2(:, theContrastLevel, N, :)),1);
        plot(spatialSupportDegs, meanSpatialProfile, 'rs-');
        axis 'square'
        set(gca, 'YLim', CLim);
        drawnow
    end
    
end

function [stimulusSpatialModulation, stimulusSpatialModulationOrtho, spatialSupportDegs, noiseNorm] = ...
    generateStimulusSpatialModulation(stimulusSizeDegs, noiseNorm, stimulusType, oriDegs, contrastLevels, noiseInstances)

    coneApertureMicrons = 2; micronsPerDegree = 300;
    stimulusPixelSizeArcMin = 0.75*coneApertureMicrons / micronsPerDegree * 60;
    stimulusWidthArcMin = stimulusSizeDegs * 60;
    
    % Grating params
    gratingParams.oriDegs = oriDegs;
    gratingParams.sigmaArcMin = stimulusSizeDegs/7*60;
    gratingParams.contrastLevels = contrastLevels;
    % Noise params
    noiseParams.steepness = 100;

    switch (stimulusType)
        case 'low frequency'
            gratingParams.sfCPD = 4;
            noiseParams.spectrumShape = 'highPassCornerFreq';
            noiseParams.cornerFrequencyCPD = 10;
        case 'high frequency'
            gratingParams.sfCPD = 11;
            noiseParams.spectrumShape = 'lowPassCornerFreq';
            noiseParams.cornerFrequencyCPD = 5;
        otherwise
            error('Unknown stimulus type: ''%s''.', stimulusType);
    end
    
    
    figNo = 1;
    [stimStruct, noiseNorm] = ...
        generateGratingInNoiseSpatialModulationPattern(...
        stimulusWidthArcMin, stimulusPixelSizeArcMin,  ...
        gratingParams, noiseParams, noiseNorm, noiseInstances, figNo);
    
    spatialSupportDegs = stimStruct.spatialSupportDegs;
    
    % Normalize to [0 .. 1]
    stimulusSpatialModulation(:,:,:,:,1) = (squeeze(stimStruct.image(:,:,:,:)) + 1)/2;
    stimulusSpatialModulation(:,:,:,:,2) = stimulusSpatialModulation(:,:,:,:,1);
    stimulusSpatialModulation(:,:,:,:,3) = stimulusSpatialModulation(:,:,:,:,1);
    
    stimulusSpatialModulationOrtho(:,:,:,:,1) = (squeeze(stimStruct.imageOrtho(:,:,:,:)) + 1)/2;
    stimulusSpatialModulationOrtho(:,:,:,:,2) = stimulusSpatialModulationOrtho(:,:,:,:,1);
    stimulusSpatialModulationOrtho(:,:,:,:,3) = stimulusSpatialModulationOrtho(:,:,:,:,1);
end


function d = generateDisplay(viewingDistance, maxLuminance)
    d = displayCreate('LCD-Apple', 'viewing distance', viewingDistance);
    d = displaySet(d, 'gTable', 'linear');
    
    scaleFactor = maxLuminance/displayGet(d,'peak luminance');
    if (scaleFactor > 1)
        spd = displayGet(d,'spd');
        spd = spd * scaleFactor;
        d = displaySet(d, 'spd', spd);
    end
    
    peakLum = displayGet(d,'peak luminance');
    
end
