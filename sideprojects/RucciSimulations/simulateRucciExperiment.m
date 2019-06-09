function simulateRucciExperiment

    % Determine resources
    [rootDir,~] = fileparts(which(mfilename));
    [resourcesDir, parforWorkers, localHostName] = determineResources(rootDir);
    
    % if we are not on Leviathan do a small simulation
    smallCompute = ~contains(localHostName, 'leviathan');
    
    % Actions
    analyzeStimulusSpatioTemporalSpectra = ~true;
    generateScenes = true;
    generateOpticalImages = true;
    generateMosaicResponses = true;
    visualizeMosaicResponses = ~true;
    computeEnergyMechanismResponses = true;
    estimatePerformance = true;
    
    % Load the cone mosaic if needed
    if (generateMosaicResponses || visualizeMosaicResponses || computeEnergyMechanismResponses)
        load(fullfile(resourcesDir, 'ConeMosaic_1.0Degs_Iterations_2000_Tolerance_0.000250.mat'), 'theMosaic');
    end
    
    % Simulation parameters
    minContrast = 0.03/100;
    maxContrast = 3/100;
    nContrastLevels = 6;
    contrastLevels = logspace(log10(minContrast), log10(maxContrast), nContrastLevels);
    stimulusSizeDegs = 1.0;
    nTrials = 512;
    eyePosition = 'dynamic';
    %eyePosition = 'stabilized';
    fixationDurationSeconds = 0.8;
    warmupTimeSeconds = 0.5;
    mosaicIntegrationTimeSeconds = 2.5/1000;
    meanLuminanceCdPerM2 = 21;  % match Rucci 2007 paper, which said 21 cd/m2
    
    % Only compute responses for the first instance of noise stimulus    
    analyzedNoiseInstance = 1;
    
    
    if (smallCompute)
        % ----- ONLY FOR TESTING  -----
        nContrastLevels = 2;
        minContrast = 99/100;
        maxContrast = 100/100;
        contrastLevels = logspace(log10(minContrast), log10(maxContrast), nContrastLevels);
        nTrials = 2;
        warmupTimeSeconds = 0.1;
        fixationDurationSeconds = 0.2;
        % ----- ONLY FOR TESTING  -----
    end
    
    if (analyzeStimulusSpatioTemporalSpectra) 
        noiseInstances = 3;
        reComputeSpectralAnalyses = true;
        analyzeStabilizedAndDynamicSpectra(stimulusSizeDegs, fixationDurationSeconds, noiseInstances, reComputeSpectralAnalyses);
        return;
    end
    
    % Generate the scenes
    if (generateScenes)
        fprintf('GENERATING SCENES ...\n');
        noiseInstances = 2;         % only computing responses for 1 though
        generateAllScenes(noiseInstances, stimulusSizeDegs, meanLuminanceCdPerM2, contrastLevels, resourcesDir);  
    end
    
    % Generate the optical images
    if (generateOpticalImages)
        fprintf('GENERATING OPTICAL IMAGES ...\n');
        %Load previously computed scenes
        scenesFileName = scenesDataFileName(meanLuminanceCdPerM2, contrastLevels, resourcesDir);
        load(scenesFileName, 'nullScene', 'lowFrequencyScenes', 'highFrequencyScenes', ...
             'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho', 'contrastLevels', 'noiseInstances');
        % Display scene profiles
        visualizeSceneLuminanceProfiles(nullScene, highFrequencyScenes, highFrequencyScenesOrtho, contrastLevels, noiseInstances);
        % Compute ois
        generateAllOpticalImages(nullScene, lowFrequencyScenes, highFrequencyScenes, lowFrequencyScenesOrtho, highFrequencyScenesOrtho, contrastLevels, noiseInstances, meanLuminanceCdPerM2, resourcesDir);
    end
    
    % Generate the cone mosaic responses
    if (generateMosaicResponses)
        fprintf('GENERATING CONE MOSAIC RESPONSES ...\n');
        % Load previously computed optical images
        oisFileName = oisDataFileName(meanLuminanceCdPerM2, contrastLevels, resourcesDir);
        load(oisFileName, 'nullSceneOI', 'lowFrequencyOIs', 'highFrequencyOIs', ...
            'lowFrequencyOIsOrtho', 'highFrequencyOIsOrtho', 'contrastLevels');
     
        visualizeOpticalImages(nullSceneOI, lowFrequencyOIs, highFrequencyOIs, lowFrequencyOIsOrtho, highFrequencyOIsOrtho, contrastLevels, analyzedNoiseInstance);
        
        generateAllMosaicResponses(theMosaic, nullSceneOI, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, ...
                mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, contrastLevels, analyzedNoiseInstance, ...
                nTrials, eyePosition, parforWorkers, resourcesDir);
    end
    
    % Compute the energy mechanism responses
    if (computeEnergyMechanismResponses)
        fprintf('COMPUTING ENERGY MECHANISM RESPONSES ...\n');
        % Load spatial pooling templates from the scenesFile
        scenesFile = scenesDataFileName(meanLuminanceCdPerM2, contrastLevels, resourcesDir);
        
        % Load the stimulus-derived templates, which are used in conjuction
        % with the coneMosaic to derive the spatial pooling weights
        load(scenesFile, 'lowFrequencyTemplate', 'lowFrequencyTemplateOrtho', ...
            'highFrequencyTemplate', 'highFrequencyTemplateOrtho', ...
            'spatialSupportDegs', 'contrastLevels');
        
        % Generate spatial pooling kernels
        spatialPoolingKernels = generateSpatialPoolingKernels(theMosaic, ...
            lowFrequencyTemplate, lowFrequencyTemplateOrtho, ...
            highFrequencyTemplate, highFrequencyTemplateOrtho, spatialSupportDegs); 
        
        % Visualize pooling kernels together with the stimuli
        visualizePoolingKernelsAndStimuli = ~true;
        if (visualizePoolingKernelsAndStimuli) 
            % Load the scenes
            load(scenesFile, 'lowFrequencyScenes', 'highFrequencyScenes', ...
                'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho');
        
            [~,maxContrastIndex] = max(contrastLevels); noiseInstance = 1;
            visualizeSpatialPoolingKernelsAndStimuli(spatialPoolingKernels, ...
                lowFrequencyScenes{maxContrastIndex,noiseInstance}, ...
                lowFrequencyScenesOrtho{maxContrastIndex,noiseInstance}, ...
                highFrequencyScenes{maxContrastIndex,noiseInstance}, ...
                highFrequencyScenesOrtho{maxContrastIndex,noiseInstance}, ...
                spatialSupportDegs);
        end
        
        % Compute responses of energy mechanisms  to the high frequency orthogonal orientation stimuli
        stimDescriptor = 'highFrequency';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir);
        
        % Compute responses of energy mechanisms to the high frequency standard orientation stimuli
        stimDescriptor = 'highFrequencyOrtho';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir);
        
        % Compute responses of energy mechanisms to the low frequency standard orientation stimuli
        stimDescriptor = 'lowFrequency';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir);
        
        % Compute responses of energy mechanisms  to the low frequency orthogonal orientation stimuli
        stimDescriptor = 'lowFrequencyOrtho';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir); 
     end
    
    
    if (estimatePerformance) 
        % Estimate performance for the high frequency stimulus
        stimDescriptor = 'highFrequency';
        figNo = 3000;
        estimatePerformanceForStimulus(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, contrastLevels, resourcesDir,  figNo);
        
        % Estimate performance for the low frequency stimulus
        stimDescriptor = 'lowFrequency';
        figNo = 4000;
        estimatePerformanceForStimulus(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, contrastLevels, resourcesDir,  figNo);
    end
    
    
    if (visualizeMosaicResponses)
        trialNoToVisualize = 1; 
        contrastLevel = max(contrastLevels);  
        %figNo = 1000;
        %visualizeAllResponses('zeroContrast', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, eyePosition, trialNoToVisualize, resourcesDir, figNo);
        
        %figNo = 1001;
        %visualizeAllResponses('highFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, eyePosition, trialNoToVisualize, resourcesDir, figNo);
        
        figNo = 1002;
        visualizeAllResponses('highFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, eyePosition, trialNoToVisualize, resourcesDir, figNo);
        
        %figNo = 2001;
        %visualizeAllResponses('lowFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, eyePosition, trialNoToVisualize, resourcesDir, figNo);
        
        %figNo = 2002;
        %visualizeAllResponses('lowFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, eyePosition, trialNoToVisualize, resourcesDir, figNo);
    end 
end

function [resourcesDir, parforWorkers, localHostName] = determineResources(rootDir)
    s = GetComputerInfo;
    localHostName = lower(s.networkName);
    if (contains(localHostName, 'manta'))
        dropboxRoot = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)';
        resourcesDir = fullfile(dropboxRoot, 'IBIO_analysis', 'IBIOColorDetect', 'SideProjects', 'RucciSimulations');
        parforWorkers = 4;
    elseif (contains(localHostName, 'leviathan'))
        parforWorkers = 14;
        dropboxRoot = '/media/dropbox_disk/Dropbox (Aguirre-Brainard Lab)';
        resourcesDir = fullfile(dropboxRoot, 'IBIO_analysis', 'IBIOColorDetect', 'SideProjects', 'RucciSimulations');
    else
        parforWorkers = 2;
        fprintf(2,'Unknown dropbox directory for computer named: ''%s''.\n', s.localHostName);
        resourcesDir = fullfile(rootDir, 'Resources');
    end
end
