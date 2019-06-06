function simulateRucciExperiment

    [rootDir,~] = fileparts(which(mfilename));
    
    s = GetComputerInfo;
    localHostName = lower(s.networkName);
    if (contains(localHostName, 'manta'))
        dropboxRoot = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)';
        parforWorkers = 2;
    elseif (contains(localHostName, 'leviathan'))
        parforWorkers = 12;
        dropboxRoot = '/media/dropbox_disk/Dropbox (Aguirre-Brainard Lab)';
    else
        error('Unknown computer name: ''%s''.', s.localHostName);
    end
    
    resourcesDir = fullfile(dropboxRoot, 'IBIO_analysis', 'IBIOColorDetect', 'SideProjects', 'RucciSimulations');
    
    generateScenes = true;
    generateOpticalImages = true;
    generateMosaicResponses = true;
    visualizeMosaicResponses = true;
    computeEnergyMechanismResponses = true;
    classifyMosaicResponses = true;
    
    if (generateMosaicResponses || visualizeMosaicResponses || computeEnergyMechanismResponses)
        % Load the mosaic
        load(fullfile(resourcesDir, 'ConeMosaic_1.0Degs_Iterations_2000_Tolerance_0.000250.mat'), 'theMosaic');
    end
    
    % Contrast levels (exploring ....)
    contrastLevels = [1.0 0.3 0.1 0.03];
    nTrials = 512;
    fixationDurationSeconds = 0.8;
    warmupTimeSeconds = 0.4;
    mosaicIntegrationTimeSeconds = 2.5/1000;
    meanLuminanceCdPerM2 = 21;  % match Rucci 2007 paper, which said 21 cd/m2
    % Only compute responses for the first instance of noise stimulus    
    analyzedNoiseInstance = 1;
    
    if (contains(localHostName, 'manta'))
        % ----- ONLY FOR TESTING  -----
        contrastLevels = [1 0.5];
        nTrials = 4;
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
        load(scenesFile, 'nullScene', 'lowFrequencyScenes', 'highFrequencyScenes', ...
             'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho', 'contrastLevels', 'noiseInstances');
        % Display scene profiles
        visualizeSceneLuminanceProfiles(nullScene, lowFrequencyScenes, highFrequencyScenes, contrastLevels, noiseInstances);
        % Compute ois
        generateAllOpticalImages(nullScene, lowFrequencyScenes, highFrequencyScenes, lowFrequencyScenesOrtho, highFrequencyScenesOrtho, contrastLevels, noiseInstances, meanLuminanceCdPerM2, resourcesDir);
    end
    
    if (generateMosaicResponses)
        % Load previously computed optical images
        oisFile = fullfile(resourcesDir, sprintf('ois_luminance_%2.1f.mat', meanLuminanceCdPerM2));
        load(oisFile, 'nullSceneOI', 'lowFrequencyOIs', 'highFrequencyOIs', ...
         'lowFrequencyOIsOrtho', 'highFrequencyOIsOrtho', 'contrastLevels');
     
        visualizeOpticalImages(nullSceneOI, lowFrequencyOIs, highFrequencyOIs, lowFrequencyOIsOrtho, highFrequencyOIsOrtho, contrastLevels, analyzedNoiseInstance);
        
        generateAllMosaicResponses(theMosaic, nullSceneOI, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, ...
                mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, contrastLevels, analyzedNoiseInstance, ...
                nTrials, parforWorkers, resourcesDir);
    end
    
    if (computeEnergyMechanismResponses)
        % Load spatial pooling templates from the scenesFile
        scenesFile = fullfile(resourcesDir, sprintf('scenes_luminance_%2.1f.mat', meanLuminanceCdPerM2));
        load(scenesFile, 'lowFrequencyScenes', 'highFrequencyScenes', ...
            'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho', ...
            'lowFrequencyTemplate', 'lowFrequencyTemplateOrtho', ...
            'highFrequencyTemplate', 'highFrequencyTemplateOrtho', ...
            'spatialSupportDegs', 'contrastLevels');
        
        % Generate spatial pooling kernels
        [~,maxContrastIndex] = max(contrastLevels);
        noiseInstance = 1;
        spatialPoolingKernels = generateSpatialPoolingKernels(theMosaic, ...
            lowFrequencyTemplate, lowFrequencyTemplateOrtho, ...
            highFrequencyTemplate, highFrequencyTemplateOrtho, ...
            lowFrequencyScenes{maxContrastIndex,noiseInstance}, ...
            lowFrequencyScenesOrtho{maxContrastIndex,noiseInstance}, ...
            highFrequencyScenes{maxContrastIndex,noiseInstance}, ...
            highFrequencyScenesOrtho{maxContrastIndex,noiseInstance}, ...
            spatialSupportDegs); 
        
        
        % Compute outputs of energy mechanisms (tuned to standard and orthogonal orientation) to the high frequency standard orientation stimuli
        %stimDescriptor = 'highFrequency';
        %computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir);
        
        % Compute outputs of energy mechanisms (tuned to standard and orthogonal orientation) to the high frequency orthogonal orientation stimuli
        stimDescriptor = 'highFrequencyOrtho';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir);
        
        
        % Compute outputs of energy mechanisms (tuned to standard and orthogonal orientation) to the low frequency standard orientation stimuli
        %stimDescriptor = 'lowFrequency';
        %computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir);
        
        % Compute outputs of energy mechanisms (tuned to standard and orthogonal orientation) to the low frequency orthogonal orientation stimuli
       % stimDescriptor = 'lowFrequencyOrtho';
       % computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir); 
    end
    
    
    if (classifyMosaicResponses) 
        stimDescriptor = 'highFrequency'; figNo = 3000;
        visualizeEnergyResponses(stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, resourcesDir, figNo);
        
        stimDescriptor = 'highFrequencyOrtho'; figNo = figNo + 1;
        visualizeEnergyResponses(stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, resourcesDir, figNo);
        
        stimDescriptor = 'lowFrequency'; figNo = figNo + 1;
        visualizeEnergyResponses(stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, resourcesDir, figNo);
        
        stimDescriptor = 'lowFrequencyOrtho'; figNo = figNo + 1;
        visualizeEnergyResponses(stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, resourcesDir, figNo);
        
    end
    
    
    if (visualizeMosaicResponses)
        trialNoToVisualize = 1;
        contrastLevel = 1.0;  
        
        figNo = 1000;
        %visualizeAllResponses('zeroContrast', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        figNo = 1001;
        visualizeAllResponses('highFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        figNo = 1002;
        visualizeAllResponses('highFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        figNo = 2001;
        %visualizeAllResponses('lowFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        figNo = 2002;
        %visualizeAllResponses('lowFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
    end 
end