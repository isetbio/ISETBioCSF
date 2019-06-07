function simulateRucciExperiment

    % Determine resources
    [rootDir,~] = fileparts(which(mfilename));
    [resourcesDir, parforWorkers] = determineResources(rootDir);
    
    % Actions
    smallCompute = true;
    generateScenes = true;
    generateOpticalImages = true;
    generateMosaicResponses = true;
    visualizeMosaicResponses = ~true;
    computeEnergyMechanismResponses = true;
    classifyMosaicResponses = true;
    
    % Load the cone mosaic if needed
    if (generateMosaicResponses || visualizeMosaicResponses || computeEnergyMechanismResponses)
        load(fullfile(resourcesDir, 'ConeMosaic_1.0Degs_Iterations_2000_Tolerance_0.000250.mat'), 'theMosaic');
    end
    
    % Simulation parameters
    minContrast = 1/100;
    maxContrast = 20/100;
    nContrastLevels = 10;
    contrastLevels = logspace(log10(minContrast), log10(maxContrast), nContrastLevels);
    nTrials = 512;
    fixationDurationSeconds = 0.8;
    warmupTimeSeconds = 0.4;
    mosaicIntegrationTimeSeconds = 2.5/1000;
    meanLuminanceCdPerM2 = 21;  % match Rucci 2007 paper, which said 21 cd/m2
    
    % Only compute responses for the first instance of noise stimulus    
    analyzedNoiseInstance = 1;
    
    
    if (smallCompute)
        % ----- ONLY FOR TESTING  -----
        contrastLevels = [1 0.3 0.1 0.03 0.01];
        nTrials = 4;
        warmupTimeSeconds = 0.1;
        fixationDurationSeconds = 0.3;
        % ----- ONLY FOR TESTING  -----
    end
    
    % Generate the scenes
    if (generateScenes)
        fprintf('GENERATING SCENES ...\n');
        noiseInstances = 2;         % only computing responses for 1 though
        stimulusSizeDegs = 1.0;     % small enough to allow faster computations
        generateAllScenes(noiseInstances, stimulusSizeDegs, meanLuminanceCdPerM2, contrastLevels, resourcesDir);  
    end
    
    % Generate the optical images
    if (generateOpticalImages)
        fprintf('GENERATING OPTICAL IMAGES ...\n');
        %Load previously computed scenes
        scenesFile = fullfile(resourcesDir, sprintf('scenes_luminance_%2.1f.mat', meanLuminanceCdPerM2));
        load(scenesFile, 'nullScene', 'lowFrequencyScenes', 'highFrequencyScenes', ...
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
        oisFile = fullfile(resourcesDir, sprintf('ois_luminance_%2.1f.mat', meanLuminanceCdPerM2));
        load(oisFile, 'nullSceneOI', 'lowFrequencyOIs', 'highFrequencyOIs', ...
            'lowFrequencyOIsOrtho', 'highFrequencyOIsOrtho', 'contrastLevels');
     
        visualizeOpticalImages(nullSceneOI, lowFrequencyOIs, highFrequencyOIs, lowFrequencyOIsOrtho, highFrequencyOIsOrtho, contrastLevels, analyzedNoiseInstance);
        
        generateAllMosaicResponses(theMosaic, nullSceneOI, lowFrequencyOIs, highFrequencyOIs, ...
                lowFrequencyOIsOrtho, highFrequencyOIsOrtho, ...
                mosaicIntegrationTimeSeconds, fixationDurationSeconds, warmupTimeSeconds, contrastLevels, analyzedNoiseInstance, ...
                nTrials, parforWorkers, resourcesDir);
    end
    
    % Compute the energy mechanism responses
    if (computeEnergyMechanismResponses)
        fprintf('COMPUTING ENERGY MECHANISM RESPONSES ...\n');
        % Load spatial pooling templates from the scenesFile
        scenesFile = fullfile(resourcesDir, sprintf('scenes_luminance_%2.1f.mat', meanLuminanceCdPerM2));
        
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
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir);
        
        % Compute responses of energy mechanisms to the high frequency standard orientation stimuli
        stimDescriptor = 'highFrequencyOrtho';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir);
        
        % Compute responses of energy mechanisms to the low frequency standard orientation stimuli
        stimDescriptor = 'lowFrequency';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir);
        
        % Compute responses of energy mechanisms  to the low frequency orthogonal orientation stimuli
        stimDescriptor = 'lowFrequencyOrtho';
        computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir); 
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
        %figNo = 1000;
        %visualizeAllResponses('zeroContrast', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        %figNo = 1001;
        %visualizeAllResponses('highFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        figNo = 1002;
        visualizeAllResponses('highFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        %figNo = 2001;
        %visualizeAllResponses('lowFrequency', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
        
        %figNo = 2002;
        %visualizeAllResponses('lowFrequencyOrtho', theMosaic, contrastLevel, analyzedNoiseInstance, nTrials, trialNoToVisualize, resourcesDir, figNo);
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
        parforWorkers = 12;
        dropboxRoot = '/media/dropbox_disk/Dropbox (Aguirre-Brainard Lab)';
        resourcesDir = fullfile(dropboxRoot, 'IBIO_analysis', 'IBIOColorDetect', 'SideProjects', 'RucciSimulations');
    else
        parforWorkers = 2;
        fprintf(2,'Unknown dropbox directory for computer named: ''%s''.\n', s.localHostName);
        resourcesDir = fullfile(rootDir, 'Resources');
    end
end
