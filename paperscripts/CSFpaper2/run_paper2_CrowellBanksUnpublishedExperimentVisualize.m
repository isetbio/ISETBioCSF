function run_paper2_CrowellBanksUnpublishedExperimentVisualize
% Visualize responses from experiment 3 of Crowell & Banks
%exit
% Syntax:
%   run_paper2_CrowellBanksUnpublishedExperiment
%
% Description:
%    Compute simulated responses to experiment 3 (spatial summation experiment) 
%    of Crowell & Banks (unpublished data) only for 100% contrast and
%    visualize mosaic responses
% 
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History
%    09/20/19  npc  Wrote it.

    
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to compute responses
    computeMosaic = ~true;
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = true;
    visualizePerformance = ~true;
    
    % Pupil diameter used
    pupilDiamMm = 2.5;                                  % Crowell & Banks employed a 2.5 mm artificial pupil
    
    % Mean luminance used
    luminanceCdM2 = 100;                                % Crowell & Banks employed 100 cd/m2 stimuli
    
    % Cycles/deg examined
    cyclesPerDegreeExamined = [14]; %[5 14 28]; % [2.5 5 14 28];    % Crowell & Banks employed 1.75, 5, 14, and 28 c/deg.
    
    % Patch sizes examined
    patchSize2SigmaCycles = [1.7 0.8 0.4]; % [3.3 1.7 0.8 0.4];  % Gabor patch sizes (2 x sigma in degrees) employed by Crowell & Banks
    
    % Performance for threshold
    thresholdCriterionFraction = 0.75;                  % Performance threshold employed by Crowell & Banks was 75%
    
    % Response instances to compute
    nTrainingSamples = 1032-1028;                         % 1032 to signify the Crowell & Banks runs
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamics of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % How long the stimulus is. We might be changing the duration. 100 ms is the default
    stimulusDurationInSeconds = 100/1000;
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    % Compute photocurrent responses
    computePhotocurrents = true;
      
    observerTypesExamined = {'computational'};                % {'ideal', 'computational'};
    computationalObserverClassifier = 'svmV1FilterBank';  % Choose from : 'svmV1FilterEnsemble', 'svmV1FilterBank'
    
    % Assemble conditions list to be examined (different patch sizes)
    condIndex = 0;
    for patchSizeIndex = 1:numel(patchSize2SigmaCycles)
        for observerTypeIndex = 1:numel(observerTypesExamined)
            observerType = observerTypesExamined{observerTypeIndex};
            if (strcmp(observerType, 'ideal'))
                % the ideal observer, isomerizations, no eye movements
                observerLegend = 'ideal observer';
                emPathType = 'frozen0';
                centeredEMPaths =  'atStimulusModulationMidPoint';
                performanceSignal = 'isomerizations';
                performanceClassifier = 'mlpt';
                spatialPoolingKernelParams.type = 'V1QuadraturePair';
                spatialPoolingKernelParams.activationFunction = 'energy';
            elseif (strcmp(observerType, 'computational'))
                % the computational observer, pCurrent + fix. eye movements
                observerLegend = computationalObserverClassifier;
                emPathType = 'randomNoSaccades';
                centeredEMPaths =  'atStimulusModulationMidPoint';
                performanceSignal = 'photocurrents';
                performanceClassifier = computationalObserverClassifier;
                spatialPoolingKernelParams.type = 'V1QuadraturePair';
                spatialPoolingKernelParams.activationFunction = 'energy';
            else
                error('Unknown observer type: ''%s''.', observerType);
            end
            
            for sfIndex = 1:numel(cyclesPerDegreeExamined)   
                condIndex = condIndex + 1;
                examinedCond(condIndex).label = sprintf('%s, %2.1f cycles', observerLegend, patchSize2SigmaCycles(patchSizeIndex));
                examinedCond(condIndex).emPathType = emPathType;
                examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
                examinedCond(condIndex).performanceSignal = performanceSignal;
                examinedCond(condIndex).performanceClassifier = performanceClassifier;
                examinedCond(condIndex).spatialPoolingKernelParams.type = spatialPoolingKernelParams.type;
                
                examinedCond(condIndex).patchSize2SigmaCycles = patchSize2SigmaCycles(patchSizeIndex);
                examinedCond(condIndex).cyclesPerDegreeExamined = cyclesPerDegreeExamined(sfIndex);
                examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = spatialPoolingKernelParams.activationFunction;
                
                switch (cyclesPerDegreeExamined(sfIndex))
                    case 2.5
                        examinedCond(condIndex).minimumMosaicFOVdegs = 3.94; % specify [] for mosaic that is matched to the stimulus
                    case 5
                        examinedCond(condIndex).minimumMosaicFOVdegs = 1.97; % specify [] for mosaic that is matched to the stimulus
                    case 14
                        examinedCond(condIndex).minimumMosaicFOVdegs = 0.98; 
                    case 28
                        examinedCond(condIndex).minimumMosaicFOVdegs = 0.492; 
                    otherwise
                        error('minimumMosaicFOVdegs not specified');
                end
            end % sfIndex
        end % observerTypeIndex
    end % patchSizeIndex
    
    % Go
    examinedLegends = {};

    for condIndex = 1:numel(examinedCond)
        % Get default params
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
        
        % Customize params
        params.luminancesExamined = luminanceCdM2;
        params.thresholdCriterionFraction = thresholdCriterionFraction;
        params.nTrainingSamples = nTrainingSamples;
        
        % Update params
        cond = examinedCond(condIndex);
        
        % Legend
        examinedLegends{numel(examinedLegends)+1} = cond.label;
        
        % Eye movement typw
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMPaths;
        
        % Signal
        params.performanceSignal = cond.performanceSignal;
        
        % Classifier
        params.performanceClassifier = cond.performanceClassifier;
        params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
        params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
            
        % Stimulus patch size
        params.patchSize2SigmaCycles = cond.patchSize2SigmaCycles;
        
        % Also adjust stimulus pixels
        params.imagePixels = max([256 round(0.5*params.imagePixels * params.patchSize2SigmaCycles / 3.3)*2]);

        % Stimulus SF
        params.cyclesPerDegreeExamined = cond.cyclesPerDegreeExamined;
        
        % Mosaic size
        params.minimumMosaicFOVdegs = cond.minimumMosaicFOVdegs;
        
        % Keep the optical image size matched to the cone mosaic, not the
        % stimulus
        params.opticalImagePadSizeDegs = cond.minimumMosaicFOVdegs;
        
         % Only the max contrast level
        params.lowContrast = 1.0; 
        params.highContrast =  1.0; 
        params.nContrastsPerDirection = 1;
        
        % Adjust the 'svmV1FilterEnsemble' pooling params to match the
        % current 'patchSize2SigmaCycles'
        if (strcmp(params.performanceClassifier, 'svmV1FilterEnsemble'))
            ensembleFilterParams = struct(...
                        'spatialPositionsNum',  1, ...   % 1 results in a 3x3 grid, 2 in 5x5
                        'spatialPositionOffsetDegs', 0.025, ... 
                        'cyclesPerRFs', params.patchSize2SigmaCycles*1.41, ...           % each template contains 5 cycles of the stimulus
                        'orientations', 0);
            fNames = fieldnames(ensembleFilterParams);
            for fNameIndex = 1:numel(fNames)
                    fName = fNames{fNameIndex};
                    params.spatialPoolingKernelParams.(fName) = ensembleFilterParams.(fName);
            end
        end
        
        % Update params
        params = getRemainingDefaultParams(params, ...
            computePhotocurrents, computeResponses, computeMosaic, ...
            visualizeResponses, findPerformance, visualizePerformance); 

            
        % Go !
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
end



function params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, computeMosaic, visualizeResponses, findPerformance, visualizePerformance)
                         
    % Simulation steps to perform
    params.computeMosaic = computeMosaic; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
    params.computePhotocurrentResponseInstances = computePhotocurrents && computeResponses;
    params.visualizeResponses = visualizeResponses;
    params.visualizeOuterSegmentFilters = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = true;
    params.visualizeSpatialPoolingScheme = true;
    params.visualizeResponsesWithSpatialPoolingSchemeInVideo = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = visualizePerformance;
    params.deleteResponseInstances = ~true;
end

