function run_paper2_CrowellBanksUnpublishedExperiment
% Simulate experiment 3 of Crowell & Banks
%
% Syntax:
%   run_paper2_CrowellBanksUnpublishedExperiment
%
% Description:
%    Simulate experiment 3 (spatial summation experiment) of Crowell &
%    Banks (unpublished data)
% 
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

    plotHumanData();
    
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to compute responses
    computeMosaic = ~true;
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = true;
    visualizePerformance = true;
    
    % Type of inference engine to employ
    observerTypesExamined = {'computational'};                % {'ideal', 'computational'};
    computationalObserverClassifier = 'svmV1FilterEnsemble';  % Choose from : 'svmV1FilterEnsemble', 'svmV1FilterBank'
   
    % Pupil diameter used
    pupilDiamMm = 2.5;                                  % Crowell & Banks employed a 2.5 mm artificial pupil
    
    % Mean luminance used
    luminanceCdM2 = 100;                                % Crowell & Banks employed 100 cd/m2 stimuli
    
    % Cycles/deg examined
    cyclesPerDegreeExamined = [5 14 28]; % [2.5 5 14 28];    % Crowell & Banks employed 1.75, 5, 14, and 28 c/deg.
    
    % Patch sizes examined
    patchSize2SigmaCycles = [3.3 1.7]; % [3.3 1.7 0.8 0.4];  % Gabor patch sizes (2 x sigma in degrees) employed by Crowell & Banks
    
    % Performance for threshold
    thresholdCriterionFraction = 0.75;                  % Performance threshold employed by Crowell & Banks was 75%
    
    % Response instances to compute
    nTrainingSamples = 1032;                         % 1032 to signify the Crowell & Banks runs
    
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
      
    % Assemble conditions list to be examined (different patch sizes)
    condIndex = 0;
    for patchSizeIndex = 1:numel(patchSize2SigmaCycles)
        for observerTypeIndex = 1:numel(observerTypesExamined)
            observerType = observerTypesExamined{observerTypeIndex};
            if (strcmp(observerType, 'ideal'))
                % the ideal observer: based on isomerizations with no eye movements
                observerLegend = 'ideal observer';
                emPathType = 'frozen0';
                centeredEMPaths =  'atStimulusModulationMidPoint';
                performanceSignal = 'isomerizations';
                performanceClassifier = 'mlpt';
                spatialPoolingKernelParams.type = 'V1QuadraturePair';
                spatialPoolingKernelParams.activationFunction = 'energy';
                
            elseif (strcmp(observerType, 'computational'))
                % the computational observer: based on pCurrent with fixational eye movements
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
        params.luminances = luminanceCdM2;
        params.thresholdCriterionFraction = thresholdCriterionFraction;
        params.nTrainingSamples = nTrainingSamples;
        
        % Update params based on condition currently examined
        cond = examinedCond(condIndex);
        
        % Legend
        examinedLegends{numel(examinedLegends)+1} = cond.label;
        
        % Eye movement type
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMPaths;
        
        % Signal on which inference engine is acting on
        params.performanceSignal = cond.performanceSignal;
        
        % Classifier type (mltp, svm, etc)
        params.performanceClassifier = cond.performanceClassifier;
        params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
        params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
            
        % Stimulus patch size
        params.patchSize2SigmaCycles = cond.patchSize2SigmaCycles;
        
        % Also adjust stimulus pixels based on patch size
        params.imagePixels = max([256 round(0.5*params.imagePixels * params.patchSize2SigmaCycles / 3.3)*2]);
        
        % Stimulus SF
        params.cyclesPerDegreeExamined = cond.cyclesPerDegreeExamined;
        
        % Mosaic size
        params.minimumMosaicFOVdegs = cond.minimumMosaicFOVdegs;
        
        % Keep the optical image size matched to the cone mosaic, not the
        % stimulus
        params.opticalImagePadSizeDegs = cond.minimumMosaicFOVdegs;
        
        % Adjust the 'svmV1FilterEnsemble' pooling params to match the
        % current 'patchSize2SigmaCycles'
        if (strcmp(params.performanceClassifier, 'svmV1FilterEnsemble'))
            ensembleFilterParams = struct(...
                        'spatialPositionsNum',  1, ...   % 1 results in a 3x3 grid, 2 in 5x5
                        'spatialPositionOffsetDegs', 0.025, ... 
                        'cyclesPerRFs', params.patchSize2SigmaCycles*2.5, ...           % each template contains 5 cycles of the stimulus
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
    
    if (visualizePerformance)
        condIndex = 0;
        lumIndex = 1;
        theLegends = {};
        
        data.patchSize2SigmaCycles = patchSize2SigmaCycles;
        data.observerTypesExamined = observerTypesExamined;
        data.cyclesPerDegreeExamined = cyclesPerDegreeExamined;
        
        
        for patchSizeIndex  = 1:numel(patchSize2SigmaCycles)
            for observerTypeIndex = 1:numel(observerTypesExamined)
                for sfIndex = 1:numel(cyclesPerDegreeExamined)   
                    condIndex = condIndex + 1;
                    if (sfIndex == 1) 
                        theLegends{numel(theLegends)+1} = examinedLegends{condIndex};
                    end
                    figData = theFigData{condIndex};
                    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
                    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
                    data.cpd(patchSizeIndex, observerTypeIndex, sfIndex) = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
                    data.contrastSensitivity(patchSizeIndex, observerTypeIndex, sfIndex) = 1./(thresholdContrasts*referenceContrast);
                end % sfIndex
            end % observerType
        end
        
        data.theLegends = theLegends;
        data.performanceClassifier = params.performanceClassifier;
        dataFilename = sprintf('data_%s.mat', strrep(data.performanceClassifier, ' ', '_'));
        save(dataFilename, 'data');
        
        % Generate the figure
        generateCrowellBanksFigure(data);
    end
    
end

function params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, computeMosaic, visualizeResponses, findPerformance, visualizePerformance)
                         
    % Simulation steps to perform
    params.computeMosaic = computeMosaic; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
    params.computePhotocurrentResponseInstances = computePhotocurrents && computeResponses;
    params.visualizeResponses = visualizeResponses;
    params.visualizeResponsesWithSpatialPoolingSchemeInVideo = ~true;
    params.visualizeOuterSegmentFilters = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = visualizePerformance;
    params.deleteResponseInstances = ~true;
end

function plotHumanData()

load('CrowellBanksSubjects.mat');
figure(555)
subplot(1,2,1);
plot(CrowellBanksEx3MSB_33cycles.x, CrowellBanksEx3MSB_33cycles.y, 'bs-', 'LineWidth', 1.5);
hold on;
plot(CrowellBanksEx3MSB_17cycles.x, CrowellBanksEx3MSB_17cycles.y, 'bo-','LineWidth', 1.5);
plot(CrowellBanksEx3MSB_08cycles.x, CrowellBanksEx3MSB_08cycles.y, 'bd-','LineWidth', 1.5);
plot(CrowellBanksEx3MSB_04cycles.x, CrowellBanksEx3MSB_04cycles.y, 'b^-','LineWidth', 1.5);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XLim', [1 100], 'YLim', [0.9 1000]);
legend({'3.3', '1.7', '0.8', '0.4'});
grid on
title('MSB');

subplot(1,2,2);
plot(CrowellBanksEx3JAC_33cycles.x, CrowellBanksEx3JAC_33cycles.y, 'bs-', 'LineWidth', 1.5);
hold on;
plot(CrowellBanksEx3JAC_17cycles.x, CrowellBanksEx3JAC_17cycles.y, 'bo-', 'LineWidth', 1.5);
plot(CrowellBanksEx3JAC_08cycles.x, CrowellBanksEx3JAC_08cycles.y, 'bd-', 'LineWidth', 1.5);
plot(CrowellBanksEx3JAC_04cycles.x, CrowellBanksEx3JAC_04cycles.y, 'b^-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XLim', [1 100], 'YLim', [0.9 1000]);
legend({'3.3', '1.7', '0.8', '0.4'});
title('JAC');
grid on
end

