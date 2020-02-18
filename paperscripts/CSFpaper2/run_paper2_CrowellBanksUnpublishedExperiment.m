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

    %plotHumanData();
    
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to compute responses
    computeMosaic = ~true;
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = ~true;
    visualizePerformance = true;
    makeSummaryFigure = true;
    
    
    % Type of inference engine to employ
    observerTypeExamined = 'computational'; % 'ideal';
    computationalObserverClassifier = 'svmV1FilterEnsemble';  % Choose from : 'svmV1FilterEnsemble', 'svmV1FilterBank'
   
    % Pupil diameter used
    pupilDiamMm = 2.5;                                  % Crowell & Banks employed a 2.5 mm artificial pupil
    % Mean luminance used
    luminanceCdM2 = 100;                                % Crowell & Banks employed 100 cd/m2 stimuli
    
    % Patch sizes examined
    patchSize2SigmaCycles = 1.7;  % [3.3 1.7 0.8 0.4];  % Gabor patch sizes (2 x sigma in degrees) employed by Crowell & Banks
    
    switch (patchSize2SigmaCycles)
        case 3.3
            if (strcmp(observerTypeExamined, 'ideal'))
                cyclesPerDegExamined = [2.5 5 14 28 39 47];
            else
                cyclesPerDegExamined = [2.5 5 14 28];
            end
        case 1.7
            if (strcmp(observerTypeExamined, 'ideal'))
                cyclesPerDegExamined = [2.5 5 14 28 39 47];
            else
                cyclesPerDegExamined = [2.5 5 14 28];
            end
        case 0.8
            if (strcmp(observerTypeExamined, 'ideal'))
               cyclesPerDegExamined = [2.5 5 14 28 39 47];
            else
               cyclesPerDegExamined = [2.5 5 14 28];
            end
        case 0.4
            if (strcmp(observerTypeExamined, 'ideal'))
                cyclesPerDegExamined = [2.5 5 14 28 39 47];
            else
                cyclesPerDegExamined = [2.5 5 14];  
            end
    end
    
    
    
    % Performance for threshold
    thresholdCriterionFraction = 0.75;                  % Performance threshold employed by Crowell & Banks was 75%
    
    % Response instances to compute
    nTrainingSamples = 500;                            % 1032 to signify the Crowell & Banks runs
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamics of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % How long the stimulus is. We might be changing the duration. 100 ms is the default
    stimulusDurationInSeconds = 100/1000;
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    
      
    % Assemble conditions list to be examined (different patch sizes)
    if (strcmp(observerTypeExamined, 'ideal'))
        % the ideal observer: based on isomerizations with no eye movements
        observerLegend = 'ideal observer';
        emPathType = 'frozen0';
        centeredEMPaths =  'atStimulusModulationMidPoint';
        performanceSignal = 'isomerizations';
        performanceClassifier = 'mlpt';
        spatialPoolingKernelParams.type = 'V1QuadraturePair';
        spatialPoolingKernelParams.activationFunction = 'energy';
        % Compute photocurrent responses
        computePhotocurrents = ~true;
    elseif (strcmp(observerTypeExamined, 'computational'))
        % the computational observer: based on pCurrent with fixational eye movements
        observerLegend = computationalObserverClassifier;
        emPathType = 'randomNoSaccades';
        centeredEMPaths =  'atStimulusModulationMidPoint';
        performanceSignal = 'photocurrents';
        performanceClassifier = computationalObserverClassifier;
        spatialPoolingKernelParams.type = 'V1QuadraturePair';
        spatialPoolingKernelParams.activationFunction = 'energy';
        % Compute photocurrent responses
        computePhotocurrents = true;
    else
        error('Unknown observer type: ''%s''.', observerType);
    end

    condIndex = 0;
    for sfIndex = 1:numel(cyclesPerDegExamined)
        sf = cyclesPerDegExamined(sfIndex);
        patchSize = patchSize2SigmaCycles;

        condIndex = condIndex + 1;
        examinedCond(condIndex).label = 'computational obs.';
        examinedCond(condIndex).emPathType = emPathType;
        examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
        examinedCond(condIndex).performanceSignal = performanceSignal;
        examinedCond(condIndex).performanceClassifier = performanceClassifier;
        examinedCond(condIndex).spatialPoolingKernelParams.type = spatialPoolingKernelParams.type;
        % Compute photocurrent responses
        examinedCond(condIndex).computePhotocurrents = computePhotocurrents;

        examinedCond(condIndex).patchSize2SigmaCycles = patchSize;
        examinedCond(condIndex).cyclesPerDegreeExamined = sf;
        examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = spatialPoolingKernelParams.activationFunction;

        switch (sf)
            case 2.5
                if (examinedCond(condIndex).patchSize2SigmaCycles < 0.8)
                        examinedCond(condIndex).minimumMosaicFOVdegs = 0.98;
                    elseif (examinedCond(condIndex).patchSize2SigmaCycles < 1.7)
                        examinedCond(condIndex).minimumMosaicFOVdegs = 1.97;
                    else
                        examinedCond(condIndex).minimumMosaicFOVdegs = 3.94; % specify [] for mosaic that is matched to the stimulus
                    end

            case 5
                examinedCond(condIndex).minimumMosaicFOVdegs = 1.97; % specify [] for mosaic that is matched to the stimulus
            case 14
                examinedCond(condIndex).minimumMosaicFOVdegs = 0.98;
            case {28, 39, 47}
                examinedCond(condIndex).minimumMosaicFOVdegs = 0.492;
            otherwise
                error('minimumMosaicFOVdegs not specified');
        end
    end % sfPatchIndex
  
    
    % Go
    examinedLegends = {};

    for condIndex = 1:numel(examinedCond)
        % Get default params
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
        
        % Customize params
        params.luminancesExamined = luminanceCdM2;
        params.thresholdCriterionFraction = thresholdCriterionFraction;
        params.nTrainingSamples = nTrainingSamples;
        
        % Update params based on condition currently examined
        cond = examinedCond(condIndex);
        
        % Legend
        examinedLegends{1} = cond.label;
        
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
        params.imagePixels = max([192 round(0.5*params.imagePixels * params.patchSize2SigmaCycles / 3.3)*2]);
        
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
                        'cyclesPerRFs', params.patchSize2SigmaCycles*1.41, ...           % each template contains 5 cycles of the stimulus
                        'orientations', 0);
            fNames = fieldnames(ensembleFilterParams);
            for fNameIndex = 1:numel(fNames)
                    fName = fNames{fNameIndex};
                    params.spatialPoolingKernelParams.(fName) = ensembleFilterParams.(fName);
            end
        end
        
 
        params.parforWorkersNumForClassification = 2;
            
        % Update params
        params = getRemainingDefaultParams(params, ...
            cond.computePhotocurrents, computeResponses, computeMosaic, ...
            visualizeResponses, findPerformance, visualizePerformance); 

            
        % Go !
        [~,~, tmpFigData] = run_BanksPhotocurrentEyeMovementConditions(params);
        theFigData{1}.banksEtAlReplicate.cyclesPerDegree(condIndex) = tmpFigData.banksEtAlReplicate.cyclesPerDegree;
        theFigData{1}.banksEtAlReplicate.mlptThresholds(condIndex) = tmpFigData.banksEtAlReplicate.mlptThresholds;
    end % condIndex
    
    if (makeSummaryFigure)
        variedParamName = sprintf('CrowellBanksSpatialSummation_%1.1fcycles_', patchSize2SigmaCycles);
        
        theRatioLims = [0.07 1.15];
        theRatioTicks = [0.05 0.1 0.2 0.5 1.0];
        formatLabel = 'ComparedToBanksSubjects';
  
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'showSubjectData', ~true, ...
            'showCrowellBanksSubjectDataInsteadOfBanks87SubjectData', true, ...
            'CrowellBanksSubjectDataPatchSizeCycles', patchSize2SigmaCycles, ...
            'showSubjectMeanData', ~true, ...
            'plotFirstConditionInGray', ~true, ...
            'plotRatiosOfOtherConditionsToFirst', ~true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ... 
            'theLegendPosition', [0.40,0.88,0.58,0.08], ...    % custom legend position and size
            'paperDir', 'CSFpaper2', ...                        % sub-directory where figure will be exported
            'figureHasFinalSize', true ...                      % publication-ready size
            );
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

