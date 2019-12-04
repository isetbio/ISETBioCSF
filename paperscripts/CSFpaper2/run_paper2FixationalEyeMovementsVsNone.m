function run_paper2FixationalEyeMovementsVsNone
% Compute and contrast performance for fixational eye movements vs no EM.
%
% Syntax:
%   run_paper2FixationalEyeMovementsVsNone
%
% Description:
%    Compute and contrast performance for fixational eye movements vs no EM
%
%    The computation is done via the ecc-based cone efficiency & macular pigment
%    mosaic and the default Thibos subject. We use a 5 msec integration
%    time but we will examine the effect of shorter integration times. 
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
%    01/31/19  npc  Wrote it.

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Whether to compute responses
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = true;
    visualizePerformance = ~true;
    
    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
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
    
    performanceSignal = 'photocurrents'; % 'isomerizations';  % 'photocurrents', 'isomerizations';
    
    showDataFromLinearPooling = ~true;
    showDataFromQuadraturePooling = ~showDataFromLinearPooling;

    % Assemble conditions list to be examined
    % Init condition index
    condIndex = 0;
    
    % Recompute responses in which we have eye movements to fix issue with
    % mean currents
    recomputeWithFixedCurrents = ~true;
    
    
    % Reference CSF
    if (~recomputeWithFixedCurrents)
        condIndex = condIndex+1;
        examinedCond(condIndex).label = 'SVM-Template-L, noEM';
        examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
        examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
        examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
        examinedCond(condIndex).performanceSignal = performanceSignal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).centeredEMPaths = true;



        if (~computeResponses) && (showDataFromQuadraturePooling) 
            condIndex = condIndex+1;
            examinedCond(condIndex).label = 'SVM-Template-E, noEM';
            examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
            examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
            examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
            examinedCond(condIndex).performanceSignal = performanceSignal;
            examinedCond(condIndex).emPathType = 'frozen0';
            examinedCond(condIndex).centeredEMPaths = true;
        end

        if (showDataFromLinearPooling)
            condIndex = condIndex+1;
            examinedCond(condIndex).label = 'SVM-Template-L, drift';
            examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
            examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
            examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
            examinedCond(condIndex).performanceSignal = performanceSignal;
            examinedCond(condIndex).emPathType = 'randomNoSaccades';
            examinedCond(condIndex).centeredEMPaths = 'atStimulusModulationMidPoint';
        end
    %     
        if (~computeResponses) && (showDataFromQuadraturePooling)
            condIndex = condIndex+1;
            examinedCond(condIndex).label = 'SVM-Template-E, drift';
            examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
            examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
            examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
            examinedCond(condIndex).performanceSignal = performanceSignal;
            examinedCond(condIndex).emPathType = 'randomNoSaccades';
            examinedCond(condIndex).centeredEMPaths = 'atStimulusModulationMidPoint';
        end
    
    else
        computeResponses = true;
        condIndex = condIndex+1;
        examinedCond(condIndex).label = 'SVM-Template-L, drift';
        examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
        examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
        examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
        examinedCond(condIndex).performanceSignal = performanceSignal;
        examinedCond(condIndex).emPathType = 'randomNoSaccades';
        examinedCond(condIndex).centeredEMPaths = 'atStimulusModulationMidPoint';
            
    end
    
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        % Get default params
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
        
        % Use 1030 vs 1024 trials to differentiate results from when the
        % mosaic size is matched to the stimulus (see below params.minimumMosaicFOVdegs)
        params.nTrainingSamples = 1030;
        
        % Default spatial frequencies
        params.cyclesPerDegreeExamined = [32 50 60]; % [4 8 12 16 24 32 50 60];
        
        if (recomputeWithFixedCurrents)
            % High spatial frequencies only
            params.cyclesPerDegreeExamined = [32 50 60];
        end
        
        % Do not use mosaics smaller than 0.5 degs 
        %params.minimumMosaicFOVdegs = -0.328; % 0.492 IS NOT GOOD. TRY: 0.328 , 0.246, 0.158 TRY THIS TO SEE IF WE DO BETTER AT 60 C/DEG WITH ISOMERIZATIONS
        
        
        % Update params
        cond = examinedCond(condIndex);
        params.performanceClassifier = cond.performanceClassifier;
        params.performanceSignal = cond.performanceSignal;
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMPaths;
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
    
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        
        % Update params
        params = getRemainingDefaultParams(params, ...
            computePhotocurrents, computeResponses, ...
            visualizeResponses, findPerformance, visualizePerformance); 
        
        % Go !
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
        
        theRatioLims = [0.01 1.2];
        theRatioTicks = [0.01 0.05 0.1 0.2 0.5 1.0];
        if (showDataFromQuadraturePooling)
            variedParamName = 'QuadraturePoolingEyeMovements';
        else
            variedParamName = 'LinearPoolingEyeMovements';
        end
        
        formatLabel = 'ComparedToBanksSubjects';
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'showSubjectData', ~true, ...
            'showSubjectMeanData', ~true, ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ... 
            'theLegendPosition', [0.38,0.905,0.59,0.08], ...   % custom legend position and size
            'paperDir', 'CSFpaper2', ...                        % sub-directory where figure will be exported
            'figureHasFinalSize', true ...                      % publication-ready size
            );
    end
end

function params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance)
                         
    % Simulation steps to perform
    params.computeMosaic = ~true; 
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
