function run_paper2EffectOfPulseDuration
% Compute and contrast performance at the level of photocurrents for
% different pulse durations
%
% Syntax:
%   run_paper2EffectOfPulseDuration
%
% Description:
%    Compute and contrast performance at the different pulse durations
%
%    The computation is done via the ecc-based cone efficiency & macular pigment
%    mosaic and the default Thibos subject. We use a 5 msec integration
%    time because there are no fixational eye movements. 
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
%    01/08/19  npc  Wrote it.

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
    visualizePerformance = true;
    
    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamcs of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % Assemble conditions list to be examined
    % Init condition index
    condIndex = 0;
    
    performanceSignal = 'isomerizations'; % 'photocurrents';
    performanceClassifier = 'svmV1FilterBank';
    spatialPoolingKernelType = 'V1CosUnit';            % choose between 'V1CosUnit' and 'V1QuadraturePair';
    spatialPoolingKernelActivationFunction = 'linear'; % choose between 'linear' and 'energy';
    emPathType = 'frozen0';                            % choose between 'frozen0' and 'randomNoSaccades';
    centeredEMPaths = true;                            % choose between true and 'atStimulusModulationMidPoint';
    
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = '25 ms';
    examinedCond(condIndex).frameRate = 40;
    examinedCond(condIndex).stimulusDurationInSeconds = 25/1000;
    examinedCond(condIndex).performanceClassifier = performanceClassifier;
    examinedCond(condIndex).spatialPoolingKernelParams.type = spatialPoolingKernelType;
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = spatialPoolingKernelActivationFunction;
    examinedCond(condIndex).performanceSignal = performanceSignal;
    examinedCond(condIndex).emPathType = emPathType;
    examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
    
        
    condIndex = condIndex+1;
    examinedCond(condIndex).label = '50 ms';
    examinedCond(condIndex).frameRate = 40;
    examinedCond(condIndex).stimulusDurationInSeconds = 50/1000;
    examinedCond(condIndex).performanceClassifier = performanceClassifier;
    examinedCond(condIndex).spatialPoolingKernelParams.type = spatialPoolingKernelType;
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = spatialPoolingKernelActivationFunction;
    examinedCond(condIndex).performanceSignal = performanceSignal;
    examinedCond(condIndex).emPathType = emPathType;
    examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
    
    if (~computeResponses)
        condIndex = condIndex+1;
        examinedCond(condIndex).label = '100 ms';
        examinedCond(condIndex).frameRate = 10;
        examinedCond(condIndex).stimulusDurationInSeconds = 100/1000;
        examinedCond(condIndex).performanceClassifier = performanceClassifier;
        examinedCond(condIndex).spatialPoolingKernelParams.type = spatialPoolingKernelType;
        examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = spatialPoolingKernelActivationFunction;
        examinedCond(condIndex).performanceSignal = performanceSignal;
        examinedCond(condIndex).emPathType = emPathType;
        examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
    end
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = '200 ms';
    examinedCond(condIndex).frameRate = 25;
    examinedCond(condIndex).stimulusDurationInSeconds = 200/1000;
    examinedCond(condIndex).performanceClassifier = performanceClassifier;
    examinedCond(condIndex).spatialPoolingKernelParams.type = spatialPoolingKernelType;
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = spatialPoolingKernelActivationFunction;
    examinedCond(condIndex).performanceSignal = performanceSignal;
    examinedCond(condIndex).emPathType = emPathType;
    examinedCond(condIndex).centeredEMPaths = centeredEMPaths;

    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        
        cond = examinedCond(condIndex);
        frameRate = cond.frameRate;
        stimulusDurationInSeconds = cond.stimulusDurationInSeconds;
        
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds,  frameRate, stimulusDurationInSeconds, computationInstance);
        
        % Use 1028 vs 1024 trials to differentiate results from old photocurrent computation
        params.nTrainingSamples = 1030;
        
        % Try out for subset of SFs
        params.cyclesPerDegreeExamined = [8 16 24 32 50 60]; % [16 24 32 50 60];
        
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
        params.performanceClassifier = cond.performanceClassifier;
        params.performanceSignal = cond.performanceSignal;
        params.emPathType = 'frozen0';
        
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        
        % Update params
        computePhotocurrents = true;
        
        params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance);  
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
        variedParamName = 'SignalType';
        
        theRatioLims = [0.07 1.15];
        theRatioTicks = [0.05 0.1 0.2 0.5 1.0];
        formatLabel = 'ComparedToBanksSubjects';
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'showSubjectData', ~true, ...
            'showSubjectMeanData', ~true, ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ... 
            'theLegendPosition', [0.220,0.905,0.77,0.08], ...   % custom legend position and size
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
