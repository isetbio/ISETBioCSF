function run_adaptiveOpticsCSFunderSpatialSummation
% Compute and contrast performance at the level of isomerizations 
% under AO conditions with and without spatial pooling
%
% Syntax:
%   run_adaptiveOpticsCSFunderSpatialSummation
%
% Description:
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
%    10/10/19  npc  Wrote it.

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = ~true;
    
    % Whether to compute responses
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = true;
    visualizePerformance = ~true;
    
    emPathType = 'frozen0';
    centeredEMPaths = false;
    performanceSignal = 'isomerizations';
   
    computePhotocurrents = false;
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamcs of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % How long the stimulus is.
    % We might be changing the duration. 100 ms is the default
    stimulusDurationInSeconds = 100/1000;
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    % Assemble conditions list to be examined
    % Init condition index
    condIndex = 0;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'Human-3mm, SVM - No Spatial Summation';
    examinedCond(condIndex).performanceClassifier = 'svm';
    examinedCond(condIndex).pupilDiamMm = 3.0;
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';

    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'AO-8mm, SVM - No Spatial Summation';
    examinedCond(condIndex).performanceClassifier = 'svm';
    examinedCond(condIndex).pupilDiamMm = 8.0;
    examinedCond(condIndex).opticsModel = 'AOoptics80mmPupil';
    
%     condIndex = condIndex+1;
%     examinedCond(condIndex).label = 'AO-8mm, SVM - 1.7 ArcMin Summation';
%     examinedCond(condIndex).performanceClassifier = 'svmGaussPooledResponses';
%     examinedCond(condIndex).pupilDiamMm = 8.0;
%     examinedCond(condIndex).opticsModel = 'AOoptics80mmPupil';
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        
        cond = examinedCond(condIndex);
        params = getCSFPaper2DefaultParams(cond.pupilDiamMm, integrationTimeMilliseconds,  frameRate, stimulusDurationInSeconds, computationInstance);
        
        % Use 1024 trials
        params.nTrainingSamples = 1024;
        
        % Try out for subset of SFs
        params.cyclesPerDegreeExamined = [12 24 32 60]; % [4 8 12 16 24 32 50 60];
        
        
        params.performanceSignal = performanceSignal;
        params.emPathType = emPathType;
        params.centeredEMPaths = centeredEMPaths;
        
        
        params.opticsModel = cond.opticsModel;
        params.pupilDiamMm = cond.pupilDiamMm;
        params.performanceClassifier = cond.performanceClassifier;
        
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
        
        params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance);  
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
        variedParamName = 'GaussPooling';
        
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
            'theLegendPosition', [0.30,0.905,0.67,0.08], ...   % custom legend position and size
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
    params.visualizeResponsesWithSpatialPoolingSchemeInVideo = ~true;
    params.visualizeResponses = visualizeResponses;
    params.visualizeOuterSegmentFilters = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = visualizePerformance;
    params.deleteResponseInstances = ~true;
end