function run_paper2PupilEffectOnIsomerizations
% Contrast performance at the level of isomerizations for 2 vs 3 mm pupils
%
% Syntax:
%   run_paper2PupilEffectOnIsomerizations
%
% Description:
%    Compute and contrast performance at the level of isomerizations for 2 vs 3 mm pupils
%    This is done in the absence of eye movements using the SVM-Template-Linear inference
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
    findPerformance = ~true;
    visualizePerformance = true;
    
    
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
    examinedCond(condIndex).label = 'cone exc., 3mm pupil';
    examinedCond(condIndex).performanceClassifier = 'mlpt'; % 'svmV1FilterBank';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
    examinedCond(condIndex).performanceSignal = 'isomerizations';
    examinedCond(condIndex).pupilDiamMm = 3.0;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'cone exc., 2mm pupil';
    examinedCond(condIndex).performanceClassifier = 'mlpt'; %'svmV1FilterBank';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
    examinedCond(condIndex).performanceSignal = 'isomerizations';
    examinedCond(condIndex).pupilDiamMm = 2.0;
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        
        cond = examinedCond(condIndex);
        
        params = getCSFPaper2DefaultParams(cond.pupilDiamMm, integrationTimeMilliseconds,  frameRate, stimulusDurationInSeconds, computationInstance);
        
        % Use 1028 vs 1024 trials to differentiate results from old photocurrent computation
        params.nTrainingSamples = 1030;
        
        % Try out for subset of SFs
        params.cyclesPerDegreeExamined = [4 8 12 16 24 32 50 60];
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
        params.performanceClassifier = cond.performanceClassifier;
        params.performanceSignal = cond.performanceSignal;
        params.emPathType = 'frozen0';
        
        params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
        params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        
        % Update params
        computePhotocurrents = true;
        
        
        params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance);  
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
        variedParamName = 'PupilSizeIsomerizations';
        
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
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = visualizePerformance;
    params.deleteResponseInstances = ~true;
end
