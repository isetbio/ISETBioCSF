function run_paper2IsomerizationsVsPhotocurrents
% Compute and contrast performance at the level of isomerizations vs photocurrents.
%
% Syntax:
%   run_paper2IsomerizationsVsPhotocurrents
%
% Description:
%    Compute and contrast performance at the level of isomerizations vs photocurrents.
%    This is done in the absence of eye movements and contrasts 
%    1. isomerizations using the mlpt inference engine VS
%    2. isomerizations using the SVM-Template-Linear inference engine VS
%    3. photocurrents using the SVM-Template-Linear inference engine
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
    
    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamcs of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % Assemble conditions list to be examined
    % Init condition index
    condIndex = 0;
    
    if (~computeResponses)
        condIndex = condIndex+1;
        examinedCond(condIndex).label = 'Ideal observer, R*';
        examinedCond(condIndex).performanceClassifier = 'mlpt';
        examinedCond(condIndex).performanceSignal = 'isomerizations';

        condIndex = condIndex+1;
        examinedCond(condIndex).label = 'SVM-Template-Linear, R*';
        examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
        examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
        examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
        examinedCond(condIndex).performanceSignal = 'isomerizations';
    end
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'SVM-Template-Linear, pCurrent';
    examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
    examinedCond(condIndex).performanceSignal = 'photocurrents';
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, computationInstance);
        
        cond = examinedCond(condIndex);
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
        params.performanceClassifier = cond.performanceClassifier;
        params.performanceSignal = cond.performanceSignal;
        
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        
        % Update params
        if (strcmp(cond.label, 'Ideal observer, isomerizations'))
            computePhotocurrents = false;
        else
            computePhotocurrents = true;
        end
        
        params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance);  
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
        variedParamName = 'SignalType';
        theRatioLims = [0.05 1.0];
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
