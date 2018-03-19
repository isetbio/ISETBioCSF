function run_EyeMovementVaryConditions
% This is the script used to assess the impact of different types of eye movements on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 1;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = ~true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    
    % Optics to use
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
    
    % Response duration params
    params.frameRate = 20; %(2 frames)
    params.responseStabilizationMilliseconds = 100;
    params.responseExtinctionMilliseconds = 50;
    
   % examinedCond(1).emPathType = 'frozen0';
   % examinedCond(1).classifier = 'mlpt';
   % examinedCond(1).legend = 'no eye movements, MLPT';
    
   % examinedCond(1).emPathType = 'frozen0';
   % examinedCond(1).classifier = 'svm';
   % examinedCond(1).legend = 'no eye movements, SVM';
    
    examinedCond(1).emPathType = 'frozen0';
    examinedCond(1).classifier = 'svmV1FilterBank';
    examinedCond(1).legend = 'no eye movements, SVM (QPhE)';
    examinedCond(1).centeredEMpaths = true;
    
    examinedCond(2).emPathType = 'random';
    examinedCond(2).classifier = 'svm';
    examinedCond(2).legend = 'drifts+\mu-saccades, SVM';
    examinedCond(2).centeredEMpaths = true;
    
    examinedCond(3).emPathType = 'random';
    examinedCond(3).classifier = 'svmV1FilterBank';
    examinedCond(3).legend = 'drifts+\mu-saccades, SVM (QPhE)';
    examinedCond(3).centeredEMpaths = true;

    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = ~true;
    params.deleteResponseInstances = ~true;
    
    % When computing responses only do the first  conditions 1 and 3
    if (params.computeResponses)
        examinedCond(2) = examinedCond(3);
        examinedCond = examinedCond(1:2);
    end
    
    % Go
    examinedEyeMovementTypeLegends = {};
    for condIndex = 1:numel(examinedCond)
        cond = examinedCond(condIndex);
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMpaths;
        params.performanceClassifier = cond.classifier;
        examinedEyeMovementTypeLegends{condIndex} = cond.legend;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'EyeMovement';
        theRatioLims = [0.1 2];
        theRatioTicks = [0.1 0.2 0.5 1 2];
        generateFigureForPaper(theFigData, examinedEyeMovementTypeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' B ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end