function run_EyeMovementVaryConditions
% This is the script used to assess the impact of different types of eye movements on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
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
    
    condIndex = 0;
    
    
   % examinedCond(condIndex).emPathType = 'frozen0';
   % examinedCond(condIndex).classifier = 'mlpt';
   % examinedCond(condIndex).legend = 'no eye movements, MLPT';
    
   % examinedCond(condIndex).emPathType = 'frozen0';
   % examinedCond(condIndex).classifier = 'svm';
   % examinedCond(condIndex).legend = 'no eye movements, SVM';
    
    condIndex = condIndex+1;
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'no eye movements, SVM (QPhE)';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    
%     condIndex = condIndex+1;
%     examinedCond(condIndex).emPathType = 'frozen0';
%     examinedCond(condIndex).classifier = 'svmV1FilterBank';
%     examinedCond(condIndex).legend = 'no eye movements(2), SVM (QPhE)';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 10;
%     examinedCond(condIndex).responseStabilizationMilliseconds = 40;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    

    condIndex = condIndex+1;
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'svm';
    examinedCond(condIndex).legend = 'drifts+\mu-saccades (origin), SVM';
    examinedCond(condIndex).centeredEMpaths = ~true;
    examinedCond(condIndex).frameRate = 10;  %(10 frames/sec, so 1 frame, 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'drifts+\mu-saccades (origin), SVM (QPhE)';
    examinedCond(condIndex).centeredEMpaths = ~true;
    examinedCond(condIndex).frameRate = 10;  %(10 frames/sec, so 1 frame, 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'svm';
    examinedCond(condIndex).legend = 'drifts+\mu-saccades (random), SVM';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'drifts+\mu-saccades (random), SVM (QPhE)';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 20;  %(20 frames/sec, so 2 frames, each 50 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
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
        params.frameRate = cond.frameRate;
        params.responseStabilizationMilliseconds = cond.responseStabilizationMilliseconds;
        params.responseExtinctionMilliseconds = cond.responseExtinctionMilliseconds;
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