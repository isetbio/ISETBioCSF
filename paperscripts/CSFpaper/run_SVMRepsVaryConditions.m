function run_SVMRepsVaryConditions
% This is the script used to assess the impact of different # of trials on the SVM-based CSF
%  
    % How to split the computation
    % 16 or 32
    computationInstance = 32;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = ~true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    
    % Optics to use
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    
    params.frameRate = 10; %(10 frames/sec = one 100 msec frame)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
     % Eye movement setup
    params.emPathType = 'random';
    params.centeredEMPaths = true;
    
    % Trials to generate
    params.nTrainingSamples = 1024*16;
    
    params.lowContrast = 0.13;
    params.highContrast =  0.18;
    params.nContrastsPerDirection =  2;
    
    % Trials to use in the classifier - vary this one
    params.performanceTrialsUsed = params.nTrainingSamples;
    
    maxK = 1;
    for k = 1:maxK
        performanceTrialsUsed = params.nTrainingSamples/(2^(k-1));
        examinedCond(maxK+1-k).classifier = 'svmV1FilterBank';
        examinedCond(maxK+1-k).performanceTrialsUsed = performanceTrialsUsed;
        examinedCond(maxK+1-k).legend = sprintf('QPhE SVM, %d trials', performanceTrialsUsed);
    end
    
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
    params.visualizeSpatialPoolingScheme = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for condIndex = 1:numel(examinedCond)
        params.performanceClassifier = examinedCond(condIndex).classifier;
        params.performanceTrialsUsed = examinedCond(condIndex).performanceTrialsUsed;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'SVMTrials';
        theRatioLims = [0.05 0.5];
        theRatioTicks = [0.05 0.1 0.2 0.5];
        generateFigureForPaper(theFigData, examinedInferenceEngineLegends, variedParamName, '', ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end