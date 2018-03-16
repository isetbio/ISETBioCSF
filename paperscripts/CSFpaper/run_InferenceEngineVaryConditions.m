function run_InferenceEngineVaryConditions
% This is the script used to assess the impact of different inference engines on the CSF
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
    % Optics models we tested
    examinedInferenceEngines = {...
        'mlpt' ...
        'svm' ...
        'svmV1FilterBank' ...
    };
    examinedInferenceEngineLegends = {...
        'MLPT' ...
        'SVM' ...
        'SVM (QPhE)' ...
    };


    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for engineIndex = 1:numel(examinedInferenceEngines)
        params.performanceClassifier = examinedInferenceEngines{engineIndex};
        params.opticsModel = opticsName;
        [~,~, theFigData{engineIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'InferenceEngine';
        theRatioLims = [0.05 0.5];
        theRatioTicks = [0.05 0.1 0.2 0.5];
        generateFigureForPaper(theFigData, examinedInferenceEngineLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

