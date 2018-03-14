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
    
    examinedCond(1).emPathType = 'frozen0';
    examinedCond(1).classifier = 'mlpt';
    examinedCond(1).legend = 'eyes fixed, MLPT';
    
    examinedCond(2).emPathType = 'frozen0';
    examinedCond(2).classifier = 'svm';
    examinedCond(2).legend = 'eyes fixed, SVM';
    
    examinedCond(3).emPathType = 'frozen0';
    examinedCond(3).classifier = 'svmV1FilterBank';
    examinedCond(3).legend = 'eyes fixed, SVM (QPhE)';
    
    examinedCond(4).emPathType = 'random';
    examinedCond(4).classifier = 'svm';
    examinedCond(4).legend = 'drifts+microsaccades, SVM';
    
    examinedCond(5).emPathType = 'random';
    examinedCond(5).classifier = 'svmV1FilterBank';
    examinedCond(5).legend = 'drifts+microsaccades, SVM (QPhE)';


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
    examinedEyeMovementTypeLegends = {};
    for condIndex = 1:numel(examinedCond)
        cond = examinedCond(condIndex);
        params.emPathType = cond.emPathType;
        params.performanceClassifier = cond.classifier;
        examinedEyeMovementTypeLegends{condIndex} = cond.legend;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'EyeMovement';
        generateFigureForPaper(theFigData, examinedEyeMovementTypeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), 'figureType', 'CSF');
    end
end