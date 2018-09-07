function hFig = runConditionToVisualizePsychometricCurve
    
    computationInstance = 0;
    
    % Mosaic to use (ecc-based cone density AND cone efficiency)
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; 
    
    % Optics to use
    opticsName = 'ThibosAverageSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    
    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = true;
    
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    [~,~, theFigData{1}] = run_BanksPhotocurrentEyeMovementConditions(params);
    
    % Generate CSF paper
    variedParamName = '';
    theRatioLims = [0.05 0.5];
    theRatioTicks = [0.05  0.1 0.2 0.5];
    examinedInferenceEngineLegend = params.performanceClassifier;
    hFig = generateFigureForPaper(theFigData, examinedInferenceEngineLegend, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotUsingLargeBlueDisks', true, ...
            'showLegend', false, ...
            'plotFirstConditionInGray', ~true, ...
            'plotRatiosOfOtherConditionsToFirst', ~true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
    );
        
end

