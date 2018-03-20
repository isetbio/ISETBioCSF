function runSummaryCSF()
     
    [theParams{1}, theLegends{1}] = paramsForBanksCondition();
    [theParams{2}, theLegends{2}] = paramsForBanksGeislerOpticsEccMosaic();
    [theParams{3}, theLegends{3}] = paramsForWvfOpticsEccMosaic();
    [theParams{4}, theLegends{4}] = paramsForSVMQPhE();
    [theParams{5}, theLegends{5}] = paramsForEyeMovements();
    [theParams{6}, theLegends{6}] = paramsForPhotocurrent();
    
    theFigData = {};
    for k = 1:numel(theParams)
        [~,~, theFigData{numel(theFigData)+1}] = run_BanksPhotocurrentEyeMovementConditions(theParams{k});
    end
    
    variedParamName = 'Summary';
    theRatioLims = [0.03 6];
    theRatioTicks = [0.03  0.1 0.3 1 3 5];
    generateFigureForPaper(theFigData, theLegends, variedParamName, 'set1', ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
end

function [params, legend] = paramsForBanksCondition()
    legend = 'Banks ''87';
    mosaicName = 'originalBanks'; 
    opticsName = 'Geisler';
    computationInstance = 0;
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    params.pupilDiamMm = 3.0;
    params.luminancesExamined = 34;
    params.emPathType = 'frozen0';
    params.centeredEMpaths = true;
    
    params.cyclesPerDegreeExamined = params.cyclesPerDegreeExamined(1:end-1);
     
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    params.performanceClassifier = 'mlpt';
    params.performanceSignal = 'isomerizations';
    
    params = computeParams(params);
end

function [params, legend] = paramsForBanksGeislerOpticsEccMosaic()
    legend = 'realistic mosaic';
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    opticsName = 'Geisler';
    computationInstance = 0;
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    params.pupilDiamMm = 3.0;
    params.luminancesExamined = 34;
    params.emPathType = 'frozen0';
    params.centeredEMpaths = true;
     
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    params.performanceClassifier = 'mlpt';
    params.performanceSignal = 'isomerizations';
    
    params = computeParams(params);
end

function [params, legend] = paramsForWvfOpticsEccMosaic()
    legend = 'wvf-based optics';
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    computationInstance = 0;
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    params.pupilDiamMm = 3.0;
    params.luminancesExamined = 34;
    params.emPathType = 'frozen0';
    params.centeredEMpaths = true;
     
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    params.performanceClassifier = 'mlpt';
    params.performanceSignal = 'isomerizations';
    
    params = computeParams(params);
    
end

function [params, legend] = paramsForSVMQPhE()
    legend = 'SVM QPhE classifier';
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    computationInstance = 0;
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    params.pupilDiamMm = 3.0;
    params.luminancesExamined = 34;
    params.emPathType = 'frozen0';
    params.centeredEMpaths = true;
     
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    params.performanceClassifier = 'svmV1FilterBank';
    params.performanceSignal = 'isomerizations';
    
    params = computeParams(params);
end

function [params, legend] = paramsForEyeMovements()
    legend = 'eye movements';
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    computationInstance = 0;
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    params.pupilDiamMm = 3.0;
    params.luminancesExamined = 34;
    params.emPathType = 'random';
    params.centeredEMpaths = true;
     
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    params.performanceClassifier = 'svmV1FilterBank';
    params.performanceSignal = 'isomerizations';
    
    params = computeParams(params);
end



function [params, legend] = paramsForPhotocurrent()
    legend = 'photocurrent';
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    computationInstance = 0;
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    params.pupilDiamMm = 3.0;
    params.luminancesExamined = 34;
    params.emPathType = 'random';
    params.centeredEMpaths = true;
    
    params.frameRate = 20; %(2 frames)
    params.responseStabilizationMilliseconds = 100;
    params.responseExtinctionMilliseconds = 50;
    
    params.performanceClassifier = 'svmV1FilterBank';
    params.performanceSignal = 'photocurrents';
    
    params = computeParams(params);
end


function params = computeParams(params)
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
end