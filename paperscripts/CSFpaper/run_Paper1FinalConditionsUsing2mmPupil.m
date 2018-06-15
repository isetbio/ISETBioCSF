function run_Paper1FinalConditionsUsing2mmPupil
% This is the script used to assess how the final conditions in paper1
% with a 2 mm pupil compare to the Banks prediction.
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    

    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    
    condIndex = 0;
    
    % Original Banks computation
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Banks mosaic and optics, MLPT';
    examinedCond(condIndex).mosaicName = 'originalBanks';
    examinedCond(condIndex).opticsModel = 'Geisler';
    examinedCond(condIndex).inferenceEngine = 'mlpt';
    
    % Our best estimate of mosaic + optics, MLPT inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic and optics, MLPT';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'mlpt';
    
    % Our best estimate of mosaic + optics, SVMpool inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic and optics, SVMpool';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        cond = examinedCond(condIndex);
        mosaicName = cond.mosaicName;
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        params.opticsModel = cond.opticsModel;
        params.performanceClassifier = cond.inferenceEngine;
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        params = getRemainingDefaultParams(params, condIndex);  
        examinedLegends{numel(examinedLegends) + 1} = examinedCond(condIndex).conditionLabel;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'VariousParams';
        theRatioLims = [0.05 1.5];
        theRatioTicks = [0.05  0.1 0.2 0.5 1.0];
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, 'ComparedToBanks87', ...
            'figureType', 'CSF', ...
            'showSubjectData', true, ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

function params = getRemainingDefaultParams(params, condIndex)

    % All condsexained  with 2 mm pupil
    params.pupilDiamMm = 2.0;
            
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    if (condIndex == 1)
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50];
    else
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];
    end
    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;

    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;
    
                
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
end
