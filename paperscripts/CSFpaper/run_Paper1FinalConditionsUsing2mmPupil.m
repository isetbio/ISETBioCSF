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
    examinedCond(condIndex).signal = 'isomerizations';
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).centeredEMPaths = ~true;
    examinedCond(condIndex).frameRate = 10;
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    
    % Our best estimate of mosaic + optics, MLPT inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics, MLPT';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'mlpt';
    examinedCond(condIndex).signal = 'isomerizations';
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).centeredEMPaths = ~true;
    examinedCond(condIndex).frameRate = 10;
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    
    % Our best estimate of mosaic + optics, SVMpool inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics, SVM-Template (cos-profile)';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedCond(condIndex).signal = 'isomerizations';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).centeredEMPaths = ~true;
    examinedCond(condIndex).frameRate = 10;
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    
    % Our best estimate of mosaic + optics + eye movements, SVMpool inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics + fixationalEM, SVM-Template (quadrature energy)';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedCond(condIndex).signal = 'isomerizations';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).centeredEMPaths = true;
    examinedCond(condIndex).frameRate = 20;
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    
    
    % Our best estimate of mosaic + optics + eye movements + photocurrents, SVMpool inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics + fixationalEM + pCurrent, SVM-Template (quadrature energy)';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedCond(condIndex).signal = 'photocurrents';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).centeredEMPaths = true;
    examinedCond(condIndex).frameRate = 20;
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        cond = examinedCond(condIndex);
        mosaicName = cond.mosaicName;
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
  
        params.opticsModel = cond.opticsModel;
        params.performanceClassifier = cond.inferenceEngine;
        params.performanceSignal = cond.signal;
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMPaths;
        params.frameRate = cond.frameRate;
        params.responseStabilizationMilliseconds = cond.responseStabilizationMilliseconds;
        params.responseExtinctionMilliseconds = cond.responseExtinctionMilliseconds;
        
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        
        % Update params
        params = getRemainingDefaultParams(params, condIndex);  
        
        examinedLegends{numel(examinedLegends) + 1} = examinedCond(condIndex).conditionLabel;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'VariousParams';
        theRatioLims = [0.03 1.5];
        theRatioTicks = [0.05  0.1 0.2 0.5 1.0];
        formatLabel = 'ComparedToBanks87Photocurrents';  % 'ComparedToBanks87'
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
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
