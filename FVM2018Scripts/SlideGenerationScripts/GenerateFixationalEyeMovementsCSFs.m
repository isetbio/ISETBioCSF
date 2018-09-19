function GenerateFixationalEyeMovementsCSFs

    % Script to generate the slide with the Banks'87 ideal and human observer data
    
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    computeResponses = false;
    findPerformance = false;
    
    % Init condition index
    condIndex = 0;
    

    % ISETBio simulation using Banks conditions
    condIndex = condIndex+1;
    examinedConds(condIndex).conditionLabel = 'Banks mosaic/optics, ideal observer, noEM';
    examinedConds(condIndex).mosaicName = 'originalBanks';
    examinedConds(condIndex).opticsModel = 'Geisler';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).inferenceEngine = 'mlpt';
    examinedConds(condIndex).emPathType = 'frozen0';
    examinedConds(condIndex).centeredEMPaths = ~true;
    examinedConds(condIndex).frameRate = 10;
    examinedConds(condIndex).responseStabilizationMilliseconds =  40;
    examinedConds(condIndex).responseExtinctionMilliseconds = 40;
    
    
    % Our best estimate of mosaic + optics, MLPT inference engine
    condIndex = condIndex+1;                                
    examinedConds(condIndex).conditionLabel = 'Realistic mosaic/optics, ideal observer, noEM';
    examinedConds(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection';
    examinedConds(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).inferenceEngine = 'mlpt';
    examinedConds(condIndex).emPathType = 'frozen0';
    examinedConds(condIndex).centeredEMPaths = ~true;
    examinedConds(condIndex).frameRate = 10;
    examinedConds(condIndex).responseStabilizationMilliseconds =  40;
    examinedConds(condIndex).responseExtinctionMilliseconds = 40;
    
    
    % Our best estimate of mosaic + optics, MLPT inference engine
    condIndex = condIndex+1;                                
    examinedConds(condIndex).conditionLabel = 'Realistic mosaic/optics, SVM-Template-L, noEM';
    examinedConds(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection';
    examinedConds(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedConds(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedConds(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedConds(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).emPathType = 'frozen0';
    examinedConds(condIndex).centeredEMPaths = ~true;
    examinedConds(condIndex).frameRate = 10;
    examinedConds(condIndex).responseStabilizationMilliseconds =  40;
    examinedConds(condIndex).responseExtinctionMilliseconds = 40;
    
    if (1==2)
    % Our best estimate of mosaic + optics, SVM-Template-linear inference engine
    condIndex = condIndex+1;                                
    examinedConds(condIndex).conditionLabel = 'Realistic mosaic/optics, SVM-Template-L, fixEM';
    examinedConds(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; 
    examinedConds(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedConds(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedConds(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
    examinedConds(condIndex).emPathType = 'random';
    examinedConds(condIndex).centeredEMPaths = true;
    examinedConds(condIndex).frameRate = 20;
    examinedConds(condIndex).responseStabilizationMilliseconds = 100;
    examinedConds(condIndex).responseExtinctionMilliseconds = 50;
    end
    
    
    % Our best estimate of mosaic + optics, SVM-Template-quadrature engine
    condIndex = condIndex+1;                                
    examinedConds(condIndex).conditionLabel = 'Realistic mosaic/optics, SVM-Template-Q, fixEM';
    examinedConds(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; 
    examinedConds(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedConds(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedConds(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedConds(condIndex).emPathType = 'random';
    examinedConds(condIndex).centeredEMPaths = true;
    examinedConds(condIndex).frameRate = 20;
    examinedConds(condIndex).responseStabilizationMilliseconds = 100;
    examinedConds(condIndex).responseExtinctionMilliseconds = 50;

    
    % Go
    [examinedLegends, theFigData, pupilDiamMm] = runConditions(...
        examinedConds, computationInstance, computeResponses, findPerformance);
    
    if (makeSummaryFigure)
        variedParamName = 'FixationalEyeMovements';
        theRatioLims = [0.02 2.0];
        theRatioTicks = [0.05  0.1 0.2 0.5 1.0];
        formatLabel = '_InferenceEngineCombo'; 
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'figDataIndicesToDisplay', [1 numel(examinedConds)], ...
            'showSubjectData', (pupilDiamMm == 2), ...
            'showSubjectMeanData', false, ...
            'showLegend', ~true, ...
            'plotFirstConditionInGray', true, ...
            'showBanksPaperIOAcurves', ~true, ...
            'showOnly23CDM2IOAcurve', ~true, ...
            'plotRatiosOfOtherConditionsToFirst', false, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    

    
end

function [examinedLegends, theFigData, pupilDiamMm] = runConditions(examinedConds, computationInstance, computeResponses, findPerformance)
    examinedLegends = {};
    theFigData = {};
    
    for condIndex = 1:numel(examinedConds)
        cond = examinedConds(condIndex);
        mosaicName = cond.mosaicName;
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        
        % Update params
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
        
        params = getRemainingDefaultParams(params, condIndex, cond.conditionLabel, computeResponses, findPerformance);  
        
        % Update returned items
        pupilDiamMm(condIndex) = params.pupilDiamMm;
        examinedLegends{numel(examinedLegends) + 1} = cond.conditionLabel;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (any(diff(pupilDiamMm) ~= 0))
        pupilDiamMm
        error('PuliDiamMM are different for different conditions');
    else
        pupilDiamMm = pupilDiamMm(1);
    end
    
end


function params = getRemainingDefaultParams(params, condIndex, conditionLabel, computeResponses, findPerformance)
            
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    
    if contains(conditionLabel,'Banks mosaic/optics')
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50];
    else
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];
    end
    
    if (strcmp(conditionLabel, 'Banks mosaic/optics, MLPT, 3mm'))
        params.pupilDiamMm = 3.0;
    else
        params.pupilDiamMm = 2.0;
    end
    
    fprintf('>>>>>> \t %d pupil: %f\n', condIndex, params.pupilDiamMm);
    
                
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
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
    params.findPerformance = findPerformance;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
end

