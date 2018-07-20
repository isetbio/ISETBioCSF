function run_Paper1FinalConditions2vs3mmPupil
% This is the script used to assess how the final conditions in paper1
% with a 2 vs 3mm pupil compare to the Banks prediction.
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Init condition index
    condIndex = 0;
    
    if (1==2)
    % Our best estimate of mosaic + optics, MLPT inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics, MLPT, 3mm';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosDefaultSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'mlpt';
    examinedCond(condIndex).signal = 'isomerizations';
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).centeredEMPaths = ~true;
    examinedCond(condIndex).frameRate = 10;
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    end
    
    % Our best estimate of mosaic + optics, SVMpool inference engine
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics, MLPT, 2mm'; %'Realistic mosaic/optics, SVM-Template (cos-profile)';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosDefaultSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'mlpt';
    examinedCond(condIndex).signal = 'isomerizations';
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).centeredEMPaths = ~true;
    examinedCond(condIndex).frameRate = 10;
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    
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
        params = getRemainingDefaultParams(params, condIndex, cond.conditionLabel);  
        
        examinedLegends{numel(examinedLegends) + 1} = examinedCond(condIndex).conditionLabel;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Pupil';
        theRatioLims = [0.03 1.5];
        theRatioLims = [0.3 0.8];
        theRatioTicks = [0.3:0.1:1];
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

function params = getRemainingDefaultParams(params, condIndex, conditionLabel)
            
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    
    params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];

    if (contains(conditionLabel, '3mm'))
        params.pupilDiamMm = 3.0;
    elseif (contains(conditionLabel, '2mm'))
        params.pupilDiamMm = 2.0;
    else
        error('condition label does not contain pupil size');
    end
    
    fprintf('>>>>>> \t %d pupil: %f\n', condIndex, params.pupilDiamMm);
    
                
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
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
end
