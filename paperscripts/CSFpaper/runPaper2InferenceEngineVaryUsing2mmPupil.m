function runPaper2InferenceEngineVaryUsing2mmPupil
    
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;

    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Init condition index
    condIndex = 0;
    
    if (1==2)
    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics + fixationalEM + isomerizations, SVM-Template (quadrature)';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'svmV1FilterBank';
    examinedCond(condIndex).signal = 'isomerizations'; % 'photocurrents';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).centeredEMPaths = true;
    examinedCond(condIndex).frameRate = 20;
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    end

    condIndex = condIndex+1;
    examinedCond(condIndex).conditionLabel = 'Realistic mosaic/optics + fixationalEM + isomerizations, SVM-V1ensemble (quadrature)';
    examinedCond(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedCond(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedCond(condIndex).inferenceEngine = 'svmV1FilterEnsemble';
    examinedCond(condIndex).signal = 'isomerizations'; % 'photocurrents';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).centeredEMPaths = true;
    examinedCond(condIndex).frameRate = 20;
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    
    
    ensembleFilterParams = struct(...
            'spatialPositionsNum',  9, ...
            'cyclesPerRFs', [2.5], ... % [1.0 1.5 2.0 2.5], ...
            'orientations', [0]);
   
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
        elseif (strcmp(params.performanceClassifier,'svmV1FilterEnsemble'))
            fNames = fieldnames(ensembleFilterParams);
            for fNameIndex = 1:numel(fNames)
                fName = fNames{fNameIndex};
                params.spatialPoolingKernelParams.(fName) = ensembleFilterParams.(fName);
            end
            % Remove
            params.parforWorkersNumForClassification = min([12 2*params.parforWorkersNumForClassification])
        end
        
        % Update params
        params = getRemainingDefaultParams(params, condIndex, cond.conditionLabel); 
        % Remove
        params.findPerformance = true
        
        
        examinedLegends{numel(examinedLegends) + 1} = examinedCond(condIndex).conditionLabel;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
        
    end
    
    if (makeSummaryFigure)
        variedParamName = 'VariousParams';
        theRatioLims = [0.02 2.0];
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

function params = getRemainingDefaultParams(params, condIndex, conditionLabel)
            
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
