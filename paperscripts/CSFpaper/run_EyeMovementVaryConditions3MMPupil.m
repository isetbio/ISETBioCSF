function run_EyeMovementVaryConditions3MMPupil
% This is the script used to assess the impact of different types of eye movements on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; 
    
    % Optics to use
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    defaultSpatialPoolingKernelParams = params.spatialPoolingKernelParams;
    
    % Make spatialpooling kernel params struct for SVM-Template
    svmTemplateSpatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    svmTemplateSpatialPoolingKernelParams.type = 'V1CosUnit';
    svmTemplateSpatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
    
    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
   
    
    % Use 3 mm to get more of the high frequencies
    params.pupilDiamMm = 3.0;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    

    condIndex = 0;
    
    
    if (1==2)
    condIndex = condIndex+1;  
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).classifier = 'mlpt';
    examinedCond(condIndex).legend = 'no eye movements, MLPT';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 21 frames, each 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
    
    condIndex = condIndex+1;  
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'mlpt';
    examinedCond(condIndex).legend = 'no eye movements, MLPT';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 21 frames, each 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
   
    else 
%     condIndex = condIndex+1;  
%     examinedCond(condIndex).emPathType = 'frozen0';
%     examinedCond(condIndex).classifier = 'svm';
%     examinedCond(condIndex).legend = 'no eye movements, SVM-PCA';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 21 frames, each 100 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%     examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
%     
    condIndex = condIndex+1;  
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'no eye movements, SVM-Template';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
    
    
%     condIndex = condIndex+1;  
%     examinedCond(condIndex).emPathType = 'random';
%     examinedCond(condIndex).classifier = 'svm';
%     examinedCond(condIndex).legend = 'fixational eye movements, SVM-PCA';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 1 frames, each 100 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%     examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
%   
    condIndex = condIndex+1;  
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'fixational eye movements, SVM-Template';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;

    condIndex = condIndex+1;  
    examinedCond(condIndex).emPathType = 'random';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'fixational eye movements, SVM-QTemplate';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 40;
    examinedCond(condIndex).responseExtinctionMilliseconds = 40;
    examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    end
    
%     condIndex = condIndex+1;  
%     examinedCond(condIndex).emPathType = 'frozen0';
%     examinedCond(condIndex).classifier = 'svmV1FilterBank';
%     examinedCond(condIndex).legend = 'no eye movements, SVM-QTemplate';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%     examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
%     
%     condIndex = condIndex+1;  
%     examinedCond(condIndex).emPathType = 'frozen0';
%     examinedCond(condIndex).classifier = 'svm';
%     examinedCond(condIndex).legend = 'no eye movements, SVM-PCA - 2 frames';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 40;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 40;
%     examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    
    

    
%     condIndex = condIndex+1;  
%     examinedCond(condIndex).emPathType = 'frozen0';
%     examinedCond(condIndex).classifier = 'svmV1FilterBank';
%     examinedCond(condIndex).legend = 'no eye movements, SVM-Template';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%    examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;

    if (1==2)
    condIndex = condIndex+1;  
    examinedCond(condIndex).emPathType = 'frozen0';
    examinedCond(condIndex).classifier = 'svmV1FilterBank';
    examinedCond(condIndex).legend = 'no eye movements, SVM-QTemplate';
    examinedCond(condIndex).centeredEMpaths = true;
    examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
    examinedCond(condIndex).responseStabilizationMilliseconds = 100;
    examinedCond(condIndex).responseExtinctionMilliseconds = 50;
    examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    end
    
    
%     condIndex = condIndex+1;
%     examinedCond(condIndex).emPathType = 'random';
%     examinedCond(condIndex).classifier = 'svmV1FilterBank';
%     examinedCond(condIndex).legend = 'fixational eye movements, SVM-Template';
%     examinedCond(condIndex).centeredEMpaths = ~true;
%     examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%     examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
%     
% 
%     condIndex = condIndex+1;
%     examinedCond(condIndex).emPathType = 'random';
%     examinedCond(condIndex).classifier = 'svmV1FilterBank';
%     examinedCond(condIndex).legend = 'fixational eye movements, SVM-QTemplate';
%     examinedCond(condIndex).centeredEMpaths = ~true;
%     examinedCond(condIndex).frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%     examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;

    
     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).emPathType = 'random';
%     examinedCond(condIndex).classifier = 'svmV1FilterEnsemble';
%     examinedCond(condIndex).legend = 'drifts+\mu-sacc. (rnd), SVM (v1RF-match QPhE(1.0))';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 10;  %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 40;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 40; 
%     examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
%     examinedCond(condIndex).spatialPoolingKernelParams.spatialPositionsNum = 9;
%     examinedCond(condIndex).spatialPoolingKernelParams.cyclesPerRFs = 1;
%     examinedCond(condIndex).spatialPoolingKernelParams.orientations = 0;
%     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).emPathType = 'random';
%     examinedCond(condIndex).classifier = 'svmV1FilterEnsemble';
%     examinedCond(condIndex).legend = 'drifts+\mu-sacc. (rnd), SVM (v1RF-match QPhE(1.5))';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 10;  %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 40;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 40; 
%     examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
%     examinedCond(condIndex).spatialPoolingKernelParams.spatialPositionsNum = 9;
%     examinedCond(condIndex).spatialPoolingKernelParams.cyclesPerRFs = 1.5;
%     examinedCond(condIndex).spatialPoolingKernelParams.orientations = 0;
%     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).emPathType = 'random';
%     examinedCond(condIndex).classifier = 'svmV1FilterEnsemble';
%     examinedCond(condIndex).legend = 'drifts+\mu-sacc. (rnd), SVM (v1RF-match QPhE(2.0))';
%     examinedCond(condIndex).centeredEMpaths = true;
%     examinedCond(condIndex).frameRate = 10;  %(20 frames/sec, so 2 frames, each 50 msec long)
%     examinedCond(condIndex).responseStabilizationMilliseconds = 40;
%     examinedCond(condIndex).responseExtinctionMilliseconds = 40; 
%     examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
%     examinedCond(condIndex).spatialPoolingKernelParams.spatialPositionsNum = 9;
%     examinedCond(condIndex).spatialPoolingKernelParams.cyclesPerRFs = 2.0;
%     examinedCond(condIndex).spatialPoolingKernelParams.orientations = 0;
    
    
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
    params.visualizeMosaicWithFirstEMpath = true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    
    % Go
    examinedEyeMovementTypeLegends = {};
    for condIndex = 1:numel(examinedCond)
        cond = examinedCond(condIndex);
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMpaths;
        params.frameRate = cond.frameRate;
        params.responseStabilizationMilliseconds = cond.responseStabilizationMilliseconds;
        params.responseExtinctionMilliseconds = cond.responseExtinctionMilliseconds;
        params.performanceClassifier = cond.classifier;
        params.spatialPoolingKernelParams = cond.spatialPoolingKernelParams;

        examinedEyeMovementTypeLegends{condIndex} = cond.legend;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'EyeMovement';
        theRatioLims = [0.1 2];
        theRatioTicks = [0.1 0.2 0.5 1 2];
        generateFigureForPaper(theFigData, examinedEyeMovementTypeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end