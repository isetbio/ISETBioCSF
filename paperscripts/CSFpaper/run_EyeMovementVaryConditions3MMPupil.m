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
    classifier = 'svmV1FilterBank'; % 'svm'; % 'svmV1FilterBank'; %'mlpt';
    
    if (strcmp(classifier, 'mlpt'))
            
        condIndex = condIndex+1;  
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = 'mlpt';
        examinedCond(condIndex).legend = 'no eye movements, ideal observer';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;

        condIndex = condIndex+1;  
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'mlpt';
        examinedCond(condIndex).legend = 'fixational eye movements, ideal observer';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;
   
    elseif (strcmp(classifier, 'svm'))

        condIndex = condIndex+1;  
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = 'svm';
        examinedCond(condIndex).legend = 'no eye movements, SVM-PCA';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;


        condIndex = condIndex+1;  
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'svm';
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-PCA';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateSpatialPoolingKernelParams;

    elseif (strcmp(classifier, 'svmV1FilterBank'))
        % Started this late Sunday, Sept 9, to see if there is a difference
        % in photocurrents b/n template and templateQ for 3 mm and best PSF
        % which extend further out in SFs
%         condIndex = condIndex+1;  
%         examinedCond(condIndex).signal = 'photocurrents';
%         examinedCond(condIndex).emPathType = 'frozen0';
%         examinedCond(condIndex).classifier = 'svmV1FilterBank';
%         examinedCond(condIndex).legend = 'no eye movements, SVM-TemplateQ';
%         examinedCond(condIndex).centeredEMpaths = true;
%         examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 21 frames, each 100 msec long)
%         examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%         examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%         examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
        
 % Started this late Sunday, Sept 9, to see if there is a difference
        % in photocurrents b/n template and templateQ for 3 mm and best PSF
        % which extend further out in SFs
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = 'photocurrents';
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-TemplateQ';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 100;
        examinedCond(condIndex).responseExtinctionMilliseconds = 50;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    end
    

    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
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
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    
    % Go
    examinedEyeMovementTypeLegends = {};
    for condIndex = 1:numel(examinedCond)
        cond = examinedCond(condIndex);
        params.performanceSignal = cond.signal;
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
        theRatioLims = [0.05 2];
        theRatioTicks = [0.05 0.1 0.2 0.5 1 2];
        generateFigureForPaper(theFigData, examinedEyeMovementTypeLegends, variedParamName, sprintf('%s_%s_%s',mosaicName, opticsName, classifier), ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end