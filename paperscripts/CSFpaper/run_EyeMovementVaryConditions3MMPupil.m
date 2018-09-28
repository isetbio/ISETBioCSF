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
    
    % Make spatialpooling kernel params struct for SVM-Template with linear activation functiom
    svmTemplateLinearSpatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    svmTemplateLinearSpatialPoolingKernelParams.type = 'V1CosUnit';
    svmTemplateLinearSpatialPoolingKernelParams.activationFunction = 'linear';
    
    % Make spatialpooling kernel params struct for SVM-Template with HW
    % rectifier activation function
    svmTemplateHalfWaveRectSpatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    svmTemplateHalfWaveRectSpatialPoolingKernelParams.type = 'V1CosUnit';
    svmTemplateHalfWaveRectSpatialPoolingKernelParams.activationFunction = 'halfWaveRectifier';
    
    % Make spatialpooling kernel params struct for SVM-Template with FW
    % rectifier activation function
    svmTemplateFullWaveRectSpatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
    svmTemplateFullWaveRectSpatialPoolingKernelParams.type = 'V1CosUnit';
    svmTemplateFullWaveRectSpatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
    
    
    % Define signal to visualize
    signal = 'isomerizations';
    
    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
   
    
    % Use 3 mm to get more of the high frequencies
    params.pupilDiamMm = 3.0;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    

    condIndex = 0;
    classifier = 'svmTemplateLinear'; 'svmTemplateQuadrature'; %  %'mlpt'; 'svm'; 'svmTemplateLinear'; 'svmTemplateHalfWaveRectifier'; 'svmTemplateFullWaveRectifier'; 'svmTemplateQuadrature'
    
    if (strcmp(classifier, 'mlpt'))
            
        condIndex = condIndex+1;
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = classifier;
        examinedCond(condIndex).legend = 'no eye movements, ideal observer';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;

        condIndex = condIndex+1;
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = classifier;
        examinedCond(condIndex).legend = 'fixational eye movements, ideal observer';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
   
    elseif (strcmp(classifier, 'svm'))

        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = classifier;
        examinedCond(condIndex).legend = 'no eye movements, SVM-PCA';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 21 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;


        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = classifier;
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-PCA';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;

    elseif (strcmp(classifier, 'svmTemplateLinear'))
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
        

%         condIndex = condIndex+1;  
%         examinedCond(condIndex).signal = 'photocurrents';
%         examinedCond(condIndex).emPathType = 'random';
%         examinedCond(condIndex).classifier = 'svmV1FilterBank';
%         examinedCond(condIndex).legend = 'fixational eye movements, SVM-TemplateQ';
%         examinedCond(condIndex).centeredEMpaths = true;
%         examinedCond(condIndex).frameRate = 20; %(10 frames/sec, so 1 frames, each 100 msec long)
%         examinedCond(condIndex).responseStabilizationMilliseconds = 100;
%         examinedCond(condIndex).responseExtinctionMilliseconds = 50;
%         examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;

        % SVM-Template with linear activation function
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'zero eye movements, SVM-Template-L';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateLinearSpatialPoolingKernelParams;
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-Template-L';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateLinearSpatialPoolingKernelParams;
        
    elseif (strcmp(classifier, 'svmTemplateHalfWaveRectifier'))
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'zero eye movements, SVM-Template-HW';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateHalfWaveRectSpatialPoolingKernelParams;
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-Template-HW';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateHalfWaveRectSpatialPoolingKernelParams;
     
   elseif (strcmp(classifier, 'svmTemplateFullWaveRectifier'))    
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'zero eye movements, SVM-Template-FW';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateFullWaveRectSpatialPoolingKernelParams;
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-Template-FW';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = svmTemplateFullWaveRectSpatialPoolingKernelParams;
        
    elseif (strcmp(classifier, 'svmTemplateQuadrature'))
      condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'frozen0';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'zero movements, SVM-Template-Q';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
        
        condIndex = condIndex+1;  
        examinedCond(condIndex).signal = signal;
        examinedCond(condIndex).emPathType = 'random';
        examinedCond(condIndex).classifier = 'svmV1FilterBank';
        examinedCond(condIndex).legend = 'fixational eye movements, SVM-Template-Q';
        examinedCond(condIndex).centeredEMpaths = true;
        examinedCond(condIndex).frameRate = 10; %(10 frames/sec, so 1 frames, each 100 msec long)
        examinedCond(condIndex).responseStabilizationMilliseconds = 40;
        examinedCond(condIndex).responseExtinctionMilliseconds = 40;
        examinedCond(condIndex).spatialPoolingKernelParams = defaultSpatialPoolingKernelParams;
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