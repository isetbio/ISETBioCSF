function params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, computationInstance)

    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrectionAndMacularPigment';
    opticsName = 'ThibosAverageSubject3MMPupil';
    
    % Optics params
    params.opticsModel = opticsName;
    params.blur = true;                                 % employ optics
    params.apertureBlur = true;                         % employ cone aperture blur
    params.wavefrontSpatialSamples = 261*2+1;           % This gives us an OTF sampling of 1.003 c/deg
    params.opticalImagePadSizeDegs = [];
    params.pupilDiamMm = pupilDiamMm;                           % 3 is more appropriate for a 100 cd/m2 mean scene luminance
    
    % Mosaic params
    mosaicParams = getParamsForMosaicWithLabel(mosaicName);
    
    fNames = fieldnames(mosaicParams);
    for k = 1:numel(fNames)
        params.(fNames{k}) = mosaicParams.(fNames{k});
    end
    
    % Temporal params that are best for combo isomerizations/photocurrent responses
    params.frameRate = 10;                      %(10 frames/sec, so 1, 100 msec long frame)
    params.singleBinTemporalWindowIfPossible = false;
    params.responseStabilizationMilliseconds = 400;  % 500 msec is really needed to avoid being on settling part of the photocurrent response
    params.responseExtinctionMilliseconds = 100;
        
    params.integrationTimeMilliseconds = integrationTimeMilliseconds;
    
    if (params.integrationTimeMilliseconds == 2.5)
        params.secondsToInclude = 150/1000;
        params.secondsToIncludeOffset = 75/1000;
                                 
    elseif (params.integrationTimeMilliseconds == 4.0)
        params.secondsToInclude = 148/1000;
        params.secondsToIncludeOffset = 72/1000;
        
    elseif (params.integrationTimeMilliseconds == 5.0)
        params.secondsToInclude = 150/1000;
        params.secondsToIncludeOffset = 75/1000; 
    else
        error('integrationTimeMilliseconds is not 2.5 4 or 5: %f', integrationTimeMilliseconds);
    end
    
    
    % Fixational eye movements type: 'random'; 'randomNoSaccades' 'frozen0';
    params.emPathType = 'frozen0';  
    params.centeredEMPaths = true;
    
     % Stimulus luminance
    params.luminancesExamined = [34];
    
    % Stimulus size
    params.imagePixels = 512;           % stimuli will be 512x512 pixels
    
    % Grating orientation
    params.stimulusOrientationDegs = 0;
    
    % Conrast axis sampling
    params.coneContrastDirection = 'L+M+S';
    params.lowContrast = 0.00001*3;
    params.highContrast =  1.0;
    params.nContrastsPerDirection =  20;
    
    % Trials to generate
    trialsToCompute = 128;
    params.nTrainingSamples = trialsToCompute; 
    
    % Freeze noise for repeatable results
    params.freezeNoise = true;
    
    % Use all trials in the classifier (Specify [] to use all available trials)
    params.performanceTrialsUsed = [];
    
    % Performance signal
    params.performanceSignal = 'isomerizations';
    
    % 'SVM-Template-Linear'
    params.performanceClassifier = 'svmV1FilterBank';
    params.spatialPoolingKernelParams.type = 'V1CosUnit';
    params.spatialPoolingKernelParams.activationFunction = 'linear';
    
    % Split computations and specify RAM memory
    if (computationInstance == 0)
        % All mosaic sizes in 1 MATLAB session
        params.ramPercentageEmployed = 1.2; 
        params.cyclesPerDegreeExamined =  [2 4 8 16 32 50 60]; 
        params.parforWorkersNumForClassification = 3;
    elseif (computationInstance  == 1)
        % Largest mosaic
        params.ramPercentageEmployed = 1.2; 
        params.cyclesPerDegreeExamined =  [2];
        params.parforWorkersNumForClassification = 3;
    elseif (computationInstance  == 2)
        % Second largest mosaic
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined =  [4];
        params.parforWorkersNumForClassification = 3;
    elseif (computationInstance  == 3)
        % All other sizes
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined =  [8 16 32 50 60];
        params.parforWorkersNumForClassification = 3;
    elseif (computationInstance  == 4)
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined = [4];
        params.parforWorkersNumForClassification = 3;
    elseif (computationInstance  == 8)
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined = [8];
        params.parforWorkersNumForClassification = 6;
    elseif (computationInstance  == 16)
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined = [16];
        params.parforWorkersNumForClassification = 10;
    elseif (computationInstance  == 32)
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined = [32];
        params.parforWorkersNumForClassification = 12;
    elseif (computationInstance  == 50)
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined = [50];
        params.parforWorkersNumForClassification = 12;
    elseif (computationInstance  == 60)
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined = [60];
        params.parforWorkersNumForClassification = 12;
    else
        error('computational instance (%d) not found', computationalInstance);
    end

    % Adjust params.parforWorkersNumForClassification depending on computer
    [numberOfCores, ramSizeGBytes, ~] = determineSystemResources(false);
    if (ramSizeGBytes > 200)
        params.parforWorkersNumForClassification = min([numberOfCores 2*params.parforWorkersNumForClassification]);
    elseif (ramSizeGBytes > 100)
        params.parforWorkersNumForClassification = min([numberOfCores params.parforWorkersNumForClassification]);
    else
        params.parforWorkersNumForClassification = min([numberOfCores 1]);
    end
    
end