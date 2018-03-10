function run_ClassificationMethodVary
% This is the script used to assess the impact of different optics models on the CSF
%

    %'mlpt'% 'svm'; 'svmV1FilterBank';
    params.performanceClassifier = 'svm'; % 'svmV1FilterBank';
     
    % Mosaic type to employ
    mosaicType = 'ISETbioHexEccBasedLMSrealistic';
    
    % Optics to run
    params.opticsModel = 'ThibosBestPSFSubject3MMPupil';

    params.pupilDiamMm = 3.0;   % 3 is more appropriate for a 100 cd/m2 mean scene luminance
    
    params.luminancesExamined = [34];
       
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
    
    params.visualizeKernelTransformedSignals = true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % How to split the computation
    computationInstance = 0;        % ALL mosaics
    %computationInstance = 1;        % LARGEST mosaic
    %computationInstance = 2;        % second LARGEST mosaic
    %computationInstance = 3;        % ALL except the two LARGEST mosaic
    %computationInstance = 4;
    params = getFixedParamsForOpticsImpactExperiment(params,computationInstance, mosaicType);

    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);

end

function params = getFixedParamsForOpticsImpactExperiment(params, computationInstance, mosaicType)
    params.blur = true;
    params.apertureBlur = true;
    
    params.imagePixels = 512;
    params.wavefrontSpatialSamples = 261*2+1;     % This gives us an OTF sampling of 1.003 c/deg
    params.minimumOpticalImagefieldOfViewDegs = 1.0;
    
    
    % 'random'; 'frozen0';
    params.emPathType = 'frozen0'; %random'; %'random';     
    params.centeredEMPaths = false;
   
    % Use a subset of the trials. Specify [] to use all available trials
    params.nTrainingSamples = 1024;
  
    % Mosaic params
    mosaicParams = getParamsForMosaicWithLabel(mosaicType);
    
    fNames = fieldnames(mosaicParams);
    for k = 1:numel(fNames)
        params.(fNames{k}) = mosaicParams.(fNames{k});
    end
    
    params.integrationTimeMilliseconds =  5.0;
    params.resamplingFactor = [];   % Empty indicates that c_BanksEtAlPhotocurrentAndEyeMovements
                                    % will choose a resamplingFactor based
                                    % on the mosaic size, smaller mosaics
                                    % will get higher resamplingFactor
    
    params.maxGridAdjustmentIterations = []; % Empty indicates that c_BanksEtAlPhotocurrentAndEyeMovements
                                    % will choose maxGridAdjustmentIterations based
                                    % on the mosaic size, smaller mosaics
                                    % will get higher rmaxGridAdjustmentIterations
        
    % response params
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    % Conditions 
    params.lowContrast = 0.00001*3;
    params.highContrast =  1.0;
    params.nContrastsPerDirection =  20;
    
    % 'isomerizations', 'photocurrents'
    params.performanceSignal = 'isomerizations';
    

    params.performanceTrialsUsed = params.nTrainingSamples;
    
    % SVM-related
    params.useRBFSVMKernel = false;
    % Spatial pooling kernel parameters
    params.spatialPoolingKernelParams.type = 'V1QuadraturePair';  % Choose between 'V1CosUnit' 'V1SinUnit' 'V1QuadraturePair';
    params.spatialPoolingKernelParams.activationFunction = 'energy';  % Choose between 'energy' and 'fullWaveRectifier'
    params.spatialPoolingKernelParams.adjustForConeDensity = false;
    params.spatialPoolingKernelParams.temporalPCAcoeffs = Inf;  % Inf, results in no PCA, just the raw time series
    params.spatialPoolingKernelParams.shrinkageFactor = 1.0;  % > 1, results in expansion, < 1 results in shrinking 
    
    
    % Freeze noise for repeatable results
    params.freezeNoise = true;
    
    % Split computations and specify RAM memory
    if (computationInstance == 0)
        % All conditions in 1 MATLAB session
        params.ramPercentageEmployed = 1.2; 
        params.cyclesPerDegreeExamined =  [2 4 8 16 32 50 60]; 
    elseif (computationInstance  == 1)
        % Largest mosaic in session 1 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 1.2; 
        params.cyclesPerDegreeExamined =  [2];
    elseif (computationInstance  == 2)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined =  [4];
    elseif (computationInstance  == 3)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined =  [8 16 32 50 60];
    elseif (computationInstance  == 4)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 1.2;  
        params.cyclesPerDegreeExamined =  [50];
    else
        error('computational instance must be 0, 1 or 2');
    end
end

