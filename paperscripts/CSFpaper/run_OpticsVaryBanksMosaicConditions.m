function run_OpticsVaryBanksMosaicConditions
% This is the script used to assess the impact of different optics models on the CSF
% 
    examinedOpticsModels = {...
        'Geisler' ...
        'ThibosBestPSFSubject3MMPupil' ...
        'ThibosDefaultSubject3MMPupil' ...
        'ThibosAverageSubject3MMPupil' ...
        'ThibosDefocusedSubject3MMPupil' ...
        'ThibosVeryDefocusedSubject3MMPupil' ...
    };
    
    % Optics to run
    params.opticsModel = examinedOpticsModels{1};

    params.pupilDiamMm = 2.0;   % What was used in Banks et al 87
   % params.pupilDiamMm = 3.0;   % 3 is more appropriate for a 100 cd/m2 mean scene luminance
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = true;
    params.visualizeSpatialScheme = true;
    params.visualizeOIsequence = true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % How to split the computation
    %computationInstance = 0;        % ALL mosaics
    computationInstance = 1;        % LARGEST mosaic
    %computationInstance = 2;        % ALL except the LARGEST mosaic
    params = getFixedParamsForOpticsImpactExperiment(params,computationInstance);

    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);

end

function params = getFixedParamsForOpticsImpactExperiment(params, computationInstance)
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
    mosaicParams = getParamsForMosaicWithLabel('originalBanks');
    
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
    params.luminancesExamined =  [34];
    
    % 'isomerizations', 'photocurrents'
    params.performanceSignal = 'isomerizations';
    
    %'mlpt'% 'svmV1FilterBank';
    params.performanceClassifier = 'mlpt';
    
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
        params.ramPercentageEmployed = 1.2;  % use all the RAM
        params.cyclesPerDegreeExamined =  [2 4 8 16 32 60]; % [3 5 10 20 30 60];
    elseif (computationInstance  == 1)
        % Largest mosaic in session 1 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 1;  % use 90% of the RAM
        params.cyclesPerDegreeExamined =  [2];
    elseif (computationInstance  == 2)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 1;  % use 1/2 the RAM
        params.cyclesPerDegreeExamined =  [4 8 16 32 60];
    else
        error('computational instance must be 0, 1 or 2');
    end
end

