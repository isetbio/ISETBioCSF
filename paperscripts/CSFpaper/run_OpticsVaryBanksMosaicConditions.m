function run_OpticsVaryBanksMosaicConditions
% This is the script used to assess the impact of different optics models on the CSF
% 
    examinedOpticsModels = {...
        'Geisler' ...
        'WvfHuman' ...
        'WvfHumanMeanOTFmagMeanOTFphase' ...
        'WvfHumanSubject1' ...
        'WvfHumanSubject2' ...
        'WvfHumanSubject3' ...
        'WvfHumanSubject4' ...
        'WvfHumanSubject5' ...
    };
    
    % Optics to run
    params.opticsModel = examinedOpticsModels{1};

    params.pupilDiamMm = 2.0;   % Default is 2 (as in Banks et al 87), but 3 is more appropriate for 100 cd/m2 monitor
    
    % Simulation steps to perform
    params.computeMosaic = true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % How to split the computation
    computationInstance = 0;
    params = getFixedParamsForOpticsImpactExperiment(params,computationInstance);

    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);

end

function params = getFixedParamsForOpticsImpactExperiment(params, computationInstance)
    params.blur = true;
    params.apertureBlur = true;
    
    params.imagePixels = 512;
    params.wavefrontSpatialSamples = 601;
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
    params.lowContrast = 0.00001;
    params.highContrast =  0.5;
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
        params.cyclesPerDegreeExamined =  [5 10 20 30 60]; % [3 5 10 20 30 60];
    elseif (computationInstance  == 1)
        % Largest mosaic in session 1 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.5;  % use 90% of the RAM
        params.cyclesPerDegreeExamined =  [2.5];
    elseif (computationInstance  == 2)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.5;  % use 1/2 the RAM
        params.cyclesPerDegreeExamined =  [5 10 20 40 60];
    else
        error('computational instance must be 0, 1 or 2');
    end
end

