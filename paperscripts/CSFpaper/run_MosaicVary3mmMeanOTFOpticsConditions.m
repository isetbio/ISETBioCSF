function run_MosaicVary3mmMeanOTFOpticsConditions
% This is the script used to assess the impact of the mosaic's spatial structure on the CSF
% The optics used here is the 3mm pupil, meanOTF model
%
    examinedMosaicModels = {...
        'originalBanksNoRotation' ...
        'ISETbioHexRegNoScones' ...
        'ISETbioHexRegLMS' ...
        'ISETbioHexEccBasedLMS' ...
        'ISETbioHexEccBasedLMSrealistic' ...
    };
    
    % Mosaic to run
    params = getParamsForMosaicWithLabel(examinedMosaicModels{4});
    
    % Simulation steps to perform
    params.computeMosaic = true; 
    params.visualizeMosaic = true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % How to split the computation
    computationInstance = 3;
    params = getFixedParamsForMosaicImpactExperiment(params,computationInstance);
    
    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);
end


function params = getFixedParamsForMosaicImpactExperiment(params, computationInstance)
    % Optics for all subsequent computations
    params.opticsModel  = 'WvfHumanMeanOTFmagMeanOTFphase';
    params.pupilDiamMm  = 3.0;          % 3mm pupil is more appropriate for 100 cd/m2 monitor
    params.blur         = true;
    params.apertureBlur = true;
    
    params.imagePixels = 1024;
    
    % 'random'; 'frozen0';
    params.emPathType = 'frozen0'; %random'; %'random';     
    params.centeredEMPaths = false;
   
    % Use a subset of the trials.
    params.nTrainingSamples = 1024;
   
    % Mosaic params
    params.integrationTimeMilliseconds =  5.0;
    
    % response params
    params.responseStabilizationMilliseconds = 10;
    params.responseExtinctionMilliseconds = 50;
    
    % Conditions 
    params.lowContrast = 0.0001;
    params.highContrast = 0.3; % 0.8;
    params.nContrastsPerDirection =  18; %20;
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
        params.ramPercentageEmployed = 1.0;  % use all the RAM
        params.cyclesPerDegreeExamined =  [2.5 5 10 20 40 50];
    elseif (computationInstance  == 1)
        % Largest mosaic in session 1 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.9;  % use 90% of the RAM
        params.cyclesPerDegreeExamined =  [5 10 20 40 50];
    elseif (computationInstance  == 2)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.5;  % use 1/2 the RAM
        params.cyclesPerDegreeExamined =  [5];
    elseif (computationInstance  == 3)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.5;  % use 1/2 the RAM
        params.cyclesPerDegreeExamined =  [20];
    else
        error('computational instance must be 0, 1 or 2');
    end
end

