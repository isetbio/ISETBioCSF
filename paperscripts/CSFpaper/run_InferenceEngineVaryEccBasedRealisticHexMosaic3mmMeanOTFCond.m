function run_InferenceEngineVaryEccBasedRealisticHexMosaic3mmMeanOTFCond
% This is the script used to assess the impact of the inference engine
% The optics used here is the 3mm pupil, meanOTF model, and the mosaic is
% the ecc-based hex mosaic with realistic S-cone positions

    examinedInferenceEngines = {...
        'mlpt' ...
        'svmV1FilterBank' ...
        };
    
    % Inference engine to run
    params.performanceClassifier = examinedInferenceEngines{1};
    
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
    params = getFixedParamsForInferenceEngineImpactExperiment(params,computationInstance);
    
    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);
    
    fprintf('Analyze mosaic to make sure we have the desired LMS ratios in the S-cone region');
end

function params = getFixedParamsForInferenceEngineImpactExperiment(params,computationInstance)
    % Optics for all subsequent computations
    params.opticsModel  = 'WvfHumanMeanOTFmagMeanOTFphase';
    params.pupilDiamMm  = 3.0;  
    params.blur         = true;
    params.apertureBlur = true;
    
    params.imagePixels = 1024;
    
    % 'random'; 'frozen0';
    params.emPathType = 'frozen0'; %random'; %'random';     
    params.centeredEMPaths = false;
    
     % Use a subset of the trials.
    params.nTrainingSamples = 1024;
    
    % Mosaic params
    mosaicParams = getParamsForMosaicWithLabel('ISETbioHexEccBasedLMSrealistic');
    fNames = fieldnames(mosaicParams);
    for k = 1:numel(fNames)
        params.(fNames{k}) = mosaicParams.(fNames{k});
    end
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

