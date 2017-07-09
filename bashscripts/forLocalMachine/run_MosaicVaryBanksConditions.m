function run_MosaicVaryBanksConditions
 
    examinedMosaicModels = {...
        'Banks87' ...
    };
    
    % Mosaic to run
    params = getParamsForMosaic(examinedMosaicModels{1});
    
    % Simulation steps to perform
    params.computeMosaic = true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % How to split the computation
    computationInstance = 2;
    params = getFixedParamsForMosaicImpactExperiment(params,computationInstance);
    
    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);
end

function params = getParamsForMosaic(mosaicName)
    switch mosaicName
        case 'Banks87'
            params.coneSpacingMicrons = 3.0;
            params.innerSegmentDiameter = 3.0;    % for a circular sensor
            params.conePacking = 'hexReg';
            params.LMSRatio = [0.67 0.33 0];
            params.mosaicRotationDegs = 30;
        otherwise
            error('Unknown mosaic: ''%s''.', mosaicName);
    end
end

function params = getFixedParamsForMosaicImpactExperiment(params, computationInstance)

    params.opticsModel = 'WvfHumanMeanOTFmagMeanOTFphase';
    params.blur = true;
    params.apertureBlur = true;
    
    params.imagePixels = 1024;
    
    % 'random'; 'frozen0';
    params.emPathType = 'frozen0'; %random'; %'random';     
    params.centeredEMPaths = false;
   
    % Use a subset of the trials. Specify [] to use all available trials
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
    
    % Freeze noise for repeatable results
    params.freezeNoise = true;
    
    % Split computations and specify RAM memory
    if (computationInstance == 0)
        % All conditions in 1 MATLAB session
        params.ramPercentageEmployed = 1.0;  % use all the RAM
        params.cyclesPerDegreeExamined =  [2.5 5 10 20 40 50];
    elseif (computationInstance  == 1)
        % Largest mosaic in session 1 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.5;  % use 90% of the RAM
        params.cyclesPerDegreeExamined =  [2.5];
    elseif (computationInstance  == 2)
        % Remainin mosaics in session 2 of 2 parallel MATLAB sessions
        params.ramPercentageEmployed = 0.5;  % use 1/2 the RAM
        params.cyclesPerDegreeExamined =  [5 10 20 40 50];
    else
        error('computational instance must be 0, 1 or 2');
    end
end

