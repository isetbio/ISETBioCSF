function run_EMVaryEccLMSMosaicConditions
% This is the script used to assess the impact of different EM types on the CSF
% 
    % 'random'; 'randomNoSaccades'; 'frozen0';
    params.emPathType = 'random'; % 'randomNoSaccades'; 'frozen0';    
    params.centeredEMPaths = false;
    
    %'mlpt'% 'svmV1FilterBank';
    params.performanceClassifier = 'mlpt';
    params.performanceClassifier = 'svm'; % 'svmV1FilterBank';
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % How to split the computation
    %computationInstance = 0;        % ALL mosaics
    %computationInstance = 1;        % LARGEST mosaic
    %computationInstance = 2;        % second LARGEST mosaic
    %computationInstance = 3;        % ALL except the two LARGEST mosaic
    computationInstance = 4;         % Just for testing
     
    params = getFixedParamsForOpticsImpactExperiment(params,computationInstance);

    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);

end

