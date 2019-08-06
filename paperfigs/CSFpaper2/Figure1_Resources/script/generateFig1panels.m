function generateFig1panels()

    computationInstance = 16;
    pupilDiamMm = 3.0;
    integrationTimeMilliseconds = 5;
    stimulusDurationInSeconds = 100/1000;
    frameRate = 10; 
    params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
    params.responseExtinctionMilliseconds = 400;
    params.secondsToInclude = 0.1480*2;
    
    params.emPathType = 'randomNoSaccades';     % 'random' (with saccades),  'randomNoSaccades', or 'frozen0'
    % Only the max contrast level
    params.lowContrast = 1.0; 
    params.highContrast =  1.0; 
    params.nContrastsPerDirection = 1;
    
    %params.centeredEMPaths =  'atStimulusModulationOnset';
    %params.nTrainingSamples = 8;
    
    params.centeredEMPaths =  'atStimulusModulationMidPoint';
    params.nTrainingSamples = 4;
    
    % Simulation steps to perform
    params.computeResponses = ~true;
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computePhotocurrentResponseInstances = true;
    params.visualizeDisplay = ~true;
    params.visualizeResponses = true;
    params.visualizeResponsesWithSpatialPoolingSchemeInVideo = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeMosaicWithFirstEMpath = true;
    params.visualizeOuterSegmentFilters = true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = ~true;
    params.deleteResponseInstances = ~true;
    
    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);
end
