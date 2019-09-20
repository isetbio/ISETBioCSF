function run_ConditionToVisualizeSpatialPoolingMechanisms
% Script for visualizing the spatial pooling mechanisms.
%
% Syntax:
%   run_ConditionToVisualizeSpatialPoolingMechanisms
%
% Description:
%   This is the script used to just visualize the different spatial pooling
%   mechanisms.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History
%    07/24/19  npc  Wrote it.

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    

    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use, 2.5 is better for capturing fixationalEM dynamics
    integrationTimeMilliseconds = 5;
    
    % Stimulus duration.
    stimulusDurationInSeconds = 100/1000;
    
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    % Get default params
    params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
    
    % Modify default params
    % (i) eye movement type
    params.emPathType = 'randomNoSaccades';     % 'random' (with saccades),  'randomNoSaccades', or 'frozen0'
    params.centeredEMPaths =  'atStimulusModulationOnset';
    
    % (ii) only 4 trials
    params.nTrainingSamples = 4;
    
    % (iii) Only the max contrast level
    params.lowContrast = 1.0; 
    params.highContrast =  1.0; 
    params.nContrastsPerDirection = 1;
    
    % (iv) only 20 c/deg
    params.cyclesPerDegreeExamined = [32]; % [24 32 50 60];
    
    % (v) minimum oi size (0.5 for visualization purposes)
    params.opticalImagePadSizeDegs = 0.5; 
    
    % (vi) * * * SPATIAL POOLING SELECTION  * * * 
    % Negative sign for the 'minimumMosaicFOVdegs' means use an ensemble of pooling mechanisms
    % on a mosaic whose FOV is the absolute value of the 'minimumMosaicFOVdegs'
    params.minimumMosaicFOVdegs = -0.328; %-0.492;
    params.spatialPoolingKernelParams.type = 'V1QuadraturePair';
    params.spatialPoolingKernelParams.activationFunction = 'energy';
    
    params.performanceClassifier = 'svmV1FilterEnsemble'; % 'svmV1FilterEnsemble'; %'svmV1FilterBank' or 'svmV1FilterEnsemble'
    
    if (strcmp(params.performanceClassifier, 'svmV1FilterEnsemble'))
        ensembleFilterParams = struct(...
                            'spatialPositionsNum',  1, ...   % 1 results in a 3x3 grid of spatial pooling templates
                            'spatialPositionOffsetDegs', 0.0328, ... 
                            'cyclesPerRFs', 6, ...           % each template contains 4 cycles of the stimulus
                            'orientations', 0);
        fNames = fieldnames(ensembleFilterParams);
        for fNameIndex = 1:numel(fNames)
            fName = fNames{fNameIndex};
            params.spatialPoolingKernelParams.(fName) = ensembleFilterParams.(fName);
       end                
    end
    
    
    % Simulation steps to perform
    params.computeResponses = ~true;
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computePhotocurrentResponseInstances = true;
    params.visualizeDisplay = ~true;
    params.visualizeResponses = true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = true;
    params.visualizeSpatialPoolingScheme = true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeOuterSegmentFilters = ~true;
    params.visualizeResponsesWithSpatialPoolingSchemeInVideo = ~true;
    
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);
end


