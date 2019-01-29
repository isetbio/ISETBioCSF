function run_ConditionToVisualizePhotocurrentAndEyeMovements
% Script for visualizing photocurrent and eye movements.
%
% Syntax:
%   run_ConditionToVisualizePhotocurrentAndEyeMovements
%
% Description:
%    
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
%    01/28/19  npc  Wrote it.

% This is the script used to just visualize responses

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 16;
    

    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use, 2.5 is better for capturing fixationalEM dynamics
    integrationTimeMilliseconds = 5.0;
    
    % How long the stimulus is.
    % We might be changing the duration. 100 ms is the default
    stimulusDurationInSeconds = 100/1000;
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    
    % Get params
    params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
       
    % Modify default params
    % Eye movement types
    params.emPathType = 'randomNoSaccades';  % 'random' (with saccades),  'randomNoSaccades', or 'none'
    
    % Eye movement centering
    % 'atStimulusModulationMidPoint' (the centroid of the emPath within the stimulation time is at (0,0)), OR 
    % 'atStimulusModulationOnset' (he em position is (0,0) at stimulus onset)
    
    params.centeredEMPaths = 'atStimulusModulationMidPoint'; % 'atStimulusModulationOnset'
    
    % Only the max contrast level
    params.lowContrast = 1.0; 
    params.highContrast =  1.0; 
    params.nContrastsPerDirection = 1;
    params.nTrainingSamples = 1024;
        
    % Simulation steps to perform
    params.computeResponses = true;
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computePhotocurrentResponseInstances = true;
    params.visualizeDisplay = ~true;
    params.visualizeResponses = true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeMosaicWithFirstEMpath = true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = ~true;
    params.deleteResponseInstances = ~true;
    
    % Go
    run_BanksPhotocurrentEyeMovementConditions(params);