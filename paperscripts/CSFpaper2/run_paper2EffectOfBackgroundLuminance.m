function run_paper2EffectOfBackgroundLuminance
% Compute and contrast performance at the level of photocurrents for
% different background luminance levels
%
% Syntax:
%   run_paper2EffectOfBackgroundLuminance
%
% Description:
%    Compute and contrast performance at the different background luminance
%    levels in the absence of eye movements using a linear inference engine
%
%    The computation is done via the ecc-based cone efficiency & macular pigment
%    mosaic and the default Thibos subject. We use a 5 msec integration
%    time because there are no fixational eye movements. 
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
%    01/08/19  npc  Wrote it.

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Whether to compute responses
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = ~true;
    visualizePerformance = true;
    
    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamcs of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % How long the stimulus is. We might be changing the duration. 100 ms is the default
    stimulusDurationInSeconds = 100/1000;
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    % Compute photocurrent responses
    computePhotocurrents = true;
        
    nTrainingSamples = 1016;
    
    performanceSignal = 'photocurrents'; %'isomerizations'; % 'photocurrents';
    performanceClassifier = 'svmV1FilterBank';
    spatialPoolingKernelType = 'V1CosUnit';            % choose between 'V1CosUnit' and 'V1QuadraturePair';
    spatialPoolingKernelActivationFunction = 'linear'; % choose between 'linear' and 'energy';
    emPathType = 'frozen0';                            % choose between 'frozen0' and 'randomNoSaccades';
    
    % Assemble conditions list to be examined
    % Init condition index
    condIndex = 0;
    
       
    if (~computeResponses) && (~findPerformance)
        condIndex = condIndex+1;
        examinedCond(condIndex).label = '34 cd/m2';
        examinedCond(condIndex).backgroundLuminance = 34;
    end
    
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = '3.4 cd/m2';
    examinedCond(condIndex).backgroundLuminance = 3.4;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = '100 cd/m2';
    examinedCond(condIndex).backgroundLuminance = 100;
     
    condIndex = condIndex+1;
    examinedCond(condIndex).label = '340 cd/m2';
    examinedCond(condIndex).backgroundLuminance = 340;
%     
    use3400 = false;
    if (use3400)
        condIndex = condIndex+1;
        examinedCond(condIndex).label = '3400 cd/m2';
        examinedCond(condIndex).backgroundLuminance = 3400;
    end
    
    % Go
    examinedLegends = {};
    for condIndex = 1:numel(examinedCond)
        
        cond = examinedCond(condIndex);

        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds,  frameRate, stimulusDurationInSeconds, computationInstance);
        
        params.luminancesExamined = [cond.backgroundLuminance];
        
        % Use 1028 vs 1024 trials to differentiate results from old photocurrent computation
        params.nTrainingSamples = nTrainingSamples;
        
        % Try out for subset of SFs
        params.cyclesPerDegreeExamined = [8 16 24 32 50 60]; % [8 16 24 32 50 60]; % [16 24 32 50 60];
        
       
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
        params.performanceClassifier = performanceClassifier;
        params.performanceSignal = performanceSignal;
        params.emPathType = emPathType;
        
        
        params.spatialPoolingKernelParams.type = spatialPoolingKernelType;
        params.spatialPoolingKernelParams.activationFunction = spatialPoolingKernelActivationFunction;
        
        params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance);  
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
        variedParamName = sprintf('%sBkgndLuminance',performanceSignal);
        
        theRatioLims = [0.09 17];
        theRatioTicks = [0.1 0.3 1.0 3 10];
        theRatioLims = [0.1 5];
        theRatioTicks = [0.1 0.2 0.5 1 2 5];
        formatLabel = 'ComparedToBanksSubjects';
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'showSubjectData', ~true, ...
            'showSubjectMeanData', ~true, ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ... 
            'theLegendPosition', [0.57 0.81 0.41 0.14],  ...   % custom legend position and size
            'paperDir', 'CSFpaper2', ...                        % sub-directory where figure will be exported
            'figureHasFinalSize', true ...                      % publication-ready size
            );
        
        
    end
end

function params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, visualizeResponses, findPerformance, visualizePerformance)
                         
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
    params.computePhotocurrentResponseInstances = computePhotocurrents && computeResponses;
    params.visualizeResponses = visualizeResponses;
    params.visualizeResponsesWithSpatialPoolingSchemeInVideo = ~true;
    params.visualizeOuterSegmentFilters = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = visualizePerformance;
    params.deleteResponseInstances = ~true;
end
