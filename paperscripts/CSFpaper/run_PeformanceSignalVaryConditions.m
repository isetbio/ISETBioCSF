function run_PeformanceSignalVaryConditions
% This is the script used to assess the impact of different types of eye movements on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection';
    
    % Optics to use
    opticsName = 'ThibosAverageSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
    
    % All conds  with 2 mm pupil to compare to Banks subject data
    params.pupilDiamMm = 2.0;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    
    % Response duration params
    params.frameRate = 20; %(20 frames/sec, so 2 frames, each 50 msec long)
    params.responseStabilizationMilliseconds = 100;
    params.responseExtinctionMilliseconds = 50;
    
    % Eye movement setup
    params.emPathType = 'random';
    params.centeredEMPaths = ~true;
    
    % Performance classifier
    params.performanceClassifier = 'svmV1FilterBank';
    
    examinedSignals = {...
        'isomerizations' ...
        'photocurrents' ...
    };
    
    examinedSignalLabels = {...
        'isomerizations' ...
        'photocurrents' ...
    };
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    
    params.visualizeMosaic = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
  	for signalIndex = 1:numel(examinedSignals)
        params.performanceSignal = examinedSignals{signalIndex};
        [~,~, theFigData{signalIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'PerformanceSignal';
        theRatioLims = [0.2 1];
        theRatioTicks = [0.2 0.3 0.4 0.5 0.7 1];
        generateFigureForPaper(theFigData, examinedSignalLabels, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    
end
