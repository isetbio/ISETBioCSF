function run_PeformanceSignalVaryConditions
% This is the script used to assess the impact of different types of eye
% movements on the CSF  
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0 ; % 16 32  50 60;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection';
    
    % Optics to use
    opticsName = 'ThibosAverageSubject3MMPupil';
    %opticsName = 'ThibosDefaultSubject3MMPupil';
    
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
    params.emPathType = 'frozen0'; % 'random'; % 'frozen0'; % 
    params.centeredEMPaths = ~true;
    
    % Performance classifier
    params.performanceClassifier = 'svmV1FilterBank';
    %params.performanceClassifier = 'svmV1FilterEnsemble';
    
    if (strcmp(params.performanceClassifier,'svmV1FilterEnsemble'))
        
        ensembleFilterParams = struct(...
            'spatialPositionsNum',  9, ...
            'cyclesPerRFs', [1.0 1.5 2.0 2.5], ...
            'orientations', [0]);
        params.parforWorkersNumForClassification = 1 * params.parforWorkersNumForClassification;
        params.parforWorkersNumForClassification = 12;
        
%         ensembleFilterParams = struct(...
%             'spatialPositionsNum',  9, ...
%             'cyclesPerRFs', [2.5], ...
%             'orientations', [0]);
%         params.parforWorkersNumForClassification = 2 * params.parforWorkersNumForClassification;
%         
    end
    
    % Signals examined
    examinedSignals = {...
        'isomerizations' ...
        'photocurrents' ...
    };
    
    examinedSignalLabels = {...
        'photopigment excitations' ...
        'outer segment photocurrents' ...
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
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    examinedSignals = {examinedSignals{1:2}};
    
    % Go
  	for signalIndex = 1:numel(examinedSignals)
        params.performanceSignal = examinedSignals{signalIndex};
        
        if (strcmp(params.performanceClassifier,'svmV1FilterEnsemble'))
            fNames = fieldnames(ensembleFilterParams);
            for fNameIndex = 1:numel(fNames)
                fName = fNames{fNameIndex};
                params.spatialPoolingKernelParams.(fName) = ensembleFilterParams.(fName);
            end
        end
        [~,~, theFigData{signalIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'PerformanceSignal';
        theRatioLims = [0.2 1];
        theRatioTicks = [0.2 0.3 0.4 0.5 0.7 1];
        generateFigureForPaper(theFigData, examinedSignalLabels, variedParamName, sprintf('%s_%s_%sEM',mosaicName, opticsName, params.emPathType), ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    
end
