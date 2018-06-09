function run_Paper1FinalConditionsUsing2mmPupil
% This is the script used to assess how the final conditions in paper1
% with a 2 mm pupil compare to the Banks prediction.
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    

    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = ~true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    
    % Optics to use
    opticsName = 'ThibosAverageSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    % Pupil sizes compared
    examinedPupilSizes = {...
        2 ...
        3 ...
    };
    examinedPupilSizeLegends = {...
        '2mm pupil' ...
        '3mm pupil' ...
    };

    idx = 1:1;
    examinedPupilSizes = examinedPupilSizes(idx);
    examinedPupilSizeLegends = {examinedPupilSizeLegends{idx}};

    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];
    
    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;

    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;
    
    %params.emPathType = 'random';
    %params.centeredEMpaths = true;

    params.performanceClassifier = 'svmV1FilterBank';
    params.spatialPoolingKernelParams.type = 'V1CosUnit';
    params.spatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
                
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
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
    for pupilIndex = 1:numel(examinedPupilSizes)
        params.pupilDiamMm  = examinedPupilSizes(pupilIndex);
        [~,~, theFigData{pupilIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'PupilSize';
        theRatioLims = [0.05 0.5];
        theRatioTicks = [0.05  0.1 0.2 0.5];
        generateFigureForPaper(theFigData, examinedInferenceEngineLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

