function run_PupilSizeVaryConditions
% This is the script used to assess the impact of different pupil sizes on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'originalBanks'; 
    
    % Optics to use
    opticsName = 'Geisler';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    % Optics models we tested
    examinedPupilSizes = {...
        2 ...
        3 ...
    };
    examinedPupilSizeLegends = {...
        '2 mm pupil (34 cd/m^2)' ...
        '3 mm pupil (34 cd/m^2)' ...
    };

    idx = 1:2;
    examinedPupilSizes = examinedPupilSizes(idx);
    examinedPupilSizeLegends = examinedPupilSizeLegends(idx);
    
    % Stimulus cone contrast modulation vector
    params.opticsModel = opticsName;
    params.coneContrastDirection = 'L+M+S';
    params.cyclesPerDegreeExamined = [2 4 8 16 32 50];
    
    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeDisplay= ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
        
    % Go
    for pupilSizeIndex = 1:numel(examinedPupilSizes)
        params.pupilDiamMm = examinedPupilSizes{pupilSizeIndex};
        [~,~, theFigData{pupilSizeIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'PupilSize';
        theRatioLims = [0.5 3];
        theRatioTicks = logspace(log10(0.5), log10(2.6), 4);
        
        generateFigureForPaper(theFigData, examinedPupilSizeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' B ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end