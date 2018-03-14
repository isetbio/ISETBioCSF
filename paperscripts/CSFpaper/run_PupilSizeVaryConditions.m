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
        '2 mm pupil' ...
        '3 mm pupil' ...
    };


    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;

    if (strcmp(mosaicName, 'originalBanks')) && (strcmp(opticsName,'Geisler'))
        params.cyclesPerDegreeExamined = params.cyclesPerDegreeExamined(1:end-1);
    end
        
    % Go
    for pupilSizeIndex = 1:numel(examinedPupilSizes)
        params.pupilDiamMm = examinedPupilSizes{pupilSizeIndex};
        params.opticsModel = opticsName;
        
        [~,~, theFigData{pupilSizeIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'PupilSize';
        generateFigureForPaper(theFigData, examinedPupilSizeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), 'figureType', 'CSF');
    end
end