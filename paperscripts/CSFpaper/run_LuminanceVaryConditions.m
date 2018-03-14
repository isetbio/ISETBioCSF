function run_LuminanceVaryConditions
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
    examinedLuminances = {...
        34 ...
        340 ...
        3.4 ...
    };
    examinedPupilSizeLegends = {...
        '34  cd/m^2 (2mm pupil)' ...
        '340 cd/m^2 (2mm pupil)' ...
        '3.4 cd/m^2 (2mm pupil)' ...
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

   
    params.opticsModel = opticsName;
    params.pupilDiamMm = 2.0;
    defaultCyclesPerDegreeExamined = [2 4 8 16 32 50];
    
    % Go
    for lumIndex = 1:numel(examinedLuminances)
        params.luminancesExamined = examinedLuminances{lumIndex};

        if (lumIndex == 2)
            params.cyclesPerDegreeExamined = defaultCyclesPerDegreeExamined(1:end-1);
        else
            params.cyclesPerDegreeExamined = defaultCyclesPerDegreeExamined;
        end
    
        [~,~, theFigData{lumIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Luminance';
        generateFigureForPaper(theFigData, examinedPupilSizeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'showBanksPaperIOAcurves', true);
    end
end