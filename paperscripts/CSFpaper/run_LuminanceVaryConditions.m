function run_LuminanceVaryConditions
% This is the script used to assess the impact of different pupil sizes on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 1;
    
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

    params.opticsModel = opticsName;
    params.pupilDiamMm = 2.0;
    
    % Stimulus cone contrast modulation vector
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
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;

 
    % Go
    for lumIndex = 1:numel(examinedLuminances)
        params.luminancesExamined = examinedLuminances{lumIndex};         
        [~,~, theFigData{lumIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Luminance';
        theRatioLims = [0.1 10];
        theRatioTicks = logspace(log10(0.1), log10(10), 5);
        
        load('BanksCSF.mat', 'BanksCSF');
        for k = 1:numel(theFigData)
            theFigData{k}.BanksCSF = BanksCSF;
        end
        
        generateFigureForPaper(theFigData, examinedPupilSizeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ...
            'showBanksPaperIOAcurves', true);
    end
end