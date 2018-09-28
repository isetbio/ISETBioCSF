function run_LuminanceVaryConditions
% This is the script used to assess the impact of different pupil sizes on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all except 2 largest)
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

    % Set this to true to only show the 34 cd/m2 data
    compareToSubjectDataFor34CDM2 = ~true;
    
    if (compareToSubjectDataFor34CDM2)
        idx = 1:1;
    else
        idx = 1:3;
    end
    examinedLuminances = {examinedLuminances{idx}};
    examinedPupilSizeLegends = {examinedPupilSizeLegends{idx}};
    
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
    
    params.computeResponses = ~true;
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
    params.findPerformance =~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;

 
    % Go
    for lumIndex = 1:numel(examinedLuminances)
        params.luminancesExamined = examinedLuminances{lumIndex};
        [~,~, theFigData{lumIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Luminance';
        if (compareToSubjectDataFor34CDM2)
            theRatioLims = [0.02 2.0];
            theRatioTicks = [0.05  0.1 0.2 0.5 1.0];
            showOnly23CDM2IOAcurve = true;
        else
            theRatioLims = [0.1 10];
            theRatioTicks = [0.1 0.3 1 3 10]; 
            showOnly23CDM2IOAcurve = false;
        end
        
        load('BanksCSF.mat', 'BanksCSF');
        for k = 1:numel(theFigData)
            theFigData{k}.BanksCSF = BanksCSF;
        end
        
        generateFigureForPaper(theFigData, examinedPupilSizeLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'showSubjectData', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ...
            'showBanksPaperIOAcurves', true, ...
            'showOnly23CDM2IOAcurve', showOnly23CDM2IOAcurve);
    end
end