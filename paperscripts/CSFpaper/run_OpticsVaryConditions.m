function run_OpticsVaryConditions
% This is the script used to assess the impact of different optics models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealistic';  %'originalBanks';
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    % Optics models we tested
    examinedOpticsModels = {...
        'Geisler' ...
        'ThibosBestPSFSubject3MMPupil' ...
        'ThibosDefaultSubject3MMPupil' ...
        'ThibosAverageSubject3MMPupil' ...
        'ThibosDefocusedSubject3MMPupil' ...
        'ThibosVeryDefocusedSubject3MMPupil' ...
    };
    examinedOpticsModelLegends = {...
        'Geisler PSF' ...
        'wvf-based PSF (best subject)' ...
        'wvf-based PSF (subject A)' ...
        'wvf-based PSF (subject B)' ...
        'wvf-based PSF (subject C)' ...
        'wvf-based PSF (subject D)' ...
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
    
    % Go
    for modelIndex = 1:numel(examinedOpticsModels)
        params.opticsModel = examinedOpticsModels{modelIndex};
        [~,~, theFigData{modelIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Optics';
        theRatioLims = [0.5 12];
        theRatioTicks = [0.5 1 2 5 10];
        generateFigureForPaper(theFigData, examinedOpticsModelLegends, variedParamName, mosaicName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

