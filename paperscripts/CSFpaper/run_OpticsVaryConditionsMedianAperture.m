function run_OpticsVaryConditionsMedianAperture
% This is the script used to assess the impact of different optics models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
        
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    
    % Adjust any params we want to change from their default values
    % Optics models we tested
    examinedApertureSizes = {...
        'minimum' ...
        'mean' ...
        'mean aperture 2 deg' ...
    };
    examinedApertureSizeLegends = {...
        'minimum aperture (default, 2 um)' ...
        'mean aperture (for each mosaic)' ...
        'mean aperture (for the 2 deg mosaic)'
    };

    % We differentiate the conditions by the number of trials because there
    % is no appropriate param -> Dir
    examinedTrialsNum = [
        1024
        1028
        1020
    ];
    
    idx = [1:3];
    examinedApertureSizes = {examinedApertureSizes{idx}};
    examinedApertureSizeLegends = {examinedApertureSizeLegends{idx}};
    examinedTrialsNum = examinedTrialsNum(idx);
    
    % Optics model to use
    params.opticsModel = 'ThibosBestPSFSubject3MMPupil';
    params.opticsModel = 'ThibosAverageSubject3MMPupil';

    % Spatial frequencies to use
    %params.cyclesPerDegreeExamined = [4 8 16 32 50 60];
    
    params.coneContrastDirection = 'L+M+S';
    
    params.cyclesPerDegreeExamined = [4 8 16 32 50 60];
    
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
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
        
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for modelIndex = 1:numel(examinedApertureSizes)
        params.nTrainingSamples = examinedTrialsNum(modelIndex); 
        [~,~, theFigData{modelIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'ApertureSize';
        theRatioLims = [0.4 10];
        theRatioTicks = [0.5 1 2 5 10];
%         generateFigureForPaper(theFigData, examinedOpticsModelLegends, variedParamName, mosaicName, ...
%             'figureType', 'CSF', ...
%             'inGraphText', ' A ', ...
%             'plotFirstConditionInGray', true, ...
%             'plotRatiosOfOtherConditionsToFirst', true, ...
%             'theRatioLims', theRatioLims, ...
%             'theRatioTicks', theRatioTicks ...
%             );
        generateFigureForPaper(theFigData, examinedApertureSizeLegends, variedParamName, mosaicName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' G ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    
end
