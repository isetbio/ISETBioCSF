function run_MacularPigmentVaryConditions
% This is the script used to assess the impact of having an ecc-varying 
% macular pigment  on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    makeMosaicsFigure = ~true;
    
    % Mosaic to use
    examinedMosaicModels = {...
        'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection' ...
        'ISETbioHexEccBasedLMSrealisticEfficiencyCorrectionAndMacularPigment' ...
        };
    
    examinedMosaicLegends = {...
        'constant density MP' ...
        'ecc-based density MP'
    };


    idx = [1:2];
    examinedMosaicModels = {examinedMosaicModels{idx}};
     
    % Tun the mosaic-vary condition using the Geisler optics
    opticsName = 'Geisler';
      
    theMosaicTypesAtSpecificSF = {};
    
    % Go !
    for mosaicIndex = 1:numel(examinedMosaicModels)
        mosaicName = examinedMosaicModels{mosaicIndex};
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        
        params.opticsModel = opticsName;
       
        params.coneContrastDirection = 'L+M+S';
        %params.cyclesPerDegreeExamined = [2]; %[4 8 16 32 50];
    
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
        params.visualizeMosaic = makeMosaicsFigure;
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
        params.visualizePerformance = makeSummaryFigure;
        params.deleteResponseInstances = ~true;

        [theMosaicTypes, thePsychometricFunctions{mosaicIndex}, theFigData{mosaicIndex}] = ...
            run_BanksPhotocurrentEyeMovementConditions(params);
        
        if (makeMosaicsFigure)
            sfIndex = 1;
            
            theCurrentMosaic = theMosaicTypes.theMosaics{sfIndex};
            theMosaicTypesAtSpecificSF{numel(theMosaicTypesAtSpecificSF) + 1} = theCurrentMosaic;
            theCurrentMosaic.displayInfo();
        end
    end
    
    if (makeSummaryFigure)
        variedParamName = 'MacularPigment';
        theRatioLims = [0.3 1.2];
        theRatioTicks = [0.3 0.5 0.7 1.0 2.0];
        generateFigureForPaper(theFigData, examinedMosaicLegends, variedParamName, opticsName, ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

