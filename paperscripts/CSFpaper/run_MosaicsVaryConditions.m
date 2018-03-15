function run_MosaicsVaryConditions
% This is the script used to assess the impact of different mosaic models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    examinedMosaicModels = {...
        'originalBanks' ...
        'ISETbioHexEccBasedNoScones' ...
        'ISETbioHexEccBasedLMSrealistic' ...
        };
 
    examinedMosaicLegends = {...
        'constant LM density (Banks ''87)' ...
        'eccentricity-based LM density' ...
        'eccentricity-based LMS density' ...
    };

    % Tun the mosaic-vary condition using the Geisler optics
    opticsName = 'Geisler';
      
    % Go !
    for modelIndex = 1:numel(examinedMosaicModels)
        mosaicName = examinedMosaicModels{modelIndex};
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        
        % Special case: the Geisler optics/original Banks mosaic was run up
        % to 50 c/deg
        if (strcmp(mosaicName, 'originalBanks')) && (strcmp(opticsName,'Geisler'))
            params.cyclesPerDegreeExamined = params.cyclesPerDegreeExamined(1:end-1);
        end
        
        params.opticsModel = opticsName;
       
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

        [~,~, theFigData{modelIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Mosaic';
        theRatioLims = [0.3 1];
        theRatioTicks = [0.3 0.5 0.7 1.0];
        generateFigureForPaper(theFigData, examinedMosaicLegends, variedParamName, opticsName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' Geisler optics ', ...
            'inGraphTextPos', [1.1 3], ...
            'inGraphTextFontSize', 18, ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

