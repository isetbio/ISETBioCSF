function run_SVMRepsVaryConditions
% This is the script used to assess the impact of different # of trials on the SVM-based CSF
%  
    % Which spatial frequency to analyze
    computationInstance = 8;  %  4 (4 c/deg) 8 (8 c/deg), 16 (16 c/deg) or 32 (32 c/deg)
    performanceClassifier = 'svm'; %';     % Choose between 'svm' and 'svmV1FilterBank'
    poolingType = 'V1CosUnit' ;        % Choose between 'V1CosUnit' and 'V1QuadraturePair (only applicable if classifier = 'svmV1FilterBank')
        
    computeResponses = ~true;
    findPerformance = ~true;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealistic'; 
    
    % Optics to use
    opticsName = 'ThibosAverageSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    
    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;

    params.mosaicRotationDegs = 0;
    
    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;
 
    % Contrast range to examine
    if (computationInstance == 32)
        params.lowContrast = 0.01;
        params.highContrast =  0.3;
        params.nContrastsPerDirection =  20;
    elseif (computationInstance == 16)
        params.lowContrast = 0.001;
        params.highContrast =  0.3;
        params.nContrastsPerDirection =  18;
    elseif (computationInstance == 8)
        params.lowContrast = 0.001;
        params.highContrast = 0.1;
        params.nContrastsPerDirection =  18;
    elseif (computationInstance == 4)
        params.lowContrast = 0.001;
        params.highContrast = 0.1;
        params.nContrastsPerDirection =  18;
    else
        error('computationInstance %d not implemented.', computationInstance);
    end
    
    % Trials to use in the classifier - vary this one 
    switch (computationInstance)
        case 4
            params.nTrainingSamples = 1024*8;
            trainingSamples = params.nTrainingSamples ./ (2.^[4 3 2 1]);
        case 8
            params.nTrainingSamples = 1024*32;
            trainingSamples = params.nTrainingSamples ./ (2.^[6 5 4 3 2 1]);
        case 16
            params.nTrainingSamples = 1024*32;
            trainingSamples = params.nTrainingSamples ./ (2.^[6 5 4 3 2 1 0]);
        case 32
            params.nTrainingSamples = 1024*64;
            trainingSamples = params.nTrainingSamples ./ (2.^[7 6 5 4 3 2 1 0]);
        otherwise
            error('Training samples sequence not set for computationInstance:%d', computationInstance);
    end
    
    for k = 1:numel(trainingSamples)
        performanceTrialsUsed = trainingSamples(k);
        examinedCond(k).performanceTrialsUsed = performanceTrialsUsed;
        legends{k} = sprintf('SVM, %d trials', performanceTrialsUsed);
        legendsForPsychometricFunctions{k} = sprintf('%d trials', performanceTrialsUsed);
    end
    
    theTitle = ''; %sprintf('%2.0f c/deg, %s\n%s', computationInstance, examinedCond(1).classifier, params.emPathType);
    fixedParamName = sprintf('%.0fCPD_%s_%s', computationInstance, params.emPathType, performanceClassifier);
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    
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
    params.findPerformance = findPerformance;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for condIndex = numel(examinedCond):-1:1
        if (computeResponses && coneIndex==numel(examinedCond))
            params.computeResponses = true;
        else
            params.computeResponses = false;
        end
        params.performanceClassifier = performanceClassifier;
        if (strcmp(performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = poolingType;
            
            if (strcmp(params.spatialPoolingKernelParams.type, 'V1QuadraturePair'))
                params.spatialPoolingKernelParams.activationFunction = 'energy';
            elseif (strcmp(poolingType, 'V1CosUnit'))
                params.spatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
            else
                error('Unknwon poolingType: ''%s''.', poolingType);
            end
        end
        params.performanceTrialsUsed = examinedCond(condIndex).performanceTrialsUsed;
        [~, thePsychometricFunctions, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
       
        if (numel(thePsychometricFunctions) > 1)
            error('There were more than 1 spatial frequency point\n');
        end
        thePsychometricFunctions =  thePsychometricFunctions{:};
        thePsychometricFunctions = thePsychometricFunctions{1};
        sf0PsychometricFunctions{condIndex} = thePsychometricFunctions;
        theTrials(condIndex) = params.performanceTrialsUsed ;
    end
    
    
    if (makeSummaryFigure)
        if (computationInstance == 16)
            csLims = [20 80];
        elseif (computationInstance == 8)
            csLims = [60 140];
        else
            csLims = [1 1000];
        end
    
        generatePsychometricFunctionsPlot(sf0PsychometricFunctions, csLims, theTrials, legendsForPsychometricFunctions, theTitle, fixedParamName);
        variedParamName = 'SVMTrials';
        theRatioLims = [0.05 0.5];
        theRatioTicks = [0.05 0.1 0.2 0.5];
        generateFigureForPaper(theFigData, legends, variedParamName, '', ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

function generatePsychometricFunctionsPlot(psychometricFunctions, csLims,  theTrials, trialLegends, theTitle, fixedParamName)
    conditionsNum = numel(psychometricFunctions);
    
    hFig = figure(15); clf;
    [theAxes, theRatioAxes] = formatFigureForPaper(hFig, ...
        'figureType', 'PSYCHOMETRIC_FUNCTIONS');

    colors(1,:) = [0.5 0.5 0.5];
    if (conditionsNum>1)
            colors(2:conditionsNum,:) = brewermap(conditionsNum-1, 'Set1');
    end
    
    for cond = 1:numel(psychometricFunctions)
        psyF = psychometricFunctions{cond};
        hold(theAxes, 'on');
        plot(theAxes, psyF.x, psyF.y, 'ko-', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.5 0.5 0.5] + 0.5*(squeeze(colors(cond,:))), ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    end
    
    for cond = 1:numel(psychometricFunctions)
        psyF = psychometricFunctions{cond};
        theThresholds(cond) = psyF.thresholdContrast;
%         plot(theAxes, psyF.xFit, psyF.yFit, '-', ...
%             'Color', squeeze(colors(cond,:)), 'LineWidth', 1.5);
        plot(theAxes, psyF.x, psyF.y, 'ko', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.5 0.5 0.5] + 0.5*(squeeze(colors(cond,:))));
    end
    hold(theAxes, 'off');
    xlabel(theAxes, 'contrast', 'FontWeight', 'Bold');
    ylabel(theAxes, 'percent correct', 'FontWeight', 'Bold');
    hL = legend(theAxes, trialLegends);
    
    plot(theRatioAxes, theTrials, 1./theThresholds, 'ko-', ...
        'MarkerSize', 10, 'MarkerFaceColor', [0.75 0.75 0.75], ...
        'LineWidth', 1.5);
    xlabel(theRatioAxes, 'trials', 'FontWeight', 'Bold');
    ylabel(theRatioAxes, 'contrast sensitivity', 'FontWeight', 'Bold');
    
    inGraphTextFontSize = [];
    % Format figure
    formatFigureForPaper(hFig, ...
        'figureType', 'PSYCHOMETRIC_FUNCTIONS', ...
        'plotRatiosOfOtherConditionsToFirst', false, ...
        'theAxes', theAxes, ...
        'theRatioAxes', theRatioAxes, ...
        'theLegend', hL, ...
        'theTextFontSize', inGraphTextFontSize);
    

    csTicks = [2 5 6 7 8 9 10  20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 500 1000 2000 5000 10000];
    

    if (~isempty(theTitle))
    theText = text(2000, csLims(2)*0.97, theTitle);
    set(theText, 'FontSize', 16, ...
                    'FontWeight', 'Normal', ...
                    'BackgroundColor', [1 1 1], ...
                    'EdgeColor', [ 0 0 0], ...
                    'LineWidth', 1.0);
    end
    
    set(theAxes, 'XLim', [0.001 0.35], 'XTick', [0.01 0.03 0.1 0.3], ...
        'YLim', [0.4 1.0], 'XTick', [0.003 0.01 0.03 0.1 0.3]);
    
    set(theRatioAxes, 'XTick', [300 1000 3000 10000 30000 100000],  ...
        'YScale', 'linear', 'XLim', [400 40000], ...
        'YLim', csLims, 'YTick', csTicks);
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    variedParamName = 'TrialsNum';
    fixedParamName = strrep(fixedParamName, '\mu', 'micro');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end
