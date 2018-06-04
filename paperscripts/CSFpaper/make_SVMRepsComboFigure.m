function make_SVMRepsComboFigure
% This is the script used to assess the generate a combo SVM-type reps effect figure
%  
    % Which spatial frequency to analyze
    thePanelLabel = ' C ';         % Label for the first panel
    computationInstance = 16;  %  4 (4 c/deg) 8 (8 c/deg), 16 (16 c/deg) or 32 (32 c/deg)
    performanceClassifiersVisualized = {'svm', 'svmV1FilterBank', 'mlpt'};     % Choose between 'svm' and 'svmV1FilterBank'
    performanceClassifiersVisualizedLegends = {...
        'SVM - PCA projection', ...
        'SVM - spatial pooling'...
        'MLPT'
      }; 
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
    
    
    fixedParamName = sprintf('%.0fCPD_%s_%s_vs_%s', computationInstance, params.emPathType, performanceClassifiersVisualized{1}, performanceClassifiersVisualized{2});
    
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
    for performanceClassifierIndex = 1:numel(performanceClassifiersVisualized)
        for condIndex = numel(examinedCond):-1:1
            if (computeResponses && coneIndex==numel(examinedCond))
                params.computeResponses = true;
            else
                params.computeResponses = false;
            end
            params.performanceClassifier = performanceClassifiersVisualized{performanceClassifierIndex};
            if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
                params.spatialPoolingKernelParams.type = poolingType;

                if (strcmp(params.spatialPoolingKernelParams.type, 'V1QuadraturePair'))
                    params.spatialPoolingKernelParams.activationFunction = 'energy';
                elseif (strcmp(poolingType, 'V1CosUnit'))
                    params.spatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
                else
                    error('Unknwon poolingType: ''%s''.', poolingType);
                end
            end
            if (strcmp(params.performanceClassifier, 'mlpt')) && (condIndex == numel(examinedCond))
                paramsMLPT = params;
                paramsMLPT.lowContrast = 0.00003;
                paramsMLPT.highContrast =  1.0;
                paramsMLPT.nContrastsPerDirection =  20;
                paramsMLPT.nTrainingSamples = 1024;
                paramsMLPT.performanceTrialsUsed = 1024;
                [~, thePsychometricFunctions, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(paramsMLPT);
                thePsychometricFunctions = thePsychometricFunctions{:};
                psychometricFunctionMLPT = thePsychometricFunctions{1};
            elseif(~strcmp(params.performanceClassifier, 'mlpt'))
                params.performanceTrialsUsed = examinedCond(condIndex).performanceTrialsUsed;
                [~, thePsychometricFunctions, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
                if (numel(thePsychometricFunctions) > 1)
                    error('There were more than 1 spatial frequency point\n');
                end
                thePsychometricFunctions =  thePsychometricFunctions{:};
                thePsychometricFunctions = thePsychometricFunctions{1};
                if (performanceClassifierIndex == 1)
                    targetSFPsychometricFunctions{condIndex} = thePsychometricFunctions;
                else
                    targetSFPsychometricFunctions2{condIndex} = thePsychometricFunctions;
                end
                theTrials(condIndex) = params.performanceTrialsUsed ;
            end
        end
    end
    
    if (makeSummaryFigure)
        if (computationInstance == 16)
            csLims = [30 1500];
        elseif (computationInstance == 8)
            csLims = [30 1500];
        else
            csLims = [30 1500];
        end
    
        generatePsychometricFunctionsPlot(targetSFPsychometricFunctions, targetSFPsychometricFunctions2, psychometricFunctionMLPT, csLims, theTrials, legendsForPsychometricFunctions, performanceClassifiersVisualizedLegends, thePanelLabel, fixedParamName);
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

function generatePsychometricFunctionsPlot(psychometricFunctions, psychometricFunctions2, psychometricFunctionMLPT, csLims,  theTrials, trialLegends, classifierLegends, thePanelLabel, fixedParamName)
    conditionsNum = numel(psychometricFunctions);
    
    hFig = figure(15); clf;
    [theAxes, theAxes2] = formatFigureForPaper(hFig, ...
        'figureType', 'PSYCHOMETRIC_FUNCTIONS_2_CLASSIFIERS');

    colors(1,:) = [0.5 0.5 0.5];
    if (conditionsNum>1)
        colors(2:conditionsNum,:) = brewermap(conditionsNum-1, 'Set1');
    end
    
    % The psychometric functions for classifier #1
    for cond = 1:numel(psychometricFunctions)
        psyF = psychometricFunctions{cond};
        theThresholds(cond) = psyF.thresholdContrast;
        hold(theAxes, 'on');
        plot(theAxes, psyF.x, psyF.y, 'ko-', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.5 0.5 0.5] + 0.5*(squeeze(colors(cond,:))), ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    end
    hold(theAxes, 'off');
    xlabel(theAxes, 'contrast', 'FontWeight', 'Bold');
    ylabel(theAxes, 'percent correct', 'FontWeight', 'Bold');
    hL = legend(theAxes, trialLegends);
    title(theAxes, classifierLegends{1});
    
    % Panel label in the first panel
    inGraphTextPos = [min(psyF.x)*0.6 0.95];
    inGraphText = text(theAxes, inGraphTextPos(1), inGraphTextPos(2), thePanelLabel);
    inGraphTextFontSize = [];
    
    
    % The psychometric functions for classifier #2
    for cond = 1:numel(psychometricFunctions2)
        psyF = psychometricFunctions2{cond};
        theThresholds2(cond) = psyF.thresholdContrast;
        hold(theAxes2, 'on');
        plot(theAxes2, psyF.x, psyF.y, 'ks-', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.5 0.5 0.5] + 0.5*(squeeze(colors(cond,:))), ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    end
    hold(theAxes2, 'off');
    xlabel(theAxes2, 'contrast', 'FontWeight', 'Bold');
    ylabel(theAxes2, 'percent correct', 'FontWeight', 'Bold');
    title(theAxes2, classifierLegends{2});
    
    set(theAxes, 'XLim', [0.0005 0.15], 'XTick', [0.01 0.03 0.1 0.3], ...
        'YLim', [0.4 1.01], 'XTick', [0.001 0.003 0.01 0.03 0.1 0.3]);
    
    set(theAxes2, 'XLim', [0.0005 0.15], 'XTick', [0.01 0.03 0.1 0.3], ...
        'YLim', [0.4 1.01], 'XTick', [0.001 0.003 0.01 0.03 0.1 0.3]);
    
    % Format figure
    formatFigureForPaper(hFig, ...
        'figureType', 'PSYCHOMETRIC_FUNCTIONS_2_CLASSIFIERS', ...
        'plotRatiosOfOtherConditionsToFirst', false, ...
        'theText', inGraphText, ...
        'theTextFontSize', inGraphTextFontSize, ...
        'theAxes', theAxes, ...
        'theAxes2', theAxes2, ...
        'theLegend', hL, ...
        'theTextFontSize', inGraphTextFontSize);
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    variedParamName = 'TrialsNumPsychometricFunction';
    fixedParamName = strrep(fixedParamName, '\mu', 'micro');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
    
    
    hFig = figure(16); clf;
    theAxes = formatFigureForPaper(hFig, ...
        'figureType', 'CSF', ...
        'plotRatiosOfOtherConditionsToFirst', false);
    
    theThresholdsMLPT = psychometricFunctionMLPT.thresholdContrast;
    
    colors = brewermap(2, 'Set1');
    theTrialsLims = [400 80000];
    
    % The change in the CSFs
    plot(theAxes, theTrials, 1./theThresholds, 'ko-', ...
        'MarkerSize', 10, 'MarkerFaceColor', squeeze(colors(1,:)), ...
        'Color', squeeze(colors(1,:)), ...
        'LineWidth', 1.5);
    hold on;
    plot(theAxes, theTrials, 1./theThresholds2, 'ks-', ...
        'MarkerSize', 10, 'MarkerFaceColor',  squeeze(colors(2,:)), ...
        'Color', squeeze(colors(2,:)), ...
        'LineWidth', 1.5);
    plot(theAxes, theTrialsLims, 1/theThresholdsMLPT*[1 1], 'k--',...
         'LineWidth', 1.5);
     
    xlabel(theAxes, 'trials', 'FontWeight', 'Bold');
    ylabel(theAxes, 'contrast sensitivity', 'FontWeight', 'Bold');
    
    hL = legend(theAxes, classifierLegends);
    inGraphText = text(theAxes, 450, 1200, thePanelLabel);
    formatFigureForPaper(hFig, ...
        'figureType', 'CSF', ...
        'theAxes', theAxes, ...
        'theLegend', hL, ...
        'theText', inGraphText);
    
    set(theAxes, 'XTick', [300 1000 3000 10000 30000 100000],  ...
        'XLim', theTrialsLims, 'YLim', csLims);
    xlabel(theAxes, 'number of trials', 'FontWeight', 'bold', 'FontSize', 20);
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    variedParamName = 'TrialsNumCSF';
    fixedParamName = strrep(fixedParamName, '\mu', 'micro');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);

end
