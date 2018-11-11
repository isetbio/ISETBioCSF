function make_SVMRepsComboFigure
% This is the script used to assess the generate panels of Figures 7 and 8
%  
    % Which spatial frequency to analyze
    thePanelLabels = {' B ', ' C '};         % Label for the two psychometric function panels
    thePanelLabels = {'', ''}; 
    
    computationInstance = 32;  %  4 (4 c/deg) 8 (8 c/deg), 16 (16 c/deg) or 32 (32 c/deg)
    
    
    performanceClassifiersVisualized = {'svm', 'svmV1FilterBank', 'mlpt'};     % Choose between 'svm' and 'svmV1FilterBank'
    performanceClassifiersVisualizedLegends = {...
        'SVM-PCA', ...
        'SVM-Template'...
        'ideal observer'
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
    params.mosaicRotationDegs = 0;
    
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
            computationalInstanceLabel = '';
        case 8
            params.nTrainingSamples = 1024*32;
            trainingSamples = params.nTrainingSamples ./ (2.^[6 5 4 3 2 1]);
            computationalInstanceLabel = ' A ';
        case 16
            params.nTrainingSamples = 1024*32;
            trainingSamples = params.nTrainingSamples ./ (2.^[6 5 4 3 2 1 0]);
            computationalInstanceLabel = ' B ';
        case 32
            params.nTrainingSamples = 1024*64;
            trainingSamples = params.nTrainingSamples ./ (2.^[7 6 5 4 3 2 1 0]);
            computationalInstanceLabel = ' C ';
        otherwise
            error('Training samples sequence not set for computationInstance:%d', computationInstance);
    end
    
    computationalInstanceLabel = '';
    
    for k = 1:numel(trainingSamples)
        performanceTrialsUsed = trainingSamples(k);
        examinedCond(k).performanceTrialsUsed = performanceTrialsUsed;
        legends{k} = sprintf('SVM, %d trials', performanceTrialsUsed);
        legendsForPsychometricFunctions{k} = sprintf('N = %d', performanceTrialsUsed);
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
        if (computationInstance == 32)
            csLims = [5 800];
            contrastLims =  [0.01 0.3];
            ratioLims = [.1 .7];
            showLegend = true;
        elseif (computationInstance == 16)
            csLims = [5 800];
            contrastLims =  [0.01 0.3]/3;
            ratioLims = [.1 .7];
            showLegend = true;
        elseif (computationInstance == 8)
            csLims = [5 800];
            contrastLims =  [0.01 0.3]/6;
            ratioLims = [.1 .7];
            showLegend = true;
        else
            contrastLims =  [0.0005 0.15];
            csLims = [30 600];
            ratioLims = [.1 .7];
            showLegend = true;
        end
    
        theData = generatePsychometricFunctionsPlot(targetSFPsychometricFunctions, ...
            targetSFPsychometricFunctions2, psychometricFunctionMLPT, ...
            computationalInstanceLabel, ...
            csLims, contrastLims, ratioLims, theTrials, legendsForPsychometricFunctions, ...
            performanceClassifiersVisualizedLegends, thePanelLabels, fixedParamName, showLegend);
        
        variedParamName = 'SVMTrials';
        theRatioLims = [0.05 0.5];
        theRatioTicks = [0.05 0.1 0.2 0.5];
        generateFigureForPaper(theFigData, legends, variedParamName, '', ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks...
            );
        
        save(sprintf('theData%d.mat',computationInstance), 'theData');
    end
end

function theData = generatePsychometricFunctionsPlot(psychometricFunctions, ...
    psychometricFunctions2, psychometricFunctionMLPT, ...
    computationalInstanceLabel, ...
    csLims,  contrastLims, ratioLims, theTrials, trialLegends, ...
    classifierLegends, thePanelLabels, fixedParamName, showLegend)


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
        plot(theAxes, psyF.x, psyF.y, 'ko-', 'MarkerSize', 12, ...
            'MarkerFaceColor', squeeze(colors(cond,:)).^0.5, ...
            'MarkerEdgeColor', squeeze(colors(cond,:))*0.5, ...
            'Color', [0 0 0], 'LineWidth', 1.5);
    end
    hold(theAxes, 'off');
    xlabel(theAxes, '\it contrast', 'FontWeight', 'Normal');
    ylabel(theAxes, '\it percent correct', 'FontWeight', 'Normal');
    if (showLegend)
        hL = legend(theAxes, trialLegends, 'Interpreter', 'tex');
    else
        hL = [];
    end
    %title(theAxes, sprintf('%s (%d c/deg)', classifierLegends{1}, computationInstance));
    
    % Panel label in the first panel
%     inGraphTextPos = [min(psyF.x)*1.07 0.95];
%     inGraphText = text(theAxes, inGraphTextPos(1), inGraphTextPos(2), ' A ');
%     
    
    % The psychometric functions for classifier #2
    for cond = 1:numel(psychometricFunctions2)
        psyF = psychometricFunctions2{cond};
        theThresholds2(cond) = psyF.thresholdContrast;
        hold(theAxes2, 'on');
        plot(theAxes2, psyF.x, psyF.y, 'ko-', 'MarkerSize', 12, ...
            'MarkerFaceColor', squeeze(colors(cond,:)).^0.5, ...
            'MarkerEdgeColor', squeeze(colors(cond,:))*0.5, ...
            'Color', [0 0 0], 'LineWidth', 1.5);
    end
    hold(theAxes2, 'off');
    xlabel(theAxes2, '\it contrast', 'FontWeight', 'Normal');
    %ylabel(theAxes2, 'percent correct', 'FontWeight', 'Bold');
    
    if (showLegend)
        hL2 = legend(theAxes2, trialLegends);
    else
        hL2 = [];
    end
    %title(theAxes2, sprintf('%s (%d c/deg)', classifierLegends{2}, computationInstance));
    

    t = text(theAxes, contrastLims(1)*1.1, 0.955, thePanelLabels{1});
    t2 = text(theAxes2, contrastLims(1)*1.1, 0.955, thePanelLabels{2});
     
    % Format figure
    formatFigureForPaper(hFig, ...
        'figureType', 'PSYCHOMETRIC_FUNCTIONS_2_CLASSIFIERS', ...
        'plotRatiosOfOtherConditionsToFirst', false, ...
        'theAxes', theAxes, ...
        'theAxes2', theAxes2, ...
        'theText', t, ...
        'theText2', t2, ...
        'theLegend', hL, ...
        'theLegend2', hL2);
    
    set(theAxes, 'XLim', contrastLims, ...
        'YLim', [0.4 1.01], 'XTick', [0.001 0.003 0.01 0.03 0.1 0.3], ...
        'XTickLabel', {'0.001', '0.003', '0.01', '0.03', '0.1', '0.3'});
    
    set(theAxes2, 'XLim', contrastLims, ...
        'YLim', [0.4 1.01], 'XTick', [0.001 0.003 0.01 0.03 0.1 0.3], ...
        'XTickLabel', {'0.001', '0.003', '0.01', '0.03', '0.1', '0.3'});
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    variedParamName = 'TrialsNumPsychometricFunction';
    fixedParamName = strrep(fixedParamName, '\mu', 'micro');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
    
    
    theThresholdsMLPT = psychometricFunctionMLPT.thresholdContrast;
    
    theData{1} = struct('trialsNum', theTrials, 'contrastSensitivity', 1/theThresholdsMLPT * ones(size(theTrials)));
    theData{2} = struct('trialsNum', theTrials, 'contrastSensitivity', 1./theThresholds);
    theData{3} = struct('trialsNum', theTrials, 'contrastSensitivity', 1./theThresholds2);
    
    classifierLegends = {'ideal observer', classifierLegends{1}, classifierLegends{2}};
    theTrialsLims = [400 80000];
     
    hFig = generateCSFTrialsFigure(theData, classifierLegends, ...
       'plotRatiosOfOtherConditionsToFirst', true, ...
       'theTrialsLims', theTrialsLims, ...
       'theCSLims', csLims, ...
       'inGraphText', computationalInstanceLabel, ...
       'inGraphTextPos', [480 620], ...
       'theRatioLims', ratioLims, ...
       'theRatioTicks', [0.1 0.2 0.3 0.4 0.6], ...
       'showLegend', showLegend...
    );
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    variedParamName = 'TrialsNumCSF';
    fixedParamName = strrep(fixedParamName, '\mu', 'micro');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);

end



function hFig = generateCSFTrialsFigure(theData, variedParamLegends, varargin)

    p = inputParser;
    p.addParameter('plotFirstConditionInGray', true, @islogical);
    p.addParameter('plotRatiosOfOtherConditionsToFirst', false, @islogical);
    p.addParameter('theRatioLims', []);
    p.addParameter('theRatioTicks', []);
    p.addParameter('theTrialsLims', []);
    p.addParameter('theCSLims', []);
    p.addParameter('inGraphText', '', @ischar);
    p.addParameter('inGraphTextPos', [], @isnumeric);
    p.addParameter('inGraphTextFontSize', [], @isnumeric);
    p.addParameter('showLegend', true, @islogical);
    p.parse(varargin{:});
    
  
    plotFirstConditionInGray = p.Results.plotFirstConditionInGray;
    inGraphText = p.Results.inGraphText;
    inGraphTextPos = p.Results.inGraphTextPos;
    inGraphTextFontSize = p.Results.inGraphTextFontSize;
    plotRatiosOfOtherConditionsToFirst = p.Results.plotRatiosOfOtherConditionsToFirst;
    theRatioLims = p.Results.theRatioLims;
    theRatioTicks = p.Results.theRatioTicks;
    theTrialsLims = p.Results.theTrialsLims;
    theCSLims = p.Results.theCSLims;
    showLegend = p.Results.showLegend;
    
    % Initialize figure
    hFig = figure(20); clf;
    [theAxes, theRatioAxes] = formatFigureForPaper(hFig, ...
        'figureType', 'CSF', ...
        'plotRatiosOfOtherConditionsToFirst', plotRatiosOfOtherConditionsToFirst);
    
    
    % Displayed colors
    if (plotFirstConditionInGray)
        colors(1,:) = [0.5 0.5 0.5];
        if (numel(theData)>1)
            colors(2:numel(theData),:) = brewermap(numel(theData)-1, 'Set1');
        end
    else
        colors = brewermap(numel(theData), 'Set1');
    end
    
    faceColor = [0.3 0.3 0.3];
    
    hold(theAxes(1), 'on');
    for condIndex = 1:numel(theData)
            plot(theAxes, theData{condIndex}.trialsNum, theData{condIndex}.contrastSensitivity, 'o-', 'Color', squeeze(colors(condIndex,:)), ...
            'MarkerEdgeColor', max(0, squeeze(colors(condIndex,:)) - faceColor), ...
            'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
    end

    
    if (showLegend)
       % Add legend
       hL = legend(theAxes, variedParamLegends, 'Location', 'NorthEast');
    else
       hL = []; 
    end
    % Add text
    if (~isempty(inGraphText))
        if (isempty(inGraphTextPos))
            if (plotRatiosOfOtherConditionsToFirst)
                inGraphTextPos(1) = 550;
                inGraphTextPos(2) = 450;
            else
                inGraphTextPos(1) = 550;
                inGraphTextPos(2) = 450;
            end
        end
        t = text(theAxes, inGraphTextPos(1), inGraphTextPos(2), inGraphText);
    else
        t = '';
    end
    hold(theAxes, 'off');
    
    if (plotRatiosOfOtherConditionsToFirst)
        
        hold(theRatioAxes, 'on');
        refTrials = theData{1}.trialsNum;
        refSensitivity = theData{1}.contrastSensitivity;
        for condIndex = 2:numel(theData)
            condTrials = theData{condIndex}.trialsNum;
            condSensitivity = theData{condIndex}.contrastSensitivity;
            ratios = zeros(1,numel(refTrials ));
            for k = 1:numel(condTrials)
                idx = find(condTrials(k) == refTrials);
                if (~isempty(idx))
                    ratios(k) = condSensitivity(idx)./refSensitivity(idx);
                end
            end
            plot(theRatioAxes, refTrials, ratios, ...
                'ko-', 'Color', squeeze(colors(condIndex,:)), ...
                'MarkerEdgeColor', squeeze(colors(condIndex,:)), ...
                'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
                'MarkerSize',12,'LineWidth',2);
        end
        hold(theRatioAxes, 'off');
    end
    
    
    % Format figure
    formatFigureForPaper(hFig, ...
        'figureType', 'CSF', ...
        'plotRatiosOfOtherConditionsToFirst', plotRatiosOfOtherConditionsToFirst, ...
        'theAxes', theAxes, ...
        'theRatioAxes', theRatioAxes, ...
        'theRatioLims', theRatioLims, ...
        'theRatioTicks', theRatioTicks, ...
        'theLegend', hL, ...
        'theText', t, ...
        'theTextFontSize', inGraphTextFontSize);
    
    
    set(theAxes, 'XTick', [300 1000 3000 10000 30000 100000],  ...
        'XLim', theTrialsLims, 'YLim', theCSLims);
    set(theRatioAxes, 'XTick', [300 1000 3000 10000 30000 100000],  ...
        'XLim', theTrialsLims);
    
    xlabel(theRatioAxes, 'number of trials', 'FontWeight', 'bold', 'FontSize', 20);
    
end
