function figGenerateInferenceEngineFactors()

    close all;
    % cd to wherver this script resides
    [localDir,~] = fileparts(which(mfilename()));
    cd(localDir)
    
    inferenceEngines = {...
        'mlpt' ...
        'svmV1FilterBank' ...
        'svm' ...
    };

    inferenceEngineLabels = {'MLPT', 'SVM-V1 filter bank', 'SVM (60 xtPCAs)'};


    % Trials computed
    nTrainingSamples = 1024;
    
    % Perhaps use a subset of the trials.
    performanceTrialsUsed = 1024;
    
    useRBFSVMKernel = false;
    
    % Spatial pooling kernel parameters
    spatialPoolingKernelParams.type = 'V1QuadraturePair';  % Choose between 'V1CosUnit' 'V1SinUnit' 'V1QuadraturePair';
    spatialPoolingKernelParams.activationFunction = 'energy';  % Choose between 'energy' and 'fullWaveRectifier'
    spatialPoolingKernelParams.adjustForConeDensity = false;
    spatialPoolingKernelParams.temporalPCAcoeffs = Inf;  % Inf, results in no PCA, just the raw time series
    spatialPoolingKernelParams.shrinkageFactor = 1.0;  % > 1, results in expansion, < 1 results in shrinking 
   
    
    mosaicParams = getParamsForMosaicWithLabel('ISETbioHexEccBasedLMS');
    
    % 'random'; 'frozen0';
    emPathType = 'frozen0'; %random'; %'random';     
    centeredEMPaths = false;
    
    opticsModel = 'WvfHumanMeanOTFmagMeanOTFphase';
    pupilDiamMm = 3.0;
    
    ramPercentageEmployed = 1.0;  % use all the RAM
    freezeNoise = true;
        
    cyclesPerDegreeExamined =  [5 10 20 40 50];
    luminancesExamined =  [34]; 

    responseStabilizationMilliseconds = 10;
    responseExtinctionMilliseconds = 50;
    integrationTimeMilliseconds =  5.0;
    lowContrast = 0.0001;
    highContrast = 0.3;
    nContrastsPerDirection =  18;
    
    % What to do ?
    visualizeSpatialScheme = ~true;
    findPerformance = ~true;
    visualizePerformance = true;
    visualizeTransformedSignals = ~true;
    
    plottedSFindices = [];
    
    for inferenceEngineIndex = 1:numel(inferenceEngines)
        performanceClassifier = inferenceEngines{inferenceEngineIndex};
        performanceSignal = 'isomerizations';
        [d, psychometricFunctions] = getData();
        for sfIndex = 1:numel(psychometricFunctions)
            psychoFuncionData{inferenceEngineIndex, sfIndex} = psychometricFunctions{sfIndex}{1};
        end

        sfTest{inferenceEngineIndex} = d.cyclesPerDegree;
        for luminanceIndex = 1:size(d.mlptThresholds,1)
            inferenceCSF(luminanceIndex,inferenceEngineIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
        end
        
        if (inferenceEngineIndex > 1)
            if (any(sfTest{inferenceEngineIndex}-sfTest{1})~=0)
                error('sfs do not match');
            end
        end
        ratioCSF{inferenceEngineIndex} = squeeze(inferenceCSF(1,inferenceEngineIndex,:)) ./ squeeze(inferenceCSF(1,1,:));
    end % inferenceEngineIndex
    
    
    
    inferenceColors = [0 0 0; 1 0 0; 0.4 0.7 1];
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', 2, ...
           'rowsNum', 2, ...
           'heightMargin',   0.09, ...
           'widthMargin',    0.07, ...
           'leftMargin',     0.07, ...
           'rightMargin',    0.005, ...
           'bottomMargin',   0.07, ...
           'topMargin',      0.02);
       
    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 960 880]);

    %% Psychometric functions
    xTicks = -4:0.5:0.5;
    yTicks = 0.0:0.1:0.9; 
    xTickLabels = xTicks;
    yTickLabels = yTicks;
    theContrastLims = [-4.05 -0.45];
    theFractionCorrectLims = [0.19 1.0];
    sfs = sfTest{1};
    
    sfIndex = 1;
    subplot('Position', subplotPosVectors(1,1).v);
    legends = {};
    hold on
    for inferenceEngineIndex = 1:numel(inferenceEngines)
        
        if (inferenceEngineIndex == 1)
            % reference inference engine
            markerType  = 's';
            markerSize = 13;
        else
            markerType = 'o';
            markerSize = 10;
        end
        legends = cat(2, legends, sprintf('%s', inferenceEngineLabels{inferenceEngineIndex}));
        d = psychoFuncionData{inferenceEngineIndex, sfIndex};
        plot(log10(d.x), d.y,  markerType, 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', squeeze(inferenceColors(inferenceEngineIndex,:)));
    end
    for inferenceEngineIndex = 1:numel(inferenceEngines)
        d = psychoFuncionData{inferenceEngineIndex, sfIndex};
        plot(log10(d.xFit),d.yFit,'r','LineWidth', 2.0, 'Color', squeeze(inferenceColors(inferenceEngineIndex,:)));
        plot(log10(d.thresholdContrast*[1 1]), [0 d.criterionFraction],'--', 'LineWidth', 1.5, 'Color', squeeze(inferenceColors(inferenceEngineIndex,:)));
    end % inferenceEngineIndex
    
    titleLabel = sprintf('%2.1fcpd', sfs(sfIndex));
    
    hold off;
    panelLabelPosAndColor{1} = '(A)';
    panelLabelPosAndColor{2} = [-4.75 0.99];
    panelLabelPosAndColor{3} = [0 0 0];
    finishPsychometricPlot(gca, 'log10(contrast)', 'fraction correct', theContrastLims, theFractionCorrectLims, xTicks, yTicks, xTickLabels, yTickLabels, legends, 'SouthWest', panelLabelPosAndColor, titleLabel);
   

    sfIndex = 4;
    subplot('Position', subplotPosVectors(1,2).v);
    
    hold on
    for inferenceEngineIndex = 1:numel(inferenceEngines)
        sfs = sfTest{inferenceEngineIndex};
        if (inferenceEngineIndex == 1)
            % reference inference engine
            markerType  = 's';
            markerSize = 13;
        else
            markerType = 'o';
            markerSize = 10;
        end
        d = psychoFuncionData{inferenceEngineIndex, sfIndex};
        plot(log10(d.x), d.y,  markerType, 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', squeeze(inferenceColors(inferenceEngineIndex,:)));
    end
    for inferenceEngineIndex = 1:numel(inferenceEngines)
        d = psychoFuncionData{inferenceEngineIndex, sfIndex};
        plot(log10(d.xFit),d.yFit,'r','LineWidth', 2.0, 'Color', squeeze(inferenceColors(inferenceEngineIndex,:)));
        plot(log10(d.thresholdContrast*[1 1]), [0 d.criterionFraction],'--', 'LineWidth', 1.5, 'Color', squeeze(inferenceColors(inferenceEngineIndex,:)));
    end % inferenceEngineIndex
    
    titleLabel = sprintf('%2.1fcpd', sfs(sfIndex));
    hold off;
    panelLabelPosAndColor{1} = '(B)';
    panelLabelPosAndColor{2} = [-4.75 0.99];
    panelLabelPosAndColor{3} = [0 0 0];
    finishPsychometricPlot(gca, 'log10(contrast)', '',theContrastLims, theFractionCorrectLims, xTicks, yTicks, xTickLabels, yTickLabels, legends, 'SouthWest', panelLabelPosAndColor, titleLabel);
    
    
    sfRange = [2 71];
    csfRange = [2.1 5000];
    xTickLabels = [1 3 10 30 100];
    xTicks      = xTickLabels;
    yTickLabels = [1 3 10 30 100 300 1000];
    yTicks      = yTickLabels;
    xLim = sfRange;
    yLim = csfRange;
    
    %% CSFs
    subplot('Position', subplotPosVectors(2,1).v);
        
    legends = {};
    hold on;
    for inferenceEngineIndex = 1:numel(inferenceEngines)
        
        if (inferenceEngineIndex == 1)
            % reference inference engine
            markerType  = 's-';
            markerSize = 13;
        else
            markerType = 'o-';
            markerSize = 10;
        end
        title = ' ';
        legends = cat(2, legends, sprintf('%s', inferenceEngineLabels{inferenceEngineIndex}));
        plot(sfTest{1}, squeeze(inferenceCSF(1,inferenceEngineIndex,:)), markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', [0 0 0], ...
            'MarkerFaceColor', squeeze(inferenceColors(inferenceEngineIndex,:)), ...
            'LineWidth', 1.5, ...
            'Color', squeeze(inferenceColors(inferenceEngineIndex,:)));
    end % inferenceEngineIndex
    hold off;
    panelLabelPosAndColor{1} = '(C)';
    panelLabelPosAndColor{2} = [1.0 3600];
    panelLabelPosAndColor{3} = [0 0 0];
                
    finishPlot(gca, 'spatial frequency (cpd)', 'contrast sensitivity', ...
        'log', 'log', xLim, yLim, xTicks, yTicks,  ...
        xTickLabels, yTickLabels, panelLabelPosAndColor, legends, 'SouthWest', true, false, [1 1 1]);
    
    %% CSF ratios
    subplot('Position', subplotPosVectors(2,2).v);
    hold on
    legends = {};
    for inferenceEngineIndex = 2:numel(inferenceEngines)
        
        if (inferenceEngineIndex == 1)
            % reference inference engine
            markerType  = 's-';
            markerSize = 13;
        else
            markerType = 'o-';
            markerSize = 10;
        end
        
        plot(sfTest{1}, ratioCSF{inferenceEngineIndex}, markerType, ...
            'MarkerSize', markerSize, 'MarkerFaceColor', squeeze(inferenceColors(inferenceEngineIndex,:)), 'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5, ...
            'Color', squeeze(inferenceColors(inferenceEngineIndex,:)));
        legends = cat(2, legends, sprintf('%s: %s', inferenceEngineLabels{inferenceEngineIndex}, inferenceEngineLabels{1}));
    end % inferenenceEngineIndex
    
    yTickLabels =  [0.03 0.05 0.1 0.2];
    yTicks      = yTickLabels;
    yLim        = [0.021 0.3];
    
    panelLabelPosAndColor{1} = '(D)';
    panelLabelPosAndColor{2} = [1.0 0.27];
    panelLabelPosAndColor{3} = [0 0 0];
    finishPlot(gca, 'spatial frequency (cpd)', 'contrast sensitivity ratio', ...
        'log', 'log', xLim, yLim, xTicks, yTicks, ...
        xTickLabels, yTickLabels, panelLabelPosAndColor, legends, 'SouthWest', true, false, [1 1 1]);
    
    drawnow
    cd(localDir);
    NicePlot.exportFigToPDF(sprintf('InferenceEngine.pdf'), hFig, 300);
    
    
    function [d,psychometricFunctions]  = getData()
        [d, ~, ~, ~, ~, psychometricFunctions] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', opticsModel, ...
            'pupilDiamMm', pupilDiamMm, ...
            'cyclesPerDegree', cyclesPerDegreeExamined, ...
            'luminances', luminancesExamined, ...
            'nTrainingSamples', nTrainingSamples, ...
            'nContrastsPerDirection', nContrastsPerDirection, ...
            'lowContrast', lowContrast, ...
            'highContrast', highContrast, ...
            'ramPercentageEmployed', ramPercentageEmployed, ...
            'emPathType', emPathType, ...
            'centeredEMPaths', centeredEMPaths, ...
            'responseStabilizationMilliseconds', responseStabilizationMilliseconds, ...
            'responseExtinctionMilliseconds', responseExtinctionMilliseconds, ...
            'freezeNoise', freezeNoise, ...
            'integrationTime', integrationTimeMilliseconds/1000, ...
            'coneSpacingMicrons', mosaicParams.coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', mosaicParams.innerSegmentDiameter, ...
            'conePacking', mosaicParams.conePacking, ...
            'LMSRatio', mosaicParams.LMSRatio, ...
            'mosaicRotationDegs', mosaicParams.mosaicRotationDegs, ...
            'computeMosaic', false, ...
            'visualizeMosaic', false, ...
            'computeResponses', false, ...
            'visualizeResponses', false, ...
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'findPerformance', findPerformance, ...
            'visualizePerformance', visualizePerformance, ...
            'visualizeKernelTransformedSignals', visualizeTransformedSignals, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier, ...
            'useRBFSVMKernel', useRBFSVMKernel, ...
            'performanceTrialsUsed', nTrainingSamples, ...
            'spatialPoolingKernelParams', spatialPoolingKernelParams ...
            );
    end

end


function finishPsychometricPlot(gca, theXlabel, theYlabel, XLim, YLim, xTicks, yTicks, xTickLabels, yTickLabels, theLegends, legendLocation, panelLabelPosAndColor, titleLabel)

    set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'normal', 'FontAngle', 'normal', 'TickDir', 'both');
    set(gca, 'XTick', xTicks, 'XTickLabel', xTickLabels, 'YTick', yTicks, 'yTickLabel', yTickLabels);
    set(gca, 'XLim',XLim,'YLim',YLim);
    
    grid on
    box off
    if (~isempty(theXlabel))
        t = xlabel(theXlabel);
        set(t, 'FontWeight', 'bold', 'FontAngle', 'normal');
    end
    if (~isempty(theYlabel))
        t = ylabel(theYlabel);
        set(t, 'FontWeight', 'bold', 'FontAngle', 'normal');
    end
    
    axis 'square'
    if (~isempty(theLegends))
        hL = legend(theLegends, 'Location', legendLocation);
        set(hL, 'FontSize', 13, 'FontWeight', 'Bold', 'FontAngle', 'italic', 'FontName', 'Menlo', 'Box', 'off');
    end
    title(titleLabel);
    hTitle=get(gca,'title');
    set(hTitle, 'FontSize', 18);
    
    % Label plot
    panelLabel = panelLabelPosAndColor{1};
    panelLabelPos = panelLabelPosAndColor{2};
    text(panelLabelPos(1), panelLabelPos(2), panelLabel, 'Color', panelLabelPosAndColor{3}, 'FontSize', 30, 'FontWeight', 'bold');
    
end

function finishPlot(gca, theXlabel, theYlabel, theXscale, theYscale, theXLim, theYLim, xTicks, yTicks, xTickLabels, yTickLabels, panelLabelPosAndColor, theLegends, legendLocation, showGrid, showBox, backgroundColor)
    set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'normal', 'FontAngle', 'normal', 'TickDir', 'both');
    set(gca, 'XTick', xTicks, 'XTickLabel', xTickLabels, 'YTick', yTicks, 'yTickLabel', yTickLabels);
    if (~isempty(theXlabel))
        t = xlabel(theXlabel);
        set(t, 'FontWeight', 'bold', 'FontAngle', 'normal');
    end
    if (~isempty(theYlabel))
        t = ylabel(theYlabel);
        set(t, 'FontWeight', 'bold', 'FontAngle', 'normal');
    end
    if (showGrid); grid on; else; grid off; end
    if (showBox); box on; else; box off; end
    set(gca, 'XScale',theXscale,'YScale',theYscale, 'LineWidth', 1.0);
    set(gca, 'XLim',theXLim,'YLim',theYLim);
    set(gca, 'Color', backgroundColor);
    axis square
    if (~isempty(theLegends))
        hL = legend(theLegends, 'Location', legendLocation);
        set(hL, 'FontSize', 13, 'FontWeight', 'Bold', 'FontAngle', 'italic', 'FontName', 'Menlo', 'Box', 'off');
    end
    hTitle=get(gca,'title');
    set(hTitle, 'FontSize', 18);
    
    % Label plot
    panelLabel = panelLabelPosAndColor{1};
    panelLabelPos = panelLabelPosAndColor{2};
    text(panelLabelPos(1), panelLabelPos(2), panelLabel, 'Color', panelLabelPosAndColor{3}, 'FontSize', 30, 'FontWeight', 'bold');
end

