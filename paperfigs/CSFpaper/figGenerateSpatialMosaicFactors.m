function figGenerateSpatialMosaicFactors()
    
    close all;

    mosaicModels       = {'originalBanksNoRotation',  'ISETbioHexRegNoScones',  'ISETbioHexRegLMS',     'ISETbioHexEccBasedLMS'};
    mosaicModelsFullLabels = {'Banks''87, 3um/3um',   'ISETbio 2um/1.6um',    'ISETbio 2um/1.6um',  'ISETbio ecc-based/1.6um'};
    mosaicModelsLabels = {'Banks''87',   'ISETbio 2um/1.6um (LM)',    'ISETbio 2um/1.6um (LMS)',  'ISETbio ecc-based/1.6um (LMS)'};
    mosaicConeTypes = {'LM', 'LM', 'LMS', 'LMS'};
    
    % 'WvfHuman', 'Geisler'
    opticsModel = 'WvfHumanMeanOTFmagMeanOTFphase';
    pupilDiamMm = 3.0;
    
    % 'singleExposure'; 'timeSeries5msec'
    temporalAnalysis = 'timeSeries5msec';
    
    % 'random'; 'frozen0';
    emPathType = 'frozen0'; %random'; %'random';     
    centeredEMPaths = false;
    
    % 'isomerizations', 'photocurrents'
    performanceSignal = 'isomerizations';
    
    % Use a subset of the trials. Specify [] to use all available trials
    nTrainingSamples = 1024;
    
    % 'mlpt', 'svm', 'svmV1FilterBank'
    performanceClassifier = 'mlpt'; %'svmV1FilterBank'; %'mlpt'% 'svmV1FilterBank';
    
    % Spatial pooling kernel parameters
    spatialPoolingKernelParams.type = 'V1QuadraturePair';  % Choose between 'V1CosUnit' 'V1SinUnit' 'V1QuadraturePair';
    spatialPoolingKernelParams.activationFunction = 'energy';  % Choose between 'energy' and 'fullWaveRectifier'
    spatialPoolingKernelParams.adjustForConeDensity = false;
    spatialPoolingKernelParams.temporalPCAcoeffs = Inf;  % Inf, results in no PCA, just the raw time series
    spatialPoolingKernelParams.shrinkageFactor = 1.0;  % > 1, results in expansion, < 1 results in shrinking
    useRBFSVMKernel = false;
    ramPercentageEmployed = 1.0;  % use all the RAM
    freezeNoise = true;
    
    
    cyclesPerDegreeExamined =  [5 10 20 40 50];
    luminancesExamined =  [34]; 
    plotLuminanceLines = [false true false];

    switch temporalAnalysis 
        case 'singleExposure'
            % For 100 ms simulation
            responseStabilizationMilliseconds = 0;
            responseExtinctionMilliseconds = 100;
            integrationTimeMilliseconds =  100;
            lowContrast = 0.001;
            highContrast = 0.3;
            nContrastsPerDirection =  12;    
            
        case 'timeSeries5msec'
            % For 7.0 milliseconds simulation
            responseStabilizationMilliseconds = 10;
            responseExtinctionMilliseconds = 50;
            integrationTimeMilliseconds =  5.0;
            lowContrast = 0.0001;
            highContrast = 0.3;
            nContrastsPerDirection =  18;
    end
   
    % What to do ?
    visualizeSpatialScheme = ~true;
    findPerformance = ~true;
    visualizePerformance = true;
    visualizeTransformedSignals = ~true;
    
    
    % Retrieve reference mosaic
    for mosaicIndex = 1:numel(mosaicModels)
        params = getParamsForMosaicWithLabel(mosaicModels{mosaicIndex});
        coneSpacingMicrons = params.coneSpacingMicrons;
        innerSegmentDiameter = params.innerSegmentDiameter;
        conePacking = params.conePacking;
        LMSRatio = params.LMSRatio;
        mosaicRotationDegs = params.mosaicRotationDegs;
        mosaicLegend{mosaicIndex} = sprintf('%s (s=%2.1fum, a=%2.1fum)', mosaicModelsLabels{mosaicIndex}, coneSpacingMicrons, innerSegmentDiameter);
        [d, rParams{mosaicIndex}] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
        sfTest{mosaicIndex} = d.cyclesPerDegree;
        for luminanceIndex = 1:size(d.mlptThresholds,1)
            mosaicCSF(luminanceIndex,mosaicIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
        end
        if (mosaicIndex > 1)
            if (any(sfTest{mosaicIndex}-sfTest{1})~=0)
                error('sfs do not match');
            end
        end
        ratioCSF{mosaicIndex} = squeeze(mosaicCSF(1,mosaicIndex,:)) ./ squeeze(mosaicCSF(1,1,:));
    end % mosaicIndex
    
    
    sfRange = [2 71];
    csfRange = [2.1 5000];
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 960 1400]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', 2, ...
           'rowsNum', 3, ...
           'heightMargin',   0.05, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.005, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.01);
       
    
    %% Extract the mosaic for some frequency
    sfIndex = 1;
    xTickLabels = [-200:20:200];
    xTicks      = xTickLabels*1e-6;
    yTicks      = xTicks;
    yTickLabels = xTickLabels;
    xLim = 45*[-1 1]*1e-6;
    yLim = xLim;
    backgroundColor = [0.1 0.1 0.1];
    for mosaicIndex = 1:numel(mosaicModels)
        switch (mosaicIndex)
            case 1
                subplot('Position', subplotPosVectors(1,1).v);
            case 2
                subplot('Position', subplotPosVectors(1,2).v);
            case 3
                subplot('Position', subplotPosVectors(2,1).v);
            case 4
                subplot('Position', subplotPosVectors(2,2).v);
        end % switch
        plotMosaic(rParams{mosaicIndex}, sfIndex, sprintf('%s', mosaicModelsFullLabels{mosaicIndex}), mosaicConeTypes{mosaicIndex}); 
        switch (mosaicIndex)
            case 1
                finishPlot(gca, '', 'space (microns)', ...
                    'linear', 'linear', xLim, yLim, ...
                    xTicks, yTicks, {}, yTickLabels, [], '', false, true, backgroundColor);
            case 2
                finishPlot(gca, '', '', ...
                    'linear', 'linear', xLim, yLim, ...
                    xTicks, yTicks, {}, {}, [], '', false,  true,  backgroundColor);
            case 3
                finishPlot(gca, 'space (microns)', 'space (microns)', ...
                    'linear', 'linear', xLim, yLim, ...
                    xTicks, yTicks, xTickLabels, yTickLabels, [], '', false, true,  backgroundColor);
            case 4
                finishPlot(gca, 'space (microns)', '', ...
                    'linear', 'linear', xLim, yLim, ...
                    xTicks, yTicks, xTickLabels, {}, [], '', false, true,  backgroundColor);
                subplot('Position', subplotPosVectors(2,2).v);
        end % switch
    end % mosaicIndex

    %% CSFs
    mosaicColors = [0 0 0; 1 0 0; 0.4 0.7 1; 0 0.8 0; 0 0.3 0.8];
    subplot('Position', subplotPosVectors(3,1).v);
    %addBanksEtAlReferenceLines(rParams{1}, plotLuminanceLines);
    hold on;
    legends = {}; %{'Banks et al. (87) - 2mm'};
    for mosaicModelIndex = 1:size(mosaicCSF,2)
        if (mosaicModelIndex == 1)
            % reference mosaic
            markerType  = 's-';
            markerSize = 13;
        else
            markerType = 'o-';
            markerSize = 8;
        end
        plot(sfTest{1}, squeeze(mosaicCSF(1,mosaicModelIndex,:)), markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', [0 0 0], ...
            'MarkerFaceColor', squeeze(mosaicColors(mosaicModelIndex,:)), ...
            'LineWidth', 1.5, ...
            'Color', squeeze(mosaicColors(mosaicModelIndex,:)));
        legends = cat(2, legends, mosaicModelsLabels{mosaicModelIndex});
    end % mosaicIndex
    
    hold off;
    xTickLabels = [1 3 10 30 100];
    xTicks      = xTickLabels;
    yTickLabels = [1 3 10 30 100 300 1000 3000];
    yTicks      = yTickLabels;
    xLim = sfRange;
    yLim = csfRange;
    
    finishPlot(gca, 'spatial frequency (cpd)', 'contrast sensitivity', ...
        'log', 'log', xLim, yLim, xTicks, yTicks, ...
        xTickLabels, yTickLabels, legends, 'SouthWest', true, false, [1 1 1]);

    
    %% CSF ratios
    subplot('Position', subplotPosVectors(3,2).v);
    hold on
    legends = {};
    for mosaicModelIndex = 2:size(mosaicCSF,2)
        if (mosaicModelIndex == 1)
            % reference mosaic
            markerType  = 's-';
            markerSize = 13;
        else
            markerType = 'o-';
            markerSize = 8;
        end
        plot(sfTest{1}, squeeze(mosaicCSF(1,mosaicModelIndex,:)) ./ squeeze(mosaicCSF(1,1,:)), markerType, ...
            'MarkerSize', markerSize, 'MarkerFaceColor', squeeze(mosaicColors(mosaicModelIndex,:)), 'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5, ...
            'Color', squeeze(mosaicColors(mosaicModelIndex,:)));
        legends = cat(2, legends, sprintf('%s: %s', mosaicModelsLabels{mosaicModelIndex}, mosaicModelsLabels{1}));
    end
    
    yTickLabels =  [0.125 0.25 0.5 1 2];
    yTicks      = yTickLabels;
    yLim        = [0.120 2.0];
    
    finishPlot(gca, 'spatial frequency (cpd)', 'contrast sensitivity ratio', ...
        'log', 'log', xLim, yLim, xTicks, yTicks, ...
        xTickLabels, yTickLabels, legends, 'SouthWest', true, false, [1 1 1]);
    
    drawnow
    NicePlot.exportFigToPDF(sprintf('ConeMosaicImpactPupilMM%2.1f.pdf',pupilDiamMm), hFig, 300);
     
%     conditionsToVisualize = [18];
%     cyclesPerDeg = cyclesPerDegreeExamined(sfIndex);
%     params = getParamsForMosaicWithLabel(mosaicModels{3});
%     coneSpacingMicrons = params.coneSpacingMicrons;
%     innerSegmentDiameter = params.innerSegmentDiameter;
%     conePacking = params.conePacking;
%     LMSRatio = params.LMSRatio;
%     mosaicRotationDegs = params.mosaicRotationDegs;
%     noiseFreeResponses = getNoiseFreeResponses(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, cyclesPerDeg, conditionsToVisualize);
%     noiseFreeIsomerizations = noiseFreeResponses{1,1}.stimData{conditionsToVisualize}.noiseFreeIsomerizations;
%     [~,idx] = max(noiseFreeIsomerizations(:));
%     [~, timeBinOfMaxResponse] = ind2sub(size(noiseFreeIsomerizations), idx);
%     noiseFreeIsomerizations = noiseFreeIsomerizations(:, timeBinOfMaxResponse);
   
%     figure(hFig);
%     theAxes = axes('Position',[0.06,0.305,0.2,0.2]);
%     plotMosaicActivation(theAxes, rParamsTest2, sfIndex, noiseFreeIsomerizations, 'isomerizations', cyclesPerDeg);


    function [d, the_rParams] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs)
        [d, the_rParams, noiseFreeResponse, theOI] = c_BanksEtAlPhotocurrentAndEyeMovements(...
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
            'coneSpacingMicrons', coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', innerSegmentDiameter, ...
            'conePacking', conePacking, ...
            'LMSRatio', LMSRatio, ...
            'mosaicRotationDegs', mosaicRotationDegs, ...
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

    function noiseFreeResponses = getNoiseFreeResponses(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, cyclesPerDeg, conditionsToVisualize)
        [~, ~, noiseFreeResponses, ~] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', opticsModel, ...
            'pupilDiamMm', pupilDiamMm, ...
            'cyclesPerDegree', cyclesPerDeg, ...
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
            'coneSpacingMicrons', coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', innerSegmentDiameter, ...
            'conePacking', conePacking, ...
            'LMSRatio', LMSRatio, ...
            'mosaicRotationDegs', mosaicRotationDegs, ...
            'computeMosaic', false, ...
            'visualizeMosaic', false, ...
            'computeResponses', false, ...
            'generatePlots', false, ...
            'visualizeResponses', true, ...
            'visualizedConditionIndices', conditionsToVisualize, ...
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

function plotMosaic(rParamsAllConds, sfIndex, titleString, coneTypes)
    lumIndex = 1;
    rParams = rParamsAllConds{sfIndex,lumIndex};
    fprintf('Loading a previously saved cone mosaic\n');
    coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    visualizedAperture = 'lightCollectingArea'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
    apertureShape = 'disks';                 % choose between 'hexagons' and 'disks'
    theMosaic.visualizeGrid('axesHandle', gca, ...
        'visualizedConeAperture', visualizedAperture, ...
        'apertureShape', apertureShape);    
    conesNum = numel(find(theMosaic.pattern(:) > 1));
    title(gca, sprintf('%s (%d %s cones)', titleString, conesNum, coneTypes));
end

function plotMosaicActivation(axesHandle, rParamsAllConds, sfIndex, activationMap, signalName, cpd)
    lumIndex = 1;
    rParams = rParamsAllConds{sfIndex,lumIndex};
    fprintf('Loading a previously saved cone mosaic\n');
    coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    if (strcmp(signalName, 'isomerizations'))
        % Make it rate in R*/sec
        activationMap = activationMap / theMosaic.integrationTime;
    end
    theMosaic.renderActivationMap(axesHandle, activationMap,...
        'mapType', 'modulated disks', ...
        'colorMap', gray(1024), ...
        'signalName', signalName, ...
        'xRange', [-120 120], ...
        'yRange', [-120 120], ...
        'signalRange', [0.09 1.61]*1000 ...
        );
    title(axesHandle, sprintf('%2.2f c/deg', cpd));
    
end

function addBanksEtAlReferenceLines(rParamsAllConds, plotLuminanceLine)
    sfIndex = 1;
    lumIndex = 1;
    rParams = rParamsAllConds{sfIndex,lumIndex};
    luminanceColors = [0 1 0; 0 0 1; 1 0 0; 0 0 0];
    luminanceColors = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
    % Add unshifted version for reference
    banksFactor = 1;
    [A,B,C,D,E] = LoadDigitizedBanksFigure2;
    if (~rParams.oiParams.blur && ~rParams.mosaicParams.apertureBlur)
        plot(A(:,1),A(:,2),'k:','LineWidth',0.5);
        plot(A(:,1),A(:,2)*banksFactor,'r-','LineWidth',2);
    elseif (~rParams.oiParams.blur)
        plot(B(:,1),B(:,2),'k:','LineWidth',0.5);
        plot(B(:,1),B(:,2)*banksFactor,'r','LineWidth',2);
    else
        
        if (plotLuminanceLine(1))
            % 3.4 cd/m2
            plot(E(:,1),E(:,2)*banksFactor,'-','LineWidth',1.0, 'Color', squeeze(luminanceColors(1,:)));
        end
        
        if (plotLuminanceLine(2))
            % 34 cd/m2
            plot(D(:,1),D(:,2)*banksFactor,'-','LineWidth',1.0, 'Color', squeeze(luminanceColors(2,:)));
        end
        
        if (plotLuminanceLine(3))
            % 340 cd/m2
            plot(C(:,1),C(:,2)*banksFactor,'-','LineWidth',1.0, 'Color', squeeze(luminanceColors(3,:)));
        end
    end
end


function finishPlot(gca, theXlabel, theYlabel, theXscale, theYscale, theXLim, theYLim, xTicks, yTicks, xTickLabels, yTickLabels, theLegends, legendLocation, showGrid, showBox, backgroundColor)
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
    if (showGrid); grid on; else grid off; end;
    if (showBox); box on; else box off; end;
    set(gca, 'XScale',theXscale,'YScale',theYscale, 'LineWidth', 1.0);
    set(gca, 'XLim',theXLim,'YLim',theYLim);
    set(gca, 'Color', backgroundColor);
    axis square
    if (~isempty(theLegends))
        hL = legend(theLegends, 'Location', legendLocation);
        set(hL, 'FontSize', 14, 'FontWeight', 'Bold', 'FontAngle', 'italic', 'FontName', 'Menlo', 'Box', 'off');
    end
    hTitle=get(gca,'title');
    set(hTitle, 'FontSize', 16);
end
