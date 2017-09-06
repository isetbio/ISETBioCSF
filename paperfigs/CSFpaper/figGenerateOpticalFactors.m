function figGenerateOpticalFactors()
   
    close all;
    
    pupilDiamMm = 2;
    %pupilDiamMm = 3;
    
    if (pupilDiamMm == 2)
        opticsModels = {'Geisler', 'WvfHuman', 'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanSubject1', 'WvfHumanSubject2', 'WvfHumanSubject3', 'WvfHumanSubject4', 'WvfHumanSubject5'};
        opticsModelsLabels = {'Geisler', 'wvfMeanCoeff', 'wvfMeanOTF  ', 'wvfSubject1 ', 'wvfSubject2 ', 'wvfSubject3 ', 'wvfSubject4 ', 'wvfSubject5 '};
    else
        opticsModels = {'Geisler',  'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanSubject1', 'WvfHumanSubject2'};%, 'WvfHumanSubject3'};
        opticsModelsLabels = {'Geisler',  'wvfMeanOTF  ', 'wvfSubject1 ', 'wvfSubject2 '}; % , 'wvfSubject3 '};
    end
    
    blur = true;
    apertureBlur = false;
    
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
    
    cyclesPerDegreeExamined =  [2.5 5 10 20 40 50];
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
    visualizeKernelTransformedSignals = ~true;
    
    % 'originalBanks'; 'defaultIsetbio';  'fullIsetbioNoScones'; 'fullIsetbioWithScones'
    mosaicTest = 'originalBanks';
    luminanceIndex = 1;

    for opticsModelIndex = 1:numel(opticsModels)
        opticsModel = opticsModels{opticsModelIndex};
        % Retrieve reference mosaic ('originalBanks') CSF data
        params = getParamsForMosaicWithLabel(mosaicTest);
        coneSpacingMicrons = params.coneSpacingMicrons;
        innerSegmentDiameter = params.innerSegmentDiameter;
        conePacking = params.conePacking;
        LMSRatio = params.LMSRatio;
        mosaicRotationDegs = params.mosaicRotationDegs;
    
        if (innerSegmentDiameter==3)
            fprintf(2,'Setting inner segment diam to 2.6587 - optics computations performed before changing the code\n');
            innerSegmentDiameter = 2.6587
        end
        
        [d, the_rParams{opticsModelIndex}, ~, theOIs{opticsModelIndex}] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
        sfRef = d.cyclesPerDegree;
        opticsCSF(opticsModelIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1)); 
    end % opticsModelIndex
    
        
    opticsColors = [0 0 0; 1 0 0; 0.4 0.7 1];
    opticsColors = cat(1, opticsColors, 0.85*[0.95 0.95 0.4] + 0.15*(1-[0.2 0.95 0.9]));
    opticsColors = cat(1, opticsColors, 0.7*[0.75 0.75 0.4] + 0.3*(1-[0.2 0.75 0.75]));
    opticsColors = cat(1, opticsColors, 0.55*[0.55 0.55 0.4] + 0.45*(1-[0.2 0.55 0.55]));
    opticsColors = cat(1, opticsColors, 0.4*[0.35 0.35 0.4] + 0.6*(1-[0.2 0.35 0.36]));
    opticsColors = cat(1, opticsColors, 0.25*[0.1 0.1 0.4] + 0.75*(1-[0.2 0.1 0.1]));
    
    sfRange = [2 61];
    wavelengthsList = [550 450 650 700];
    
    hFig = figure(1); clf;
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
    
    %% OTF slices
    xTicks = [1 3 10 30 100];
    yTicks = 0:0.2:1;
    xTickLabels = xTicks;
    yTickLabels = yTicks;
    xLim = sfRange;
    yLim = [-0.02 1.0];
    backgroundColor = [1 1 1];
    
    subplot('Position', subplotPosVectors(1,1).v);
    plotOTFs(theOIs, opticsModelsLabels, opticsColors, sfRange, wavelengthsList(1));
    finishPlot(gca, '', 'modulation transfer function', 'log', 'linear', ...
        xLim, yLim, xTicks, yTicks, xTickLabels, yTickLabels, {}, '', true, false, backgroundColor);
    
    subplot('Position', subplotPosVectors(1,2).v);
    plotOTFs(theOIs, opticsModelsLabels, opticsColors, sfRange, wavelengthsList(2));
    finishPlot(gca, '', '', 'log', 'linear', ...
        xLim, yLim, xTicks, yTicks, xTickLabels, yTickLabels, {}, '', true, false, backgroundColor);
    
    subplot('Position', subplotPosVectors(2,1).v);
    plotOTFs(theOIs, opticsModelsLabels, opticsColors, sfRange, wavelengthsList(3));
    finishPlot(gca, 'spatial frequency (cpd)', 'modulation transfer function', 'log', 'linear', ...
        xLim, yLim, xTicks, yTicks, xTickLabels, yTickLabels, {}, '', true, false, backgroundColor);
    
    subplot('Position', subplotPosVectors(2,2).v);
    plotOTFs(theOIs, opticsModelsLabels, opticsColors, sfRange, wavelengthsList(4));
    finishPlot(gca, 'spatial frequency (cpd)', '', 'log', 'linear', ...
        xLim, yLim, xTicks, yTicks, xTickLabels, yTickLabels, {}, '', true, false, backgroundColor);
    
    %% CSFs
    subplot('Position', subplotPosVectors(3,1).v);
    addBanksEtAlReferenceLines(the_rParams{1}, plotLuminanceLines);
    xTickLabels = [1 3 10 30 100];
    xTicks      = xTickLabels;
    yTickLabels = [1 3 10 30 100 300 1000 3000];
    yTicks      = yTickLabels;
    xLim = sfRange;
    yLim = [2.1 3000];
    
    
    hold on;
    refCSF = opticsCSF(1,:);
    legends = {'Banks et al. (87)'};
    for opticsModelIndex = 1:size(opticsCSF,1)
        if (opticsModelIndex <= 3)
            markerType  = 's';
            markerSize = 13;
        else
            markerType = 'o';
            markerSize = 8;
        end

        plot(sfRef, opticsCSF(opticsModelIndex,:), markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', squeeze(opticsColors(opticsModelIndex,:)));
        legends = cat(2, legends, opticsModelsLabels{opticsModelIndex});
    end
    hold off;
    finishPlot(gca, 'spatial frequency (cpd)', 'contrast sensitivity', 'log', 'log', ...
        xLim, yLim, xTicks, yTicks, xTickLabels, yTickLabels, ...
        legends, 'SouthWest', true, false, backgroundColor);
    
    %% CSF ratios
    subplot('Position', subplotPosVectors(3,2).v);
    hold on
    legends = {};
    for opticsModelIndex = 2:size(opticsCSF,1)
        if (opticsModelIndex <= 3)
            markerType  = 's-';
            markerSize = 13;
        else
            markerType = 'o-';
            markerSize = 8;
        end
   
        plot(sfRef, opticsCSF(opticsModelIndex,:) ./ refCSF, markerType, 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', squeeze(opticsColors(opticsModelIndex,:)));
        legends = cat(2, legends, sprintf('%s: %s',opticsModelsLabels {opticsModelIndex}, opticsModels{1}));
    end
    hold off
    yTickLabels =  [0.125 0.25 0.5 1 2];
    yTicks      = yTickLabels;
    yLim        = [0.21 4.0];
    
    finishPlot(gca, 'spatial frequency (cpd)', 'contrast sensitivity ratio', 'log', 'log', ...
        xLim, yLim, xTicks, yTicks, xTickLabels, yTickLabels, ...
        legends, 'SouthWest', true, false, backgroundColor);
    drawnow;
    NicePlot.exportFigToPDF(sprintf('OpticsImpactPupilMM%2.1f.pdf',pupilDiamMm), hFig, 300);
    
    
    function [d, the_rParams, noiseFreeResponse, theOI] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs)
        [d, the_rParams, noiseFreeResponse, theOI] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', opticsModel, ...
            'pupilDiamMm', pupilDiamMm, ...
            'blur', blur, ...
            'apertureBlur', apertureBlur, ...
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
            'visualizeOptics', true, ...
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'findPerformance', findPerformance, ...
            'visualizePerformance', visualizePerformance, ...
            'visualizeKernelTransformedSignals', visualizeKernelTransformedSignals, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier, ...
            'useRBFSVMKernel', useRBFSVMKernel, ...
            'performanceTrialsUsed', nTrainingSamples, ...
            'spatialPoolingKernelParams', spatialPoolingKernelParams ...
            );
    end

end

function plotOTFs(theOIs, opticsModelLabels, opticsColors, sfRange, targetWavelength)
    
    lineStyle = '-';
    hold on
    legends = {};
    for opticsModelIndex = 1:numel(opticsModelLabels)
        legends = cat(2, legends, opticsModelLabels{opticsModelIndex});
        [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOIs{opticsModelIndex}, targetWavelength);
        centerPosition = floor(size(otf,1)/2) + 1;
        otfSlice = squeeze(otf(centerPosition, centerPosition:end));
        % make otf_fx(1) = 0.1, so that the zero data point is plotted
        otf_fx(centerPosition) = 0.1;
        lineColor = squeeze(opticsColors(opticsModelIndex,:));
        plot(otf_fx(centerPosition:end), otfSlice, '-', 'LineStyle', lineStyle, 'LineWidth', 1.5, 'Color',lineColor);
    end
    hold off
    xlim(sfRange); ylim([-0.03 1]);
    title(sprintf('%d nm', targetWavelength));
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
    set(hTitle, 'FontSize', 14);
end
