function figGenerateOpticalFactors()
   
    opticsModels = {'Geisler', 'WvfHuman', 'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanSubject1'};
    
    % 'singleExposure'; 'timeSeries5msec'
    temporalAnalysis = 'timeSeries5msec';
    
    % 'random'; 'frozen0';
    emPathType = 'frozen0'; %random'; %'random';     
    centeredEMPaths = false;
    
    % 'isomerizations', 'photocurrents'
    performanceSignal ='isomerizations';
    
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
        [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, referenceMosaicLegend] = paramsForComparativeMosaicAnalysis(mosaicTest);
        [d, the_rParams{opticsModelIndex}, ~, theOIs{opticsModelIndex}] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
        sfRef = d.cyclesPerDegree;
        opticsCSF(opticsModelIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1)); 
    end % opticsModelIndex
    
        
    opticsColors = [0 0 0; 1 0 0; 0 0 1; 0 0.6 0; 0.5 0 0.7; 1.0 0.5 0.3; 0.3 0.5 1.0; 0.5 0.7 0.5];
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1360 450]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', 3, ...
           'rowsNum', 1, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.06, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.12, ...
           'topMargin',      0.02);
       
    
    subplot('Position', subplotPosVectors(1,2).v);
    addBanksEtAlReferenceLines(the_rParams{1}, plotLuminanceLines);
    hold on;
    refCSF = opticsCSF(1,:);
    legends = {'Banks et al'};
    for opticsModelIndex = 1:size(opticsCSF,1)
        plot(sfRef, opticsCSF(opticsModelIndex,:), 'o', 'MarkerSize', 10, 'MarkerFaceColor', squeeze(opticsColors(opticsModelIndex,:)));
        legends = cat(2, legends, sprintf('%s model',opticsModels{opticsModelIndex}));
    end
    hL = legend(legends, 'Location', 'SouthWest');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16);
    xlim([1 100]); ylim([1 3000]);
    grid on; box off;
    
    subplot('Position', subplotPosVectors(1,3).v);
    hold on
    legends = {};
    for opticsModelIndex = 2:size(opticsCSF,1)
        plot(sfRef, log10(opticsCSF(opticsModelIndex,:) ./ refCSF), 'o-', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', squeeze(opticsColors(opticsModelIndex,:)), 'Color', squeeze(opticsColors(opticsModelIndex,:)), 'MarkerSize', 12);
        legends = cat(2, legends, sprintf('%s : %s',opticsModels{opticsModelIndex}, opticsModels{1}));
    end
    hL = legend(legends, 'Location', 'NorthWest');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlim([1 100]); ylim(log10([0.25 4]));
    yTicks = log10([0.25 0.5 1 2 4]);
    set(gca,'XScale','log','YScale','linear', 'YTick', yTicks, 'YTickLabel', yTicks, 'FontSize', 16);
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('log contrast sensitivity ratio', 'FontSize' ,16, 'FontWeight', 'bold');
    grid on; box off;
    
    
    subplot('Position', subplotPosVectors(1,1).v);
    plotOTFs(theOIs, opticsModels, opticsColors);
    drawnow;
    
    function [d, the_rParams, noiseFreeResponse, theOI] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs)
        [d, the_rParams, noiseFreeResponse, theOI] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', opticsModel, ...
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
            'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter), ...
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

function plotOTFs(theOIs, opticsModels, opticsColors)
    wavelengthsList = [450 550 650];
    lineStyles = {'--', '-', '-.'};
    hold on
    legends = {};
    for opticsModelIndex = 1:numel(opticsModels)
        legends = cat(2, legends, sprintf('%s',opticsModels{opticsModelIndex}));
        for wIndex = 1:numel(wavelengthsList)
            [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOIs{opticsModelIndex}, wavelengthsList(wIndex));
            centerPosition = floor(size(otf,1)/2) + 1;
            otfSlice = squeeze(otf(centerPosition, centerPosition:end));
            % make otf_fx(1) = 0.1, so that the zero data point is plotted
            otf_fx(centerPosition) = 1;
            lineColor = squeeze(opticsColors(opticsModelIndex,:));
            plot(otf_fx(centerPosition:end), otfSlice, '-', 'LineStyle', lineStyles{wIndex}, 'LineWidth', 1.5, 'Color',lineColor);
        end
    end
    hold off
    xlim([1 100]); ylim([0 1]);
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('OTF mag', 'FontSize' ,16, 'FontWeight', 'bold');
    grid on; box off;
    set(gca,'XScale','log','YScale','linear',  'FontSize', 16);
end
