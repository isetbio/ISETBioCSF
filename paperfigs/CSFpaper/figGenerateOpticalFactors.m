function figGenerateOpticalFactors()
   
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
    freezeNoise = ~true;
    
    
    cyclesPerDegreeExamined =  [2.5 5 10 20 40];
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
            responseExtinctionMilliseconds = 150;
            integrationTimeMilliseconds =  5.0;
            lowContrast = 0.0005;
            highContrast = 0.8;
            nContrastsPerDirection =  18;
    end
   
    % What to do ?
    visualizeSpatialScheme = ~true;
    findPerformance = ~true;
    visualizePerformance = true;
    visualizeTransformedSignals = ~true;
    
    % 'originalBanks'; 'defaultIsetbio';  'fullIsetbioNoScones'; 'fullIsetbioWithScones'
    mosaicRef = 'originalBanks';
    mosaicTest = 'fullIsetbioWithScones';
    
    luminanceIndex = 1;
    

    opticsModels = {'Geisler', 'WvfHuman'};
    
    for opticsModelIndex = 1:numel(opticsModels)
        opticsModel = opticsModels{opticsModelIndex};
        
        % Retrieve reference mosaic ('originalBanks') CSF data
        [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, referenceMosaicLegend] = paramsForComparativeMosaicAnalysis(mosaicRef);
        [d, rParamsRef{opticsModelIndex}] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
        sfRef = d.cyclesPerDegree;
        banksMosaicCSF(opticsModelIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));

    
        % Retrieve 'fullIsetbioWithScones' CSF data
        [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, mosaicTestLegend] = paramsForComparativeMosaicAnalysis(mosaicTest);
        [d, rParamsTest{opticsModelIndex}] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
        test1Color = [0.8 0.3 0.9];
    
        sfTest = d.cyclesPerDegree;
        if (any(sfTest-sfRef)~=0)
            error('sfs do not match');
        end
        testCSF(opticsModelIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end % opticsModelIndex
    
    % Retrieve reference mosaic ('originalBanks') CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, referenceMosaicLegend] = paramsForComparativeMosaicAnalysis('defaultIsetbio');
    [d, rParamsTest2{1}] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
    sfRef = d.cyclesPerDegree;
    defaultISETbioCSF(1,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));

        
        
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 992 1178]);
    subplot(2,2,1);
    addBanksEtAlReferenceLines(rParamsTest{opticsModelIndex}, plotLuminanceLines);
    hold on;
    refCSF = banksMosaicCSF(1,:);
    plot(sfRef, banksMosaicCSF(1,:), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    plot(sfRef, banksMosaicCSF(2,:), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 12);
    hL = legend({'Banks study', 'Banks mosaic: Geisler optics', 'Banks mosaic: wvf optics'}, 'Location', 'NorthOutside');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16);
    xlim([1 100]); ylim([1 3000]);
    grid on; box off;
    
    subplot(2,2,3);
    plot(sfRef, banksMosaicCSF(2,:) ./ refCSF, 'ko-', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    xlim([1 100]); ylim([1 2.8]);
    set(gca,'XScale','log','YScale','linear', 'FontSize', 16);
    title('wvf optics: Geisler optics');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity ratio  (Banks mosaic)', 'FontSize' ,16, 'FontWeight', 'bold');
    grid on; box off;
    
    subplot(2,2,2)
    addBanksEtAlReferenceLines(rParamsTest{opticsModelIndex}, plotLuminanceLines);
    hold on;
    refCSF = testCSF(1,:);
    plot(sfRef, testCSF(1,:), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    plot(sfRef, testCSF(2,:), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 12);
    hL = legend({'Banks study', 'ISETbio mosaic: Geisler optics', 'ISETbio mosaic: wvf optics'}, 'Location', 'NorthOutside');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16);
    xlim([1 100]); ylim([1 3000]);
    grid on; box off;
    
    subplot(2,2,4);
    plot(sfRef, testCSF(2,:) ./ refCSF, 'ko-', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    xlim([1 100]); ylim([1 2.8]);
    set(gca,'XScale','log','YScale','linear', 'FontSize', 16);
    title('wvf optics: Geisler optics');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity ratio (ISETbio mosaic)', 'FontSize' ,16, 'FontWeight', 'bold');
    grid on; box off;
    
    
    
    drawnow;
    
    function [d, the_rParams] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs)
        [d, the_rParams] = c_BanksEtAlPhotocurrentAndEyeMovements(...
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
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'findPerformance', findPerformance, ...
            'visualizePerformance', visualizePerformance, ...
            'visualizeTransformedSignals', visualizeTransformedSignals, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier, ...
            'useRBFSVMKernel', useRBFSVMKernel, ...
            'performanceTrialsUsed', nTrainingSamples, ...
            'spatialPoolingKernelParams', spatialPoolingKernelParams ...
            );
    end

end

