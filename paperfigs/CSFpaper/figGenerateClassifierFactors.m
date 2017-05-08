function figGenerateClassifierFactors()
    
    % 'singleExposure'; 'timeSeries5msec'
    temporalAnalysis = 'timeSeries5msec';
    
    % 'random'; 'frozen0';
    emPathType = 'frozen0'; %random'; %'random';     
    centeredEMPaths = false;
    
    % 'isomerizations', 'photocurrents'
    performanceSignal = 'isomerizations';
    
    % Use a subset of the trials. Specify [] to use all available trials
    nTrainingSamples = 1024;
    
    ramPercentageEmployed = 1.0;  % use all the RAM
    freezeNoise = ~true;
    
    cyclesPerDegreeExamined =  [5 10 20 40];
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
    mosaicRef = 'fullIsetbioWithScones';
    
    % Retrieve reference mosaic ('originalBanks') CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs] = paramsForComparativeMosaicAnalysis(mosaicRef);
    
    % Spatial pooling kernel parameters
    spatialPoolingKernelParams.type = 'V1QuadraturePair';  % Choose between 'V1CosUnit' 'V1SinUnit' 'V1QuadraturePair';
    spatialPoolingKernelParams.activationFunction = 'energy';  % Choose between 'energy' and 'fullWaveRectifier'
    spatialPoolingKernelParams.adjustForConeDensity = false;
    spatialPoolingKernelParams.temporalPCAcoeffs = Inf;  % Inf, results in no PCA, just the raw time series
    spatialPoolingKernelParams.shrinkageFactor = 1.0;  % > 1, results in expansion, < 1 results in shrinking
    useRBFSVMKernel = false;
    
    
    [d, rParams] = getData('mlpt');
    sfRef = d.cyclesPerDegree;
    correctionForGreaterThan100msecTime = 100/120
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        referenceCSF(luminanceIndex,:) = correctionForGreaterThan100msecTime * 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    
    % Spatial pooling kernel parameters
    spatialPoolingKernelParams.type = 'V1QuadraturePair';  % Choose between 'V1CosUnit' 'V1SinUnit' 'V1QuadraturePair';
    spatialPoolingKernelParams.activationFunction = 'energy';  % Choose between 'energy' and 'fullWaveRectifier'
    spatialPoolingKernelParams.adjustForConeDensity = false;
    spatialPoolingKernelParams.temporalPCAcoeffs = Inf;  % Inf, results in no PCA, just the raw time series
    spatialPoolingKernelParams.shrinkageFactor = 1.0;  % > 1, results in expansion, < 1 results in shrinking
    useRBFSVMKernel = false;
    
    [d, rParams] = getData('svmV1FilterBank');
    sfTest = d.cyclesPerDegree;
    if (any(sfTest-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF1(luminanceIndex,:) = correctionForGreaterThan100msecTime * 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    spatialPoolingKernelParams.adjustForConeDensity = true;
    [d, rParams] = getData('svmV1FilterBank');
    sfTest = d.cyclesPerDegree;
    if (any(sfTest-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF2(luminanceIndex,:) = correctionForGreaterThan100msecTime * 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    [d, rParams] = getData('svm');
    sfTest = d.cyclesPerDegree;
    if (any(sfTest-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF3(luminanceIndex,:) = correctionForGreaterThan100msecTime * 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    % Compute CSF ratios
    ratioCSF1 = testCSF1 ./ referenceCSF;
    ratioCSF2 = testCSF2 ./ referenceCSF;
    ratioCSF3 = testCSF3 ./ referenceCSF;
    
    hFig = figure(2); clf;
    set(hFig ,'Position',[100 100 450 280], 'Color', [1 1 1]);
    subplot('Position', [0.13 0.08 0.84 0.90]);
    hold on
    plot(sfRef, ratioCSF1, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.5 0.5], 'MarkerSize', 12, 'LineWidth', 1.0);
    plot(sfRef, ratioCSF2, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 0.9], 'MarkerSize', 12, 'LineWidth', 1.0);
    plot(sfRef, ratioCSF3, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.7 0.2], 'MarkerSize', 12, 'LineWidth', 1.0);
    hL = legend({'SVM (quad filter, no cone density adjust)', 'SVM (quad filter, cone density adjust)', 'SVM'}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity ratio', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16);
    xlim([1 100]); ylim([0.01 1]);
    grid on; box off;
    drawnow;

    function [d, the_rParams] = getData(performanceClassifier)
        [d, the_rParams] = c_BanksEtAlPhotocurrentAndEyeMovements(...
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

