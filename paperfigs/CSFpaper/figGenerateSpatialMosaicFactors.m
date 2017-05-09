function figGenerateSpatialMosaicFactors()
    
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
    mosaicTest1 = 'fullIsetbioWithScones';
    %mosaicTest2 = 'fullIsetbioWithScones';
    
    % Retrieve reference mosaic ('originalBanks') CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs] = paramsForComparativeMosaicAnalysis(mosaicRef);
    [d, rParams] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
    sfRef = d.cyclesPerDegree;
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        referenceCSF(luminanceIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
   
    
   % Retrieve 'fullIsetbioWithScones' CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs] = paramsForComparativeMosaicAnalysis(mosaicTest1);
    [d, rParams] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
    sfTest = d.cyclesPerDegree;
    if (any(sfTest-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF1(luminanceIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    % Compute CSF ratios
    ratioCSF1 = testCSF1 ./ referenceCSF;
    
    hFig = figure(1); clf;
    set(hFig ,'Position',[100 100 450 1000], 'Color', [1 1 1]);
    subplot('Position', [0.13 0.37 0.84 0.60]);
    hold on
    addBanksEtAlReferenceLines(rParams, plotLuminanceLines);
    plot(sfRef, referenceCSF, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    hL = legend({'Banks et al', 'ISETbio (Banks)'}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16, 'XTickLabel', {});
    xlabel('', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontSize' ,16, 'FontWeight', 'bold');
    xlim([1 100]); ylim([1 10000]);
    grid on; box off;
    
    subplot('Position', [0.13 0.06 0.84 0.25]);
    hold on
    plot(sfRef, ratioCSF1, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.4 0.4 0.9], 'MarkerSize', 12, 'LineWidth', 1.0);
    hL = legend({'ISETbio (w/out S-cones):ISETbio (Banks)', 'ISETbio  (with S-cones):ISETbio (Banks)'}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity ratio', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','linear', 'FontSize', 16);
    xlim([1 100]); ylim([0.5 1]);
    grid on; box off;
    drawnow;

    function [d, the_rParams] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs)
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

function addBanksEtAlReferenceLines(rParams, plotLuminanceLine)
    luminanceColors = [0 1 0; 0 0 1; 1 0 0; 0 0 0];
    luminanceColors = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
    % Add unshifted version for reference
    banksFactor = 1;
    [A,B,C,D,E] = LoadDigitizedBanksFigure2;
    if (~rParams.oiParams.blur && ~rParams.mosaicParams.apertureBlur)
        plot(A(:,1),A(:,2),'k:','LineWidth',0.5);
        plot(A(:,1),A(:,2)*banksFactor,'r-','LineWidth',1.5);
    elseif (~rParams.oiParams.blur)
        plot(B(:,1),B(:,2),'k:','LineWidth',0.5);
        plot(B(:,1),B(:,2)*banksFactor,'r','LineWidth',1.5);
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

