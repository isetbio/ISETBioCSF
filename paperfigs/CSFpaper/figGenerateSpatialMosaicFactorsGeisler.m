function figGenerateSpatialMosaicFactors()
    
    % 'WvfHuman', 'Geisler'
    opticsModel = 'Geisler'; % 'WvfHuman';
    
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
    
    % Retrieve reference mosaic ('originalBanks') CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, referenceMosaicLegend] = paramsForComparativeMosaicAnalysis(mosaicRef);
    [d, rParamsRef] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
    sfRef = d.cyclesPerDegree;
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        referenceCSF(luminanceIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    
    % Retrieve 'fullIsetbioWithScones' CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, mosaicTest1Legend] = paramsForComparativeMosaicAnalysis(mosaicTest1);
    [d, rParamsTest1] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
    test1Color = [0.8 0.3 0.9];
    
    sfTest1 = d.cyclesPerDegree;
    if (any(sfTest1-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF1(luminanceIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    
    if (1==2)
    % Retrieve 'defaultIsetbio' CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, mosaicTest2Legend] = paramsForComparativeMosaicAnalysis(mosaicTest2);
    [d, rParamsTest2] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);

    
    sfTest2 = d.cyclesPerDegree;
    if (any(sfTest2-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF2(luminanceIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    end
    
    % Compute CSF ratios
    ratioCSF1 = testCSF1 ./ referenceCSF;
    if (1==2)
    ratioCSF2 = testCSF2 ./ referenceCSF;
    end
    
    test2Color = [0.3 0.7 0.99];
    
    hFig = figure(1); clf;
    set(hFig ,'Position',[100 100 1000 1140], 'Color', [1 1 1]);
    subplot('Position', [0.055 0.30 0.41 0.69]);
    hold on
    addBanksEtAlReferenceLines(rParamsTest1, plotLuminanceLines);
    plot(sfRef, referenceCSF, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    plot(sfTest1, testCSF1, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', test2Color, 'MarkerSize', 12);
    %plot(sfTest2, testCSF2, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', test2Color, 'MarkerSize', 12);
    hL = legend({'Banks et al study', referenceMosaicLegend, mosaicTest1Legend}, 'Location', 'NorthOutside');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16);
    xlabel('', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontSize' ,16, 'FontWeight', 'bold');
    xlim([1 100]); ylim([1 3000]);
    grid on; box off;
    
    subplot('Position', [0.055 0.055 0.41 0.20]);
    hold on
    plot(sfRef, ratioCSF1, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', test2Color, 'MarkerSize', 12, 'LineWidth', 1.0);
    %plot(sfRef, ratioCSF2, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', test2Color, 'MarkerSize', 12, 'LineWidth', 1.0);
    hL = legend({mosaicTest1Legend}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity ratio', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','linear', 'FontSize', 16);
    xlim([1 100]); ylim([0.5 1]);
    grid on; box off;
    
    sfIndex = 3;
    subplot('Position', [0.50 0.53 0.48 0.48]);
    plotMosaic(rParamsRef, sfIndex, sprintf('%s', referenceMosaicLegend), true); 
    subplot('Position', [0.50 0.025 0.48 0.48]);
    plotMosaic(rParamsTest1, sfIndex, sprintf('%s', mosaicTest1Legend), true);

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

    function noiseFreeResponses = getNoiseFreeResponses(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, cyclesPerDeg, conditionsToVisualize)
        [~, ~, noiseFreeResponses] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', opticsModel, ...
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
            'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter), ...
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
            'visualizeTransformedSignals', visualizeTransformedSignals, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier, ...
            'useRBFSVMKernel', useRBFSVMKernel, ...
            'performanceTrialsUsed', nTrainingSamples, ...
            'spatialPoolingKernelParams', spatialPoolingKernelParams ...
            );
    end

end

function plotMosaic(rParamsAllConds, sfIndex, titleString, showXlabel)
    lumIndex = 1;
    rParams = rParamsAllConds{sfIndex,lumIndex};
    fprintf('Loading a previously saved cone mosaic\n');
    coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    theMosaic.visualizeGrid('axesHandle', gca);
    conesNum = numel(find(theMosaic.pattern(:) > 1));
    title(gca, sprintf('%s (cones: %d)', titleString, conesNum), 'FontWeight', 'normal', 'FontWeight', 'bold', 'FontSize', 16);
    if (showXlabel)
        xlabel(gca, sprintf('%2.0f microns', theMosaic.width*1e6), 'FontWeight', 'Bold', 'FontSize', 18);
    end
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



