function figGenerateSpatialMosaicFactors()
    
    % 'WvfHuman', 'Geisler'
    opticsModel = 'Geisler';
    
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
    mosaicRef = 'originalBanks';
    mosaicTest1 = 'defaultIsetbio';
    mosaicTest2 = 'fullIsetbioWithScones';
    
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
    
    % Retrieve 'defaultIsetbio' CSF data
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, mosaicTest2Legend] = paramsForComparativeMosaicAnalysis(mosaicTest2);
    [d, rParamsTest2] = getData(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs);
    test2Color = [0.3 0.7 0.99];
    
    sfTest2 = d.cyclesPerDegree;
    if (any(sfTest2-sfRef)~=0)
        error('sfs do not match');
    end
    for luminanceIndex = 1:size(d.mlptThresholds,1)
        testCSF2(luminanceIndex,:) = 1./([d.mlptThresholds(luminanceIndex,:).thresholdContrasts]*d.mlptThresholds(1).testConeContrasts(1));
    end
    
    
    % Compute CSF ratios
    ratioCSF1 = testCSF1 ./ referenceCSF;
    ratioCSF2 = testCSF2 ./ referenceCSF;
    
    hFig = figure(1); clf;
    set(hFig ,'Position',[100 100 950 1000], 'Color', [1 1 1]);
    subplot('Position', [0.055 0.30 0.41 0.69]);
    hold on
    addBanksEtAlReferenceLines(rParamsTest1, plotLuminanceLines);
    plot(sfRef, referenceCSF, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 12);
    plot(sfTest1, testCSF1, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', test1Color, 'MarkerSize', 12);
    plot(sfTest2, testCSF2, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', test2Color, 'MarkerSize', 12);
    hL = legend({'Banks et al study', referenceMosaicLegend, mosaicTest1Legend, mosaicTest2Legend}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    set(gca,'XScale','log','YScale','log', 'FontSize', 16);
    xlabel('', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontSize' ,16, 'FontWeight', 'bold');
    xlim([1 100]); ylim([1 5000]);
    grid on; box off;
    
    subplot('Position', [0.055 0.055 0.41 0.20]);
    hold on
    plot(sfRef, ratioCSF1, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', test1Color, 'MarkerSize', 12, 'LineWidth', 1.0);
    plot(sfRef, ratioCSF2, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', test2Color, 'MarkerSize', 12, 'LineWidth', 1.0);
    hL = legend({mosaicTest1Legend, mosaicTest2Legend}, 'Location', 'NorthEast');
    set(hL, 'FontSize', 13, 'FontName', 'Menlo');
    xlabel('spatial frequency (cpd)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('contrast sensitivity ratio', 'FontSize' ,16, 'FontWeight', 'bold');
    set(gca,'XScale','log','YScale','linear', 'FontSize', 16);
    xlim([1 100]); ylim([0.5 1]);
    grid on; box off;
    
    sfIndex = 3;
    subplot('Position', [0.50 0.51 0.48 0.48]);
    plotMosaic(rParamsRef, sfIndex, sprintf('%s', referenceMosaicLegend)); 
    subplot('Position', [0.50 0.01 0.48 0.48]);
    plotMosaic(rParamsTest2, sfIndex, sprintf('%s', mosaicTest2Legend));
    
    
    conditionsToVisualize = [18];
    cyclesPerDeg = cyclesPerDegreeExamined(sfIndex);
    [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, ~] = paramsForComparativeMosaicAnalysis(mosaicTest2);
    noiseFreeResponses = getNoiseFreeResponses(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, cyclesPerDeg, conditionsToVisualize);
    noiseFreeIsomerizations = noiseFreeResponses{1,1}.stimData{conditionsToVisualize}.noiseFreeIsomerizations;
    [~,idx] = max(noiseFreeIsomerizations(:));
    [~, timeBinOfMaxResponse] = ind2sub(size(noiseFreeIsomerizations), idx);
    noiseFreeIsomerizations = noiseFreeIsomerizations(:, timeBinOfMaxResponse);
   
    figure(hFig);
    theAxes = axes('Position',[0.06,0.305,0.2,0.2]);
    plotMosaicActivation(theAxes, rParamsTest2, sfIndex, noiseFreeIsomerizations, 'isomerizations', cyclesPerDeg);

    sfIndex = sfIndex + 1;
    cyclesPerDeg = cyclesPerDegreeExamined(sfIndex);
    noiseFreeResponses = getNoiseFreeResponses(coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, cyclesPerDeg, conditionsToVisualize);
    noiseFreeIsomerizations = noiseFreeResponses{1,1}.stimData{conditionsToVisualize}.noiseFreeIsomerizations;
    [~,idx] = max(noiseFreeIsomerizations(:));
    [~, timeBinOfMaxResponse] = ind2sub(size(noiseFreeIsomerizations), idx);
    noiseFreeIsomerizations = noiseFreeIsomerizations(:, timeBinOfMaxResponse);
    
    figure(hFig);
    theAxes = axes('Position',[0.26,0.305,0.2,0.2]);
    plotMosaicActivation(theAxes, rParamsTest2, sfIndex, noiseFreeIsomerizations, 'isomerizations',cyclesPerDeg);

    
    
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

function plotMosaic(rParamsAllConds, sfIndex, titleString)
    lumIndex = 1;
    rParams = rParamsAllConds{sfIndex,lumIndex};
    fprintf('Loading a previously saved cone mosaic\n');
    coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    theMosaic.visualizeGrid('axesHandle', gca);
    title(gca, sprintf('%s', titleString), 'FontWeight', 'normal', 'FontWeight', 'bold', 'FontSize', 16);
    conesNum = numel(find(theMosaic.pattern(:) > 1));
    %xlabel(gca, sprintf('%2.0f microns, cones: %d', theMosaic.width*1e6, conesNum), 'FontWeight', 'Bold', 'FontSize', 18);
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
