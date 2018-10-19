function run_GratingOrientationVaryConditions
% This is the script used to assess the impact of different grating orientations for
% the typical subject PSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
     
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Optics to use
    params.opticsModel = 'ThibosAverageSubject3MMPupil';
    
    % Optical image padding. If the FOV of the scene is less than
    % opticalImagePadSizeDegs, we will pad the oi so that its size
    % is equal to opticalImagePadSizeDegs.
    params.opticalImagePadSizeDegs = 0.5;
    
    params.cyclesPerDegreeExamined = [4 8 16 32 50 60];
    
    % Grating orientations to examine
    delta = 22.5;
    examinedOrientations = [0 delta 45 45+delta 90 90+delta 135 135+delta];
    examinedOrientationLegends = {...
        ' 0 deg' ...
        '22 deg' ...
        '45 deg' ...
        '67 deg' ...
        '90 deg' ...
        '112 deg' ...
        '135 deg' ...
        '157 deg' ...
    };
    
    params.coneContrastDirection = 'L+M+S';
    
    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;

    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;
        
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = true;
    params.visualizeDisplay = ~true;
        
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for oriIndex = 1:numel(examinedOrientations)
        params.stimulusOrientationDegs = examinedOrientations(oriIndex);
        [~,~, theFigData{oriIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'GratingOrientation';
        theRatioLims = [0.5 2];
        theRatioTicks = [0.5 1 2 5 10];
%         generateFigureForPaper(theFigData, examinedOpticsModelLegends, variedParamName, mosaicName, ...
%             'figureType', 'CSF', ...
%             'inGraphText', ' A ', ...
%             'plotFirstConditionInGray', true, ...
%             'plotRatiosOfOtherConditionsToFirst', true, ...
%             'theRatioLims', theRatioLims, ...
%             'theRatioTicks', theRatioTicks ...
%             );
        generateFigureForPaper(theFigData, examinedOrientationLegends, variedParamName, mosaicName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' G ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
        
        generateOrientationCSF(theFigData, examinedOrientations);
    end
end

function generateOrientationCSF(theFigData, examinedOrientations)
    for orientationIndex = 1:numel(examinedOrientations)
        figData = theFigData{orientationIndex};
        lumIndex = 1;
        referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
        spatialFrequency = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
        thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
        contrastSensitivity(orientationIndex,:) = 1./(thresholdContrasts*referenceContrast);
    end
    
    contrastSensitivity = log10(contrastSensitivity);
    sfColor = brewermap(numel(spatialFrequency), '*Set1');
    
    figure(111); clf;
    for sfIndex = 1:numel(spatialFrequency)
        csf = contrastSensitivity(:, sfIndex);
        %csf = csf / max(csf);
        x = zeros(1,numel(examinedOrientations)*2+1);
        y = x;
        for oriIndex = 1:numel(examinedOrientations)
            theta = examinedOrientations(oriIndex)/180*pi;
            x(oriIndex) = cos(theta)*csf(oriIndex);
            y(oriIndex) = sin(theta)*csf(oriIndex);
            x(oriIndex+numel(examinedOrientations)) = -x(oriIndex);
            y(oriIndex+numel(examinedOrientations)) = -y(oriIndex);
        end
        x(numel(examinedOrientations)*2+1) = x(1);
        y(numel(examinedOrientations)*2+1) = y(1);
        edgeColor = squeeze(sfColor(sfIndex,:))*0.5;
        faceColor = squeeze(sfColor(sfIndex,:));
        plot(x,y, 'ko-', 'Color', edgeColor, 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', faceColor);
        if (sfIndex == 1)
            hold on;
            sfLegends{1} = sprintf('%2.0f c/deg', spatialFrequency(sfIndex));
        else
            sfLegends{numel(sfLegends)+1} = sprintf('%2.0f c/deg', spatialFrequency(sfIndex));
        end
    end
    set(gca, 'XLim', [-3.5 3.5], 'YLim', [-3.5 3.5]);
    set(gca, 'XTick', [-3.5:0.5:3.5], 'YTick', [-3.5:0.5:3.5]);
    hL = legend(sfLegends);
    axis 'square';
    grid on
end
