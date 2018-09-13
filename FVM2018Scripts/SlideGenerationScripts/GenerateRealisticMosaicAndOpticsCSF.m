function GenerateRealisticMosaicAndOpticsCSF 

    % Script to generate the slide with the Banks'87 ideal and human observer data
    
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = ~true;
    makePSFfigure = true;
    
    computeResponses = false;
    findPerformance = false;
    
    % Init condition index
    condIndex = 0;
    
    
    % ISETBio simulation using Banks conditions
    condIndex = condIndex+1;
    examinedConds(condIndex).conditionLabel = 'Banks mosaic/optics, ideal observer';
    examinedConds(condIndex).mosaicName = 'originalBanks';
    examinedConds(condIndex).opticsModel = 'Geisler';
    examinedConds(condIndex).inferenceEngine = 'mlpt';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).emPathType = 'frozen0';
    examinedConds(condIndex).centeredEMPaths = ~true;
    examinedConds(condIndex).frameRate = 10;
    examinedConds(condIndex).responseStabilizationMilliseconds =  40;
    examinedConds(condIndex).responseExtinctionMilliseconds = 40;
    
    % Our best estimate of mosaic + optics, MLPT inference engine
    condIndex = condIndex+1;                                
    examinedConds(condIndex).conditionLabel = 'Realistic mosaic/optics, ideal observer';
    examinedConds(condIndex).mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
    examinedConds(condIndex).opticsModel = 'ThibosAverageSubject3MMPupil';
    examinedConds(condIndex).inferenceEngine = 'mlpt';
    examinedConds(condIndex).signal = 'isomerizations';
    examinedConds(condIndex).emPathType = 'frozen0';
    examinedConds(condIndex).centeredEMPaths = ~true;
    examinedConds(condIndex).frameRate = 10;
    examinedConds(condIndex).responseStabilizationMilliseconds = 40;
    examinedConds(condIndex).responseExtinctionMilliseconds = 40;
    
    
    
    % Go
    [examinedLegends, theFigData, pupilDiamMm] = runConditions(...
        examinedConds, computationInstance, computeResponses, findPerformance);
    
    if (makeSummaryFigure)
        variedParamName = 'RealisticMosaicAndOptics';
        theRatioLims = [0.02 2.0];
        theRatioTicks = [0.05  0.1 0.2 0.5 1.0];
        formatLabel = '_vs_BanksMosaicAndOptics'; 
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'showSubjectData', (pupilDiamMm == 2), ...
            'showSubjectMeanData', false, ...
            'showLegend', ~true, ...
            'plotFirstConditionInGray', true, ...
            'showBanksPaperIOAcurves', ~true, ...
            'showOnly23CDM2IOAcurve', ~true, ...
            'plotRatiosOfOtherConditionsToFirst', false, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    
    if (makePSFfigure)
        for k = 1:numel(examinedConds)
            examinedOpticsModels{k} = examinedConds(k).opticsModel;
            examinedOpticsModelLegends{k} = '';
        end

        visualizedWavelengths = [500 550 600];
        for iw = 1:numel(visualizedWavelengths)
            generatePSFsFigure(examinedOpticsModels, examinedOpticsModelLegends, visualizedWavelengths(iw));
        end
    end
    
end

function [examinedLegends, theFigData, pupilDiamMm] = runConditions(examinedConds, computationInstance, computeResponses, findPerformance)
    examinedLegends = {};
    theFigData = {};
    
    for condIndex = 1:numel(examinedConds)
        cond = examinedConds(condIndex);
        mosaicName = cond.mosaicName;
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        
        % Update params
        params.opticsModel = cond.opticsModel;
        params.performanceClassifier = cond.inferenceEngine;
        params.performanceSignal = cond.signal;
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMPaths;
        params.frameRate = cond.frameRate;
        params.responseStabilizationMilliseconds = cond.responseStabilizationMilliseconds;
        params.responseExtinctionMilliseconds = cond.responseExtinctionMilliseconds;
        
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        
        params = getRemainingDefaultParams(params, condIndex, cond.conditionLabel, computeResponses, findPerformance);  
        
        % Update returned items
        pupilDiamMm(condIndex) = params.pupilDiamMm;
        examinedLegends{numel(examinedLegends) + 1} = cond.conditionLabel;
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (any(diff(pupilDiamMm) ~= 0))
        pupilDiamMm
        error('PuliDiamMM are different for different conditions');
    else
        pupilDiamMm = pupilDiamMm(1);
    end
    
end


function params = getRemainingDefaultParams(params, condIndex, conditionLabel, computeResponses, findPerformance)
            
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    
    if contains(conditionLabel,'Banks mosaic/optics')
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50];
    else
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];
    end
    
    if (strcmp(conditionLabel, 'Banks mosaic/optics, MLPT, 3mm'))
        params.pupilDiamMm = 3.0;
    else
        params.pupilDiamMm = 2.0;
    end
    
    fprintf('>>>>>> \t %d pupil: %f\n', condIndex, params.pupilDiamMm);
    
                
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
end

function generatePSFsFigure(examinedOpticsModels, examinedOpticsModelLegends, visualizedWavelength)
    wavefrontSpatialSamples = 261*2+1;
    calcPupilDiameterMM = 3;
    umPerDegree = 300;
    psfRange = 4.0;
    showTranslation = false;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', numel(examinedOpticsModels), ...
                'colsNum', 1, ...
                'heightMargin', 0.015, ...
                'widthMargin', 0.02, ...
                'leftMargin', 0.01, ...
                'rightMargin', 0.001, ...
                'bottomMargin', 0.05, ...
                'topMargin', 0.001);
            
    hFig = figure(2); clf;
    formatFigureForPaper(hFig, 'figureType', 'PSFS2');
        
    for oiIndex = 1:numel(examinedOpticsModels)    
%         switch (examinedOpticsModelLegends{oiIndex})
%             case 'Geisler PSF (A)'
%                 prefix = '';
%             case 'Subject 1 PSF (B)' 
%                 prefix = '';
%             case 'Subject 2 PSF (c)'
%                 prefix = '';
%             case 'Subject 3 PSF (D)'
%                 prefix = '';
%             case 'Subject 4 PSF (E)'
%                 prefix = '';
%             case 'Subject 5 PSF (F)'
%                 prefix = '';
%         end
        prefix = [];
        
        col = 1;
        row = oiIndex;
        
        if (strcmp(examinedOpticsModels{oiIndex}, 'Geisler'))
            pupilDiamMm = 3;
            umPerDegree = 300;
            theOI = oiCreate('wvf human', pupilDiamMm,[],[], umPerDegree);
            theCustomOI = ptb.oiSetPtbOptics(theOI,'opticsModel','Geisler');
        else
            [theCustomOI, ~,~] = oiWithCustomOptics(examinedOpticsModels{oiIndex}, ...
             wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, ...
                'centeringWavelength', 550, ...  %  center PSF at 550
                'showTranslation',showTranslation);
        end
        
        optics = oiGet(theCustomOI, 'optics');
        wavelengths = opticsGet(optics, 'wave');
        [~,idx] = min(abs(wavelengths - visualizedWavelength));
        targetWavelength = wavelengths(idx);

        xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(xSfCyclesDeg,ySfCyclesDeg);
                
        waveOTF = opticsGet(optics,'otf data',targetWavelength);
            [xGridMinutes, yGridMinutes, wavePSF] = OtfToPsf(...
                xSfGridCyclesDeg,ySfGridCyclesDeg,fftshift(waveOTF), ...
                'warningInsteadOfErrorForNegativeValuedPSF', 1 ...
        );
        xMinutes = squeeze(xGridMinutes(1,:));
        yMinutes = squeeze(yGridMinutes(:,1));  
        xx = find(abs(xMinutes) <= psfRange);
        yy = find(abs(yMinutes) <= psfRange);
        wavePSF = wavePSF(yy,xx);
        xMinutes = xMinutes(xx);
        yMinutes = yMinutes(yy);
            
        pos = subplotPosVectors(row,col).v;   
        subplot('Position', pos);
        contourLevels = [0.05 0.1:0.1:1.0];
        normPSF = wavePSF/max(wavePSF(:));
        contourf(xMinutes, yMinutes, normPSF, contourLevels, 'LineWidth', 1.5);
        hold on;
        [~,midRow] = min(abs(yMinutes));
        plot(xMinutes, xMinutes(1) + 0.45*(xMinutes(end)-xMinutes(1))*squeeze(normPSF(midRow,:)), '-', 'Color', [0 0 1], 'LineWidth', 5.0);
        plot(xMinutes, xMinutes(1) + 0.45*(xMinutes(end)-xMinutes(1))*squeeze(normPSF(midRow,:)), '-', 'Color', [0.7 0.7 1], 'LineWidth', 3.0);
        plot(xMinutes, xMinutes(1) + 0.0*(xMinutes(end)-xMinutes(1))*squeeze(normPSF(midRow,:)), 'k-', 'LineWidth', 1.5);
        plot([0 0], [xMinutes(1) xMinutes(end)], 'k-', 'LineWidth', 1.0);
        plot([xMinutes(1) xMinutes(end)], [0 0], 'k-', 'LineWidth', 1.0);
        hold off;
        set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'FontSize', 14);
        set(gca, 'XTick', -10:1:10, 'YTick', -10:1:10);
        set(gca, 'YTickLabel', {});
        if (row < 3)
            set(gca, 'XTickLabel', {});
        else
            xlabel('retinal position (arc min)', 'FontWeight', 'bold');
        end
            
        inGraphTextPos = [-2.4 2.0];
        if (~isempty(prefix))
            t = text(gca, inGraphTextPos(1), inGraphTextPos(2), sprintf('%s', prefix));
        else
            t = [];
        end
        
        formatFigureForPaper(hFig, ...
                'figureType', 'PSFS2', ...
                'theAxes', gca, ...
                'theText', t, ...
                'theTextFontSize', 32, ...
                'theFigureTitle', '');
    end
    
    cmap = brewermap(1024, 'greys');
    colormap(cmap);
    drawnow;
        
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('PSFsEmployed%2.0fnm.pdf', targetWavelength));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

