function run_OpticsVaryConditions
% This is the script used to assess the impact of different optics models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealistic';
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    % Optics models we tested
    examinedOpticsModels = {...
        'Geisler' ...
        'ThibosBestPSFSubject3MMPupil' ...
        'ThibosDefaultSubject3MMPupil' ...
        'ThibosAverageSubject3MMPupil' ...
        'ThibosDefocusedSubject3MMPupil' ...
        'ThibosVeryDefocusedSubject3MMPupil' ...
    };
    examinedOpticsModelLegends = {...
        'Geisler PSF' ...
        'wvf-based PSF (best subject)' ...
        'wvf-based PSF (subject A)' ...
        'wvf-based PSF (subject B)' ...
        'wvf-based PSF (subject C)' ...
        'wvf-based PSF (subject D)' ...
    };

    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
        
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for modelIndex = 1:numel(examinedOpticsModels)
        params.opticsModel = examinedOpticsModels{modelIndex};
        [~,~, theFigData{modelIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Optics';
        theRatioLims = [0.5 12];
        theRatioTicks = [0.5 1 2 5 10];
        generateFigureForPaper(theFigData, examinedOpticsModelLegends, variedParamName, mosaicName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
        
        generatePSFsFigure(examinedOpticsModels);
    end
end

function generatePSFsFigure(examinedOpticsModels)
    wavefrontSpatialSamples = 261*2+1;
    calcPupilDiameterMM = 3;
    umPerDegree = 300;
    psfRange = 2.5;
    showTranslation = false;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 3, ...
                'colsNum', 2, ...
                'heightMargin', 0.015, ...
                'widthMargin', 0.02, ...
                'leftMargin', 0.01, ...
                'rightMargin', 0.001, ...
                'bottomMargin', 0.05, ...
                'topMargin', 0.001);
            
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'PSFS');

    prefices = {' B ', ' C ', ' D ', ' E ', ' F ', ' G '};
    for oiIndex = 1:numel(examinedOpticsModels)
        col = mod(oiIndex-1,2)+1;
        row = floor((oiIndex-1)/2)+1;
        
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
        targetWavelength = wavelengths(find(wavelengths == 550));

        xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(xSfCyclesDeg,ySfCyclesDeg);
                
        waveOTF = opticsGet(optics,'otf data',targetWavelength);
            [xGridMinutes, yGridMinutes, wavePSF] = OtfToPsf(...
                xSfGridCyclesDeg,ySfGridCyclesDeg,fftshift(waveOTF), ...
                'warningInsteadOfErrorForNegativeValuedPSF', true ...
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
        contourLevels = 0:0.1:1.0;
        contourf(xMinutes, yMinutes,wavePSF/max(wavePSF(:)), contourLevels);
        hold on;
        plot([0 0], [xMinutes(1) xMinutes(end)], 'r-', 'LineWidth', 1.0);
        plot([xMinutes(1) xMinutes(end)], [0 0], 'r-', 'LineWidth', 1.0);
        hold off;
        set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'FontSize', 14);
        set(gca, 'XTick', -10:1:10, 'YTick', -10:1:10);
        set(gca, 'YTickLabel', {});
        if (row < 3)
            set(gca, 'XTickLabel', {});
        else
            xlabel('retinal position (arc min)', 'FontWeight', 'bold');
        end
            
        inGraphTextPos = [-2.3 2.0];
        t = text(gca, inGraphTextPos(1), inGraphTextPos(2), prefices{oiIndex});
        
        formatFigureForPaper(hFig, ...
                'figureType', 'PSFS', ...
                'theAxes', gca, ...
                'theText', t, ...
                'theFigureTitle', '');
    end
    
    cmap = brewermap(1024, 'greys');
    colormap(cmap);
    drawnow;
        
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('PSFsEmployed.pdf'));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

