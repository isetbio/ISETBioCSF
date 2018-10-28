function run_OpticsVaryConditions
% This is the script used to assess the impact of different optics models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 2;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Whether to visualize the employed PSFs
    makePSFfigure = ~true;
    visualizedWavelengths = 550;
        
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; % 'ISETbioHexEccBasedLMSrealistic';
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
        'Geisler PSF (A)' ...
        'Subject 1 PSF (B)' ...
        'Subject 2 PSF (C)' ...
        'Subject 3 PSF (D)' ...
        'Subject 4 PSF (E)' ...
        'Subject 5 PSF (F)' ...
    };

    idx = [1:6];
    examinedOpticsModels = {examinedOpticsModels{idx}};
    examinedOpticsModelLegends = {examinedOpticsModelLegends{idx}};
     
    params.coneContrastDirection = 'L+M+S';

    % Denote new mosaics with higher density
    params.mosaicRotationDegs = 360;
        
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
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
        
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for modelIndex = 1:numel(examinedOpticsModels)
        params.opticsModel = examinedOpticsModels{modelIndex};
        [~,~, theFigData{modelIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Optics';
        theRatioLims = [0.4 10];
        theRatioTicks = [0.5 1 2 5 10];
%         generateFigureForPaper(theFigData, examinedOpticsModelLegends, variedParamName, mosaicName, ...
%             'figureType', 'CSF', ...
%             'inGraphText', ' A ', ...
%             'plotFirstConditionInGray', true, ...
%             'plotRatiosOfOtherConditionsToFirst', true, ...
%             'theRatioLims', theRatioLims, ...
%             'theRatioTicks', theRatioTicks ...
%             );
        generateFigureForPaper(theFigData, examinedOpticsModelLegends, variedParamName, mosaicName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' G ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    
    if (makePSFfigure)
        for iw = 1:numel(visualizedWavelengths)
            generatePSFsFigure(examinedOpticsModels, examinedOpticsModelLegends, visualizedWavelengths(iw));
        end
    end
end

function generatePSFsFigure(examinedOpticsModels, examinedOpticsModelLegends, visualizedWavelength)
    wavefrontSpatialSamples = 261*2+1;
    calcPupilDiameterMM = 3;
    umPerDegree = 300;
    psfRange = 2.5;
    showTranslation = false;
    
    rowsNum = 3;
    colsNum = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', rowsNum, ...
                'colsNum', colsNum, ...
                'heightMargin', 0.015, ...
                'widthMargin', 0.02, ...
                'leftMargin', 0.01, ...
                'rightMargin', 0.001, ...
                'bottomMargin', 0.05, ...
                'topMargin', 0.001);
            
    hFigPSF = figure(2); clf; 
    formatFigureForPaper(hFigPSF, 'figureType', 'PSFS');
    
    for oiIndex = 1:numel(examinedOpticsModels)    
        switch (examinedOpticsModelLegends{oiIndex})
            case 'Geisler PSF (A)'
                prefix = ' A ';
            case 'Subject 1 PSF (B)' 
                prefix = ' B ';
            case 'Subject 2 PSF (C)'
                prefix = ' C ';
            case 'Subject 3 PSF (F)'
                prefix = ' D ';
            case 'Subject 4 PSF (E)'
                prefix = ' E ';
            case 'Subject 5 PSF (F)'
                prefix = ' F ';
        end
        
        col = mod(oiIndex-1,colsNum)+1;
        row = floor((oiIndex-1)/colsNum)+1;
        
        switch (examinedOpticsModels{oiIndex})
            case 'Geisler'
                pupilDiamMm = 3;
                umPerDegree = 300;
                theOI = oiCreate('wvf human', pupilDiamMm,[],[], umPerDegree);
                theCustomOI = ptb.oiSetPtbOptics(theOI,'opticsModel','Geisler');
                
            otherwise
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
        contourLevels = 0:0.05:0.95;
        contourf(xMinutes, yMinutes,wavePSF/max(wavePSF(:)), contourLevels);
        hold on;
        plot([0 0], [xMinutes(1) xMinutes(end)], 'r-', 'LineWidth', 1.0);
        plot([xMinutes(1) xMinutes(end)], [0 0], 'r-', 'LineWidth', 1.0);
        hold off;
        set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'FontSize', 14);
        set(gca, 'XTick', -10:1:10, 'YTick', -10:1:10);
        set(gca, 'YTickLabel', {});
        if (row < 3)
             %set(gca, 'XTickLabel', {});
        else
             xlabel('\it position (arc min)', 'FontWeight', 'normal', 'FontSize', 24);
        end
        
        inGraphTextPos = [-2.4 2.0];
        t = text(gca, inGraphTextPos(1), inGraphTextPos(2), sprintf('%s', prefix));
        
        formatFigureForPaper(hFigPSF, ...
                'figureType', 'PSFS', ...
                'theAxes', gca, ...
                'theText', t, ...
                'theTextFontSize', 32, ...
                'theFigureTitle', '');
               
    end
    
    figure(hFigPSF)
    cmap = brewermap(1024, 'greys');
    colormap(cmap);
    drawnow;
        
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('PSFsEmployed%2.0fnm.pdf', targetWavelength));
    NicePlot.exportFigToPDF(figureName, hFigPSF, 300);
end

