function figGenerateEmployedPSFs()

    wavefrontSpatialSamples = 261*2+1;
    calcPupilDiameterMM = 3;
    umPerDegree = 300;
    psfRange = 2.5;
    showTranslation = false;
    
    thibosModels = availableCustomWvfOpticsModels();
    opticsModels{1} = 'Geisler';
    legends = {'Geisler', 'best subject', 'subject A', 'subject B', 'subject C', 'subject D'};
    for k = 6:10
        opticsModels{numel(opticsModels)+1} = thibosModels{k};
    end
    opticsModels
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 3, ...
                'colsNum', 2, ...
                'heightMargin', 0.03, ...
                'widthMargin', 0.02, ...
                'leftMargin', 0.01, ...
                'rightMargin', 0.001, ...
                'bottomMargin', 0.04, ...
                'topMargin', 0.02);
            
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'PSFS');

    
    
    
    k = 0;
 
    for oiIndex = 1:numel(opticsModels)
        
        k = k + 1;
        col = mod(k-1,2)+1;
        row = floor((k-1)/2)+1;
        
        if (strcmp(opticsModels{oiIndex}, 'Geisler'))
            pupilDiamMm = 3;
            umPerDegree = 300;
            theOI = oiCreate('wvf human', pupilDiamMm,[],[], umPerDegree);
            theCustomOI = ptb.oiSetPtbOptics(theOI,'opticsModel','Geisler');
        else
            
            [theCustomOI, ~,~] = oiWithCustomOptics(opticsModels{oiIndex}, ...
             wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, ...
                'centeringWavelength', 550, ...  %  center PSF at 550
                'showTranslation',showTranslation);
        end
        
        optics = oiGet(theCustomOI, 'optics');
        wavelengths = opticsGet(optics, 'wave');
        wavelengthsToShow = find(wavelengths == 550);
        targetWavelength = wavelengths(wavelengthsToShow);

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
        
        cmap = brewermap(1024, 'greys');
        colormap(cmap)
        
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
            axis 'equal';
            set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'FontSize', 14);
            set(gca, 'XTick', -10:1:10, 'YTick', -10:1:10);
            grid on; box on;
            
            set(gca, 'YTickLabel', {});
            
            opticsName = sprintf('%s', legends{oiIndex});
            
            if (row < 3)
                set(gca, 'XTickLabel', {});
            else
                xlabel('arc min');
            end
            
            formatFigureForPaper(hFig, ...
                'figureType', 'MOSAICS', ...
                'theAxes', gca, ...
                'theFigureTitle', opticsName);
        
            drawnow;
       
    end
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('PSFsEmployed.pdf'));
    NicePlot.exportFigToPDF(figureName, hFig, 300);

end