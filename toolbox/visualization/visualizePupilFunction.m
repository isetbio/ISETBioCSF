function visualizePupilFunction(visualizePupilFunction, paramsList)

    pupilFunction = wvfGet(theWVF, 'pupil function', targetWavelength);
    wavelengths = wvfGet(theWVF, 'calcwavelengths');
    nWaves = numel(wavelengths);
    
       
    for iWave = 1:1:nWaves
        targetWavelength = wavelengths(iWave);
        pupilFunction = wvfGet(theWVF, 'pupil function', targetWavelength);
        pupilSupport = wvfGet(theWVF, 'pupil spatial samples', 'mm', targetWavelength);

        
        hFig = figure(iWave); clf;
        set(hFig, 'Position', [10 10 1000 350], 'Color', [1 1 1]);
        
        ax1 = subplot(1,2,1);
        imagesc(pupilSupport, pupilSupport, angle(pupilFunction));
        numel(p.pupilPos)
        colormap(ax1, hsv(1024));
        axis 'image'; axis 'xy';
        set(gca, 'CLim', [-pi pi], 'XLim', 1.5*[-1 1], 'YLim', 1.5*[-1 1], 'FontSize', 14);
        phaseTicks = -pi : pi/3: pi;
        phaseTickLabels = sprintf('%2.0f deg\n', phaseTicks/pi*180);
        colorbar('Ticks', phaseTicks, 'TickLabels', phaseTickLabels);
        xlabel('pupil plane (mm)', 'FontWeight', 'bold');
        box on; grid on;
        title(sprintf('%2.0fnm', thisWave));
        
        colorbar
        drawnow;
    end
    
    if (~isempty(paramsList))
        % Export to PDF
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('PupilFunction');
        rwObject.write(fileName, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPNG', 'FigureHandle', hFig, 'FigureType', 'png');
    end
    
end

