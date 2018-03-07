function visualizePupilFunction(pupilFunctionData, paramsList)

    nWaves = numel(pupilFunctionData);
    
       
    for iWave = 1:1:nWaves
        p = pupilFunctionData{iWave};
        thisWave = p.wavelength;
        
        hFig = figure(iWave); clf;
        set(hFig, 'Position', [10 10 1000 350], 'Color', [1 1 1]);
        
        ax1 = subplot(1,2,1);
        imagesc(p.pupilPos, p.pupilPos, wrapToPi(angle(p.pupilfunc)));
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
        
        ax2 = subplot(1,2,2);
        imagesc(p.pupilPos, p.pupilPos, abs(p.pupilfunc));
        axis 'image'; axis 'xy';
        set(gca, 'XLim', 1.5*[-1 1], 'YLim', 1.5*[-1 1], 'FontSize', 14);
        xlabel('pupil plane (mm)', 'FontWeight', 'bold');
        box on; grid on;
        colormap(ax2, gray(1024));
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

