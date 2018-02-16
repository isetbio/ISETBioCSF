function visualizePSFfromOI(theOI, micronsPerDegree)
    optics = oiGet(theOI, 'optics');
    waves = opticsGet(optics, 'wave');
    
    psfSupportMicrons = opticsGet(optics,'psf support','um');
    xGridMinutes = 60*psfSupportMicrons{1}/micronsPerDegree;
    yGridMinutes = 60*psfSupportMicrons{2}/micronsPerDegree;
    xSupportMicrons = xGridMinutes(1,:);
    ySupportMicrons = yGridMinutes(:,1);

    hFig = figure(10); clf;
    set(hFig, 'Position', [10 10 1500 940], 'Color', [1 1 1]);
    
    rows = 3; cols = 5;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', rows, ...
               'colsNum', cols, ...
               'heightMargin',   0.04, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.03, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.05, ...
               'topMargin',      0.02);
    
    
    psfRange = 7;
    xx = find(abs(xSupportMicrons) <= psfRange);
    yy = find(abs(ySupportMicrons) <= psfRange);

    for waveIndex = 1:15
        row = floor((waveIndex-1)/cols)+1;
        col = mod(waveIndex-1, cols)+1;
        targetWavelength = waves(waveIndex*2);
        psf = opticsGet(optics,'psf data',targetWavelength );
        subplot('Position', subplotPosVectors(row,col).v);
        imagesc(xSupportMicrons(xx), ySupportMicrons(yy), psf(yy,xx)/max(psf(:)));
        hold on;
        plot(xSupportMicrons, xSupportMicrons*0, 'r-', 'Color', [0 0.0 0.0], 'LineWidth', 1.0);
        plot(xSupportMicrons*0, xSupportMicrons, 'r-', 'Color', [0 0.0 0.0], 'LineWidth', 1.0);
        hold off
        axis 'image'; axis 'xy';
        set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'CLim', [0 1], ...
            'XTick', [-10:2:10], 'YTick', [-10:2:10]);
        set(gca, 'FontSize', 16);
        if (row == 3) && (col == 1)
            xlabel('arc min');
            ylabel('arc min');
        else
            set(gca, 'XTickLabel', {}, 'YTickLabel',{});
        end
        
        grid on; box on;
        title(sprintf('%2.0f nm', targetWavelength));   
    end
    colormap(1-gray(1024));
    NicePlot.exportFigToPNG('PSF.png', hFig, 300);
end
