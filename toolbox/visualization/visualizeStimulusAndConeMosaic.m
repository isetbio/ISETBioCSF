function visualizeStimulusAndConeMosaic(theMosaic, thePeakOI, paramsList)
    
    % Get the cone locs in degrees
    if isa(theMosaic, 'coneMosaicHex')
        coneLocsInMeters = theMosaic.coneLocsHexGrid;
    else
        coneLocsInMeters = theMosaic.coneLocs;
    end
    mosaicFOV = theMosaic.fov;
    coneLocsInDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width  * theMosaic.fov(1);
    coneLocsInDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    
    support = oiGet(thePeakOI, 'spatial support', 'microns');
    micronsPerDegree = oiGet(thePeakOI, 'width')*1e6 / oiGet(thePeakOI, 'hfov'); % 290;

    xaxis = support(1,:,1)/micronsPerDegree;
    yaxis = support(:,1,2)/micronsPerDegree;
    illumMap = oiCalculateIlluminance(thePeakOI);
    minIllum = min(illumMap(:));
    maxIllum = max(illumMap(:));
    if (minIllum == maxIllum)
        illumRange = [minIllum*0.99 maxIllum*1.01];
    else
        illumRange = [minIllum maxIllum];
    end
    illumMap = 0.5 + 0.5*(illumMap-illumRange(1))/(illumRange(2)-illumRange(1));
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 1, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.07, ...
           'topMargin',      0.04);
       
    hFig = figure(98765); clf;
    set(hFig, 'Position', [10 10 700 780], 'Color', [1 1 1]);
            
    subplot('Position', subplotPosVectors(1,1).v);
    imagesc(xaxis, yaxis, illumMap);
    axis 'xy';
    axis 'image'
    hold on;
    plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'k.', 'MarkerSize', 8);
    hold off;
    set(gca, 'CLim', [0 1], 'YTickLabel', {}, 'XLim', mosaicFOV(1)/2*[-1 1]*1.1, 'YLim', mosaicFOV(2)/2*[-1 1]*1.1);
    set(gca, 'FontSize', 18);
    xlabel('degrees');
    title(sprintf('retinal stimulus (illuminance: %f - %f)\nand cone locations', minIllum, maxIllum), 'FontSize', 16);
    colormap(gray(1024));
    drawnow;
    
    % Save figure
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('SpatialScheme');
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');

end


