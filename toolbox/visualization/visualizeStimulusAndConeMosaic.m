function visualizeStimulusAndConeMosaic(theMosaic, oiModulated)
    
    % Get the cone locs in degrees
    if isa(theMosaic, 'coneMosaicHex')
        coneLocsInMeters = theMosaic.coneLocsHexGrid;
    else
        coneLocsInMeters = theMosaic.coneLocs;
    end
    mosaicFOV = theMosaic.fov;
    coneLocsInDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width  * theMosaic.fov(1);
    coneLocsInDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    
    support = oiGet(oiModulated, 'spatial support', 'microns');
    xaxis = support(1,:,1);
    yaxis = support(:,1,2);
    XYZ = oiGet(oiModulated, 'xyz');
    XYZ = XYZ / max(XYZ(:));
    rgbImage = xyz2srgb(XYZ);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 2, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.01);
       
    hFig = figure(98765); clf;
    set(hFig, 'Position', [10 10 1200 550]);
            
    subplot('Position', subplotPosVectors(1,1).v);
    imagesc(xaxis, yaxis, rgbImage, [0 1]);
    axis 'xy';
    axis 'image';
    hold on;
    % outline mosaic extent in green
    x = mosaicFOV(1) * [-0.5 0.5 0.5 -0.5 -0.5];
    y = mosaicFOV(2) * [-0.5 -0.5 0.5 0.5 -0.5];
    plot(x,y, 'g-', 'LineWidth', 1.5);
    hold off
    set(gca, 'XLim', mosaicFOV(1)/2*[-1 1]*1.02, 'YLim', mosaicFOV(2)/2*[-1 1]*1.02);
    set(gca, 'FontSize', 14);
    ylabel('degrees');
    title('stimulus and cone mosaic (stimulus view)');

    subplot('Position', subplotPosVectors(1,2).v);
    imagesc(xaxis, yaxis, rgbImage, [0 1]);
    axis 'xy';
    axis 'image'
    hold on;
    plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
    hold off;
    set(gca, 'XLim', mosaicFOV(1)/2*[-1 1]*1.02, 'YLim', mosaicFOV(2)/2*[-1 1]*1.02);
    set(gca, 'FontSize', 14);
    xlabel('degrees');
    title('stimulus and cone mosaic (cone mosaic view)', 'FontSize', 16);

    colormap(gray(1024));
    drawnow;
    
    NicePlot.exportFigToPDF('dd', hFig, 300);
    fprintf('Figure exported to dd\n');
end


