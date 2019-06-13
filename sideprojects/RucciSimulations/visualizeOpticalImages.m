function visualizeOpticalImages(nullSceneOI, lowFrequencyOIs, highFrequencyOIs, lowFrequencyOIsOrtho, highFrequencyOIsOrtho, contrastLevels, analyzedNoiseInstance)

    lumMapRange = [0 0.8];
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', numel(contrastLevels)+1, ...
            'colsNum', 5, ...
            'heightMargin',  0.03, ...
            'widthMargin',    0.03, ...
            'leftMargin',     0.03, ...
            'rightMargin',    0.03, ...
            'bottomMargin',   0.05, ...
            'topMargin',      0.05);
   
    hFig=figure(456); clf;
    
    ax = subplot('Position', subplotPosVectors(1,1).v);
    visualizeRetinalImage(ax, nullSceneOI, lumMapRange);
    
    for theContrastLevel = 1:numel(contrastLevels)
        ax = subplot('Position', subplotPosVectors(theContrastLevel+1, 1).v);
        visualizeRetinalImage(ax, highFrequencyOIs{theContrastLevel, analyzedNoiseInstance}, lumMapRange);
        
        ax = subplot('Position', subplotPosVectors(theContrastLevel+1, 2).v);
        visualizeRetinalImage(ax, highFrequencyOIsOrtho{theContrastLevel, analyzedNoiseInstance}, lumMapRange);
        
        ax = subplot('Position', subplotPosVectors(theContrastLevel+1, 3).v);
        visualizeRetinalImage(ax, lowFrequencyOIs{theContrastLevel, analyzedNoiseInstance}, lumMapRange);
       
        ax = subplot('Position', subplotPosVectors(theContrastLevel+1, 4).v);
        visualizeRetinalImage(ax, lowFrequencyOIsOrtho{theContrastLevel, analyzedNoiseInstance}, lumMapRange);
    end
    
    drawnow;
    
end

function visualizeRetinalImage(ax, opticalImage, lumMapRange)

    spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
    % Convert spatial support in degrees
    optics = oiGet(opticalImage, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportDegs = squeeze(spatialSupportDegs(1,:,1));
    
    % retrieve the sRGB components of the optical image (just for visualization)
    XYZmap = oiGet(opticalImage, 'xyz');
    lumMap = squeeze(XYZmap(:,:,2));
    
    imagesc(ax, spatialSupportDegs, spatialSupportDegs, lumMap/max(abs(lumMap(:))).^(1/2.25));
    axis 'image'
    axis 'xy'
    ticks = max(abs(spatialSupportDegs(:)))*[-1 0 1];
    tickLabels = sprintf('%2.1f\n', ticks);
    set(gca, 'FontSize', 14, 'CLim', lumMapRange, ...
        'XLim', max(abs(spatialSupportDegs(:)))*[-1 1], 'YLim', max(abs(spatialSupportDegs(:)))*[-1 1], ...
        'XTick', ticks, 'YTick', ticks, 'XTickLabel', tickLabels, 'YTickLabel', tickLabels);
    xlabel('spatial position (degs)');
    ylabel('spatial position (degs)');
    colormap(gray(1024));
end

        