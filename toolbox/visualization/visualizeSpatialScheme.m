function visualizeSpatialScheme(spatialParams, mosaicParams, topLevelDirParams)

    % Load the mosaic
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    
    % Get the cone locs in degrees
    if (strcmp(mosaicParams.conePacking,'rect'))
        coneLocsInMeters = theMosaic.coneLocs;
    else
        coneLocsInMeters = theMosaic.coneLocsHexGrid;
    end
    coneLocsInDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width  * theMosaic.fov(1);
    coneLocsInDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    
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
     set(hFig, 'Position', [10 10 2450 1125]);
    
     switch(spatialParams.spatialType)
        case 'Gabor'  
            % Make the stimulus spatial modulation
            
            spatialPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
            spatialModulation = spatialPattern-1;
            spatialModulation = spatialModulation/max(abs(spatialModulation(:)));
            
            xaxis = (0:(size(spatialModulation,2)-1))/size(spatialModulation,2) * spatialParams.fieldOfViewDegs;
            xaxis = xaxis - mean(xaxis);
            yaxis = (0:(size(spatialModulation,1)-1))/size(spatialModulation,1) * spatialParams.fieldOfViewDegs;
            yaxis = yaxis - mean(yaxis);
             
         otherwise
     end
    
    subplot('Position', subplotPosVectors(1,1).v);
    imagesc(xaxis, yaxis, spatialModulation);
    axis 'xy';
    axis 'image';
    hold on;
    % outline mosaic extent in green
    x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
    y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
    plot(x,y, 'g-', 'LineWidth', 1.5);
    hold off
    set(gca, 'XLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'YLim', spatialParams.fieldOfViewDegs/2*[-1 1]);
    set(gca, 'FontSize', 14);
    ylabel('degrees');
    title('stimulus and cone mosaic (stimulus view)');

    subplot('Position', subplotPosVectors(1,2).v);
    imagesc(xaxis, yaxis, spatialModulation);
    axis 'xy';
    axis 'image'
    hold on;
    plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
    hold off;
    set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02);
    set(gca, 'FontSize', 14);
    xlabel('degrees');
    title('stimulus and cone mosaic (cone mosaic view)', 'FontSize', 16);

    colormap(gray(1024));
    drawnow;
end


