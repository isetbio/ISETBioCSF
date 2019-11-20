function visualizeLargeMosaic()
    load('HexConeMosaic15DegsBest.mat', 'theMosaic');
    geomStruct = theMosaic.geometryStruct();
    
    % Visualize only part of the mosaic X(degs): [0 7], Y(degs): [-0.25 0.25]
    visualizedFOV = struct('xo', 3.5, 'yo', 0.0, 'width', 7., 'height', 0.4);
    
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    eccDegs = 0.02;
    skipCells = 5;
    
    while (max(eccDegs) < 7)
        nextPosDegs = skipCells*WatsonRGCCalc.midgetRGCRFSpacing(eccDegs(end), 'temporal meridian', 'single polarity', 'deg');
        eccDegs = cat(2,eccDegs, eccDegs(end)+nextPosDegs);
    end
    

    midgetRGCRFDensity = WatsonRGCCalc.midgetRGCRFDensity(eccDegs, 'temporal meridian', 'RFs per deg2');
    [coneSpacingDegs, coneDensity] = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, 'temporal meridian', 'Cones per deg2');
    midgetRGCRFtoConeRatio = (0.5*midgetRGCRFDensity) ./ coneDensity;
    conesPerMidgetRGC = 1./midgetRGCRFtoConeRatio;

    for eccIndex = 1:numel(eccDegs)
        ganglionPos = [eccDegs(eccIndex) 0];
        [~,coneIndex] = min(sum((bsxfun(@minus, geomStruct.coneLocs, ganglionPos)).^2,2));
        radius = conesPerMidgetRGC(eccIndex) * 0.5*geomStruct.coneApertures(coneIndex);
        RFcenter(eccIndex,1,:) = geomStruct.coneLocs(coneIndex,1) + radius*cosd(0:5:360);
        RFcenter(eccIndex,2,:) = geomStruct.coneLocs(coneIndex,2) + radius*sind(0:5:360);
        RFsurround(eccIndex,1,:) = geomStruct.coneLocs(coneIndex,1) + 5*radius*cosd(0:5:360);
        RFsurround(eccIndex,2,:) = geomStruct.coneLocs(coneIndex,2) + 5*radius*sind(0:5:360);
    end

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2640 1525], 'Color', [0.2 0.2 0.2]);
    axesHandle  = subplot('Position', [0.025 0.04+0.64 0.97 0.30]);
    theMosaic.visualizeGrid('axesHandle', axesHandle, 'ticksInVisualDegs', true, ...
        'visualizedFOV', visualizedFOV, 'backgroundColor', [0.2 0.2 0.2]);
    hold(axesHandle, 'on');
    for eccIndex = 1:numel(eccDegs)
        plotCenterSurroundRF(axesHandle, ...
            squeeze(RFcenter(eccIndex,:,:))*theMosaic.micronsPerDegree*1e-6, ...
            squeeze(RFsurround(eccIndex,:,:))*theMosaic.micronsPerDegree*1e-6 ...
            );
    end
    set(axesHandle, 'XLim', [0 2.33]*theMosaic.micronsPerDegree*1e-6, 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    xlabel(axesHandle,'');
    
    axesHandle  = subplot('Position', [0.025 0.04+0.32 0.97 0.30]);
    theMosaic.visualizeGrid('axesHandle', axesHandle, 'ticksInVisualDegs', true, ...
        'visualizedFOV', visualizedFOV, 'backgroundColor', [0.2 0.2 0.2]);
    hold(axesHandle, 'on');
    for eccIndex = 1:numel(eccDegs)
        plotCenterSurroundRF(axesHandle, ...
            squeeze(RFcenter(eccIndex,:,:))*theMosaic.micronsPerDegree*1e-6, ...
            squeeze(RFsurround(eccIndex,:,:))*theMosaic.micronsPerDegree*1e-6 ...
            );
    end
    set(axesHandle, 'XLim', [2.33 4.66]*theMosaic.micronsPerDegree*1e-6, 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    xlabel(axesHandle,'');
    
    axesHandle  = subplot('Position', [0.025 0.04 0.97 0.30]);
    theMosaic.visualizeGrid('axesHandle', axesHandle, 'ticksInVisualDegs', true, ...
        'visualizedFOV', visualizedFOV, 'backgroundColor', [0.2 0.2 0.2]);
    hold(axesHandle, 'on');
    for eccIndex = 1:numel(eccDegs)
        plotCenterSurroundRF(axesHandle, ...
            squeeze(RFcenter(eccIndex,:,:))*theMosaic.micronsPerDegree*1e-6, ...
            squeeze(RFsurround(eccIndex,:,:))*theMosaic.micronsPerDegree*1e-6 ...
            );
    end
    set(axesHandle, 'XLim', [4.66 7.0]*theMosaic.micronsPerDegree*1e-6, 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9]);
    
    NicePlot.exportFigToPDF('midgetRGC', hFig, 300);
end

function plotCenterSurroundRF(axesHandle, center, surround)
    f = 1:size(center,2);
    centerColor = [0.8 0.8 0.8];
    surroundColor = [0.9 0.9 0.9];
    
    % Plot the surround
    x = [surround(1,:) nan center(1,:)];
    y = [surround(2,:) nan center(2,:)];
    pgon = polyshape(x,y);
    plot(pgon, 'FaceColor', surroundColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 1.0, 'EdgeColor', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
    
    patch(axesHandle, 'Faces', f, 'Vertices', center', ...
        'EdgeColor', 'k', 'FaceColor', centerColor, 'FaceAlpha', 0.2, 'EdgeAlpha', 1.0, 'LineWidth', 1.5);

end
