function hFig = visualizeStimulusAndConeMosaic(theMosaic, thePeakOI, modulatedScene, paramsList, varargin)
    p = inputParser;
    p.addParameter('maxIllumValueToDisplay', [], @isnumeric);
    p.parse(varargin{:});
    maxIllumValue = p.Results.maxIllumValueToDisplay;
    
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
    micronsPerDegree = oiGet(thePeakOI, 'width')*1e6 / oiGet(thePeakOI, 'hfov');

    % Generate cone aperture outline
    apertureRadiusMicrons = 0.5*diameterForCircularApertureFromWidthForSquareAperture(theMosaic.pigment.pdWidth)*1e6;
    innerSegmentRadiusDegs = apertureRadiusMicrons/micronsPerDegree;
    apertureOutline.xDegs = cosd(0:15:360)*innerSegmentRadiusDegs;
    apertureOutline.yDegs = sind(0:15:360)*innerSegmentRadiusDegs;
    
    xaxis = support(1,:,1)/micronsPerDegree;
    yaxis = support(:,1,2)/micronsPerDegree;
    illumMap = oiCalculateIlluminance(thePeakOI);
    minIllum = min(illumMap(:));
    maxIllum = max(illumMap(:));

    maxIllumMap = max(illumMap(:));
    if (~isempty(maxIllumValue))
        maxIllumMap = maxIllumValue;
    end
    normIllumMap = illumMap/maxIllumMap;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 2, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.07, ...
           'topMargin',      0.04);
       
    hFig = figure(98765); clf;
    set(hFig, 'Position', [10 10 1400 780], 'Color', [1 1 1]);
            
    
    sceneLuminance = sceneGet(modulatedScene, 'luminance');
    spatialSupport = sceneGet(modulatedScene, 'spatial support');
    angularWidth = sceneGet(modulatedScene, 'wAngular');
    angularHeight = sceneGet(modulatedScene, 'hAngular');
    
    % Spatial support in degrees
    xSupport = squeeze(spatialSupport(1,:,1));
    ySupport = squeeze(spatialSupport(:,1,2));
    
    xSupport = xSupport/max(xSupport)*angularWidth/2;
    ySupport = ySupport/max(ySupport)*angularHeight/2;
    stimulusSupportRange = max([max(xSupport) max(ySupport)]);
    sceneLuminance = sceneLuminance / max(sceneLuminance(:));
    
    
    supportRange = max([stimulusSupportRange max(mosaicFOV)/2]);
    
    
    subplot('Position', subplotPosVectors(1,1).v);
    [m,idx] = max(sceneLuminance(:));
    [sceneMidRow, sceneMidCol] = ind2sub(size(sceneLuminance), idx);
    h = pcolor(xSupport, ySupport, sceneLuminance);
    set(h, 'EdgeColor', 'none');
    hold on;
    sceneLuminanceSlice = squeeze(sceneLuminance(sceneMidRow,:));
    plot(xSupport, sceneLuminanceSlice/max(sceneLuminanceSlice)*supportRange*0.95, 'r-', 'LineWidth', 1.5);
    plot([0 0], supportRange*[-1 1], 'k-');
    plot(supportRange*[-1 1], [0 0], 'k-');
    hold off
    axis 'image';
    set(gca, 'XLim', supportRange*[-1 1]*1.1, 'YLim', supportRange*[-1 1]*1.1, 'CLim', [0 1]);
    set(gca, 'FontSize', 18);
    colorbarTicks = 0:0.1:1;
    colorbarTicksLabels = colorbarTicks*max(sceneLuminanceSlice);
    h = colorbar('northoutside', 'Ticks', colorbarTicks, 'TickLabels', sprintf('%2.2f\n', colorbarTicksLabels));
    xlabel('degrees');
    title('modulated scene lluminance');
    
    
    
    subplot('Position', subplotPosVectors(1,2).v);
    %imagesc(xaxis, yaxis, illumMap);
    contourLevels = 0.0:0.05:1.0;
    contourf(xaxis, yaxis, normIllumMap, contourLevels, 'LineColor', [0.3 0.3 0.3]);
    axis 'xy';
    axis 'image'
    hold on;
    plot(xaxis, normIllumMap(floor(size(normIllumMap,1)/2),:)*supportRange*0.95, 'r-', 'LineWidth', 1.5);
    for coneIndex = 1:size(coneLocsInDegs,1)
        xAperture = coneLocsInDegs(coneIndex,1) + apertureOutline.xDegs;
        yAperture = coneLocsInDegs(coneIndex,2) + apertureOutline.yDegs;
        plot(xAperture, yAperture, '-', 'LineWidth', 1.0, 'Color', [0 0 0]);
    end
    plot(supportRange*[-1 1]*1.1, [0 0], 'k-');
    plot([0 0], supportRange*[-1 1]*1.1, 'k-');
    
    hold off;
    set(gca, 'CLim', [0 1], 'YTickLabel', {}, 'XLim', supportRange*[-1 1]*1.1, 'YLim', supportRange*[-1 1]*1.1);
    set(gca, 'FontSize', 18);
    colorbarTicks = contourLevels(1:2:end);
    colorbarTicksLabels = colorbarTicks*maxIllumMap;
    h = colorbar('northoutside','Ticks', colorbarTicks, 'TickLabels', sprintf('%2.2f\n', colorbarTicksLabels));
    
    xlabel('degrees');
    title(sprintf('retinal stimulus (illuminance: %f - %f)', minIllum, maxIllum), 'FontSize', 16);
    

    cMap = brewermap(1024,'*YlGnBu');
    colormap(cMap);
    drawnow;
    
    if (~isempty(paramsList))
        % Save figure
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('SpatialScheme');
        rwObject.write(fileName, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPNG', 'FigureHandle', hFig, 'FigureType', 'png');
    end
    
end


