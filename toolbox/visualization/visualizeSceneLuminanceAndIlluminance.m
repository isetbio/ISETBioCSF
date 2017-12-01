function visualizeSceneLuminanceAndIlluminance(modulatedScene, oiModulated, paramsList)

    % visualize the scene and the optical image
    visualizedPercentage = 0.1;
    sceneLuminance = sceneGet(modulatedScene, 'luminance');
    spatialSupport = sceneGet(modulatedScene, 'spatial support');
    xSupport = squeeze(spatialSupport(1,:,1));
    ySupport = squeeze(spatialSupport(:,1,2));
    idx = find((abs(squeeze(spatialSupport(:,:,1))) < max(xSupport)*visualizedPercentage) & ...
               (abs(squeeze(spatialSupport(:,:,2))) < max(ySupport)*visualizedPercentage));
    maxSceneLuminance = max(max(sceneLuminance(idx)));
    minSceneLuminance = min(min(sceneLuminance(idx)));      
    sceneLuminance = (sceneLuminance-minSceneLuminance)/(maxSceneLuminance-minSceneLuminance);
    [m,idx] = max(sceneLuminance(:));
    [sceneMidRow, sceneMidCol] = ind2sub(size(sceneLuminance), idx);
    
    
    oiIlluminance = oiGet(oiModulated,'illuminance');
    visualizedPercentageOI = visualizedPercentage * size(oiIlluminance,1)/size(sceneLuminance,1);
    
	oiSpatialSupport = oiGet(oiModulated, 'spatial support');
    oiXSupport = squeeze(oiSpatialSupport(1,:,1));
    oiYSupport = squeeze(oiSpatialSupport(:,1,2));
    idx = find((abs(squeeze(oiSpatialSupport(:,:,1))) < max(oiXSupport)*visualizedPercentageOI) & ...
               (abs(squeeze(oiSpatialSupport(:,:,2))) < max(oiYSupport)*visualizedPercentageOI));
    maxOIILuminance = max(max(oiIlluminance(idx)));
    minOIILuminance = min(min(oiIlluminance(idx)));
    oiIlluminance = (oiIlluminance-minOIILuminance)/(maxOIILuminance-minOIILuminance);
    [m,idx] = max(oiIlluminance(:));
    [oiMidRow, oiMidCol] = ind2sub(size(oiIlluminance), idx);
    
    % Draw
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 1, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.02);
    
       
    hFig = figure(99); clf;
    set(hFig, 'Position', [10 10 1290 2560], 'Color', [1 1 1]);
    subplot('Position', subplotPosVectors(1,1).v);
    pcolor(xSupport, ySupport, sceneLuminance);
    hold on;
    plot(xSupport+(xSupport(2)-xSupport(1))/2, squeeze(sceneLuminance(sceneMidRow,:))*ySupport(end)*visualizedPercentage*0.95, 'y-', 'LineWidth', 1.5);
    plot([0 0], [ySupport(1) ySupport(end)], 'w-');
    plot([xSupport(1) xSupport(end)], [0 0], 'w-');
    hold off
    axis 'image';
    set(gca, 'XLim', [-max(xSupport) max(xSupport)]*visualizedPercentage, 'XTick', [0], ...
             'YLim', [-max(ySupport) max(ySupport)]*visualizedPercentage, 'YTick', [0], 'CLim', [0 1]);
    set(gca, 'FontSize', 20);
    title('scene lluminance');
    
    subplot('Position', subplotPosVectors(2,1).v);
    pcolor(oiXSupport, oiYSupport, oiIlluminance);
    hold on;
    plot(oiXSupport+(oiXSupport(2)-oiXSupport(1))/2, squeeze(oiIlluminance(oiMidRow,:))*oiYSupport(end)*visualizedPercentageOI*0.95, 'y-','LineWidth', 1.5);
    plot([0 0], [ySupport(1) ySupport(end)], 'w-');
    plot([xSupport(1) xSupport(end)], [0 0], 'w-');
    hold off
    axis 'image';
    set(gca, 'XLim', [-max(oiXSupport) max(oiXSupport)]*visualizedPercentageOI, 'XTick', [0], ...
             'YLim', [-max(oiYSupport) max(oiYSupport)]*visualizedPercentageOI, 'YTick', [0], 'CLim', [0 1]);
    set(gca, 'FontSize', 20);
    title('optical image illuminance');
    colormap(jet(1024));
    drawnow;
    
    % Export to PDF
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('SceneLuminanceAndIlluminance');
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end

