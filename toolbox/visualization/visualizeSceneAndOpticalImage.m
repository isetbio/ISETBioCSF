function visualizeSceneAndOpticalImage(modulatedScene, oiModulated, paramsList)

    sceneNormalizingY = 1/1.4;
    oiNormalizingY = 1/1.32;
    
    sceneXYZ = sceneGet(modulatedScene, 'rgb image');
    [calFormat, cols,rows] = ImageToCalFormat(sceneXYZ);
    % Normalize the luma channel
    calFormat = calFormat * sceneNormalizingY;
    sceneRGB = (xyz2rgb(calFormat'))';
    sceneRGB = CalFormatToImage(sceneRGB, cols, rows);

    
    sceneSpatialSupport = sceneGet(modulatedScene, 'spatial support');
    sceneFOV = sceneGet(modulatedScene, 'horizontalFOV');
    
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
    
    opticalImageXYZ = oiGet(oiModulated, 'xyz');
    [calFormat, cols,rows] = ImageToCalFormat(opticalImageXYZ);
    % Normalize the luma channel
    calFormat = calFormat * oiNormalizingY;
    opticalImageRGB = (xyz2rgb(calFormat'))';
    opticalImageRGB = CalFormatToImage(opticalImageRGB, cols, rows);
    
    oiFOV = oiGet(oiModulated, 'hfov');
    oiSpatialSupport = oiGet(oiModulated, 'spatial support');
    oiXSupport = squeeze(oiSpatialSupport(1,:,1));
    oiXSupport = oiXSupport / max(abs(oiXSupport(:))) * oiFOV/2;
    oiYSupport = squeeze(oiSpatialSupport(:,1,2));
    oiYSupport = oiYSupport / max(abs(oiYSupport(:))) * oiFOV/2;
    
    xx = find(abs(oiXSupport)<= max(xSupport));
    yy = find(abs(oiYSupport)<= max(ySupport));
    tmp = ones(size(opticalImageRGB));
    tmp(yy,xx,:) = opticalImageRGB(yy,xx,:);
    opticalImageRGB = tmp;
    
    
    % Draw
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 1, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.001, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.07, ...
           'topMargin',      0.001);
      
    extraMargin = 1.05;
    extraMargin = 1.00;
    
    hFig = figure(99); clf;
    formatFigureForPaper(hFig, 'figureType', 'STIMULUS_AND_OPTICAL_IMAGE');
    
    ax = subplot('Position', subplotPosVectors(1,1).v);
    image(ax, xSupport, ySupport, sceneRGB);
    set(ax, 'CLim', [0 1]);
    
    if (max(xSupport) <= 0.1)
        ticks = [-4:0.025:4];
    elseif (max(xSupport) <= 0.2)
        ticks = [-4:0.05:4];
    elseif (max(xSupport) <= 0.4)
        ticks = [-4:0.1:4];
    elseif (max(xSupport) <= 0.6)
        ticks = [-4:0.15:4];
    elseif (max(xSupport) <= 1.0)
        ticks = [-4:0.25:4];
    elseif (max(xSupport) <= 2.0)
        ticks = [-4:0.5:4];
    else
        ticks = [-4:2:4];
    end
        
        
    set(ax, 'XLim', extraMargin*max(xSupport(:))*[-1 1], 'YLim', extraMargin*max(xSupport(:))*[-1 1]);
    set(ax, 'XTick', ticks, 'YTick', ticks);
    t = 'scene';
    formatFigureForPaper(hFig, ...
            'figureType', 'STIMULUS_AND_OPTICAL_IMAGE', ...
            'theAxes', ax, ...
            'theText', t, ...
            'theTextFontSize', [], ...
            'theFigureTitle', '');
    set(ax, 'XTickLabel', {});
   
    ax = subplot('Position', subplotPosVectors(2,1).v);
    image(ax, oiXSupport, oiYSupport, opticalImageRGB);
    set(ax, 'CLim', [0 1]);

    set(ax, 'XLim', extraMargin*max(xSupport(:))*[-1 1], 'YLim', extraMargin*max(xSupport(:))*[-1 1]);
    set(ax, 'XTick', ticks, 'YTick', ticks);
    
    t = 'optical image';
    xlabel(ax, 'degrees');
    formatFigureForPaper(hFig, ...
            'figureType', 'STIMULUS_AND_OPTICAL_IMAGE', ...
            'theAxes', ax, ...
            'theText', t, ...
            'theTextFontSize', [], ...
            'theFigureTitle', '');
    
    drawnow;
    
    if (~isempty(paramsList))
        % Export to PDF
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('SceneRGBandOpticalImageRGB');
        rwObject.write(fileName, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPNG', 'FigureHandle', hFig, 'FigureType', 'png');
    end
    
end

