function visualizeSceneOpticalImageAndMeanResponses(modulatedScene, oiModulated, theMosaic, noiseFreeIsomerizations, noiseFreePhotocurrents, paramsList)

    sceneXYZ = sceneGet(modulatedScene, 'xyz');
    oiXYZ = oiGet(oiModulated, 'xyz');
    [sceneSRGB, ~, ~] = xyz2srgb(sceneXYZ);
    [oiSRGB, ~, ~] = xyz2srgb(oiXYZ);
    
    sceneSpatialSupport = sceneGet(modulatedScene, 'spatial support');
    sceneFOV = sceneGet(modulatedScene, 'horizontalFOV');
    [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV);
    
    oiFOV = oiGet(oiModulated, 'hfov');
    oiSpatialSupport = oiGet(oiModulated, 'spatial support');
    [oiXSupport, oiYSupport] = getXYspatialSupports(oiSpatialSupport, oiFOV);
    
    
    
    noiseFreeIsomerizations = squeeze(noiseFreeIsomerizations(1,:,:)) / theMosaic.integrationTime;
    isomerizationsRange = [min(noiseFreeIsomerizations(:)) max(noiseFreeIsomerizations(:))];
    isomerizationsRange = [2800 3600];
    isomerizationsTickIncrement = 200;
    photocurrentsRange = [min(noiseFreePhotocurrents(:)) max(noiseFreePhotocurrents(:))];
    photocurrentsRange = [-86 -74];
    photocurrentsTickIncrement = 2;
    % Draw
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 4, ...
           'heightMargin',   0.001, ...
           'widthMargin',    0.025, ...
           'leftMargin',     0.06, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.01, ...
           'topMargin',      0.03);
      
    extraMargin = 1.05;
    extraMargin = 1.00;
    
    hFig = figure(99); clf;
    formatFigureForPaper(hFig, 'figureType', 'STIMULUS_OPTICAL_IMAGE_ISOMERIZATIONS_PHOTOCURRENTS');
    
    ax = subplot('Position', subplotPosVectors(1,1).v);
    image(ax, xSupport, ySupport, sceneSRGB);
    
    if (max(xSupport) <= 0.1)
        ticks = -4:0.025:4;
    elseif (max(xSupport) <= 0.2)
        ticks = -4:0.05:4;
    elseif (max(xSupport) <= 0.4)
        ticks = -4:0.1:4;
    elseif (max(xSupport) <= 0.6)
        ticks = -4:0.15:4;
    elseif (max(xSupport) <= 1.0)
        ticks = -4:0.25:4;
    elseif (max(xSupport) <= 2.0)
        ticks = -4:0.5:4;
    else
        ticks = -4:2:4;
    end
        
    XLims = extraMargin*max(xSupport(:))*[-1 1];
    YLims = extraMargin*max(ySupport(:))*[-1 1];
    
    xlabel(ax, 'degrees');
    ylabel(ax, 'degrees');
    formatFigureForPaper(hFig, ...
            'figureType', 'STIMULUS_OPTICAL_IMAGE_ISOMERIZATIONS_PHOTOCURRENTS', ...
            'theAxes', ax, ...
            'theTextFontSize', [], ...
            'theFigureTitle', 'scene');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', ticks, 'YTick', ticks);
    %set(ax, 'XTickLabel', {});
    set(ax, 'CLim', [0 1]);
    
    
    
    ax = subplot('Position', subplotPosVectors(1,2).v);
    image(ax, oiXSupport, oiYSupport, oiSRGB);
    formatFigureForPaper(hFig, ...
            'figureType', 'STIMULUS_OPTICAL_IMAGE_ISOMERIZATIONS_PHOTOCURRENTS', ...
            'theAxes', ax, ...
            'theTextFontSize', [], ...
            'theFigureTitle', 'retinal image');
    set(ax, 'XLim', XLims, 'YLim', YLims, 'XTick', ticks, 'YTick', ticks);
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    set(ax, 'CLim', [0 1]);
    

    %activationColorMap = brewermap(1024, '*Spectral');
    activationColorMap = brewermap(1024, '*Greys');
    visualizedConeAperture = 'geometricArea';
    mapType = 'modulated disks';
    
    ax = subplot('Position', subplotPosVectors(1,3).v);
    [~,idx] = max(noiseFreeIsomerizations(:));
    [peakCone, peakTimeBin] = ind2sub(size(noiseFreeIsomerizations), idx);
    t = theMosaic.timeAxis;
    fprintf('Time bin of max isomerization response: %d (%2.1f msec)\n', peakTimeBin, t(peakTimeBin)*1000);
    theMosaic.renderActivationMap(ax, noiseFreeIsomerizations(:, peakTimeBin),  ...
        'signalRange', [], ...
        'visualizedConeAperture', visualizedConeAperture, ...
        'mapType', mapType , ...
        'colorMap', activationColorMap);
    
    formatFigureForPaper(hFig, ...
            'figureType', 'STIMULUS_OPTICAL_IMAGE_ISOMERIZATIONS_PHOTOCURRENTS', ...
            'theAxes', ax, ...
            'theText', t, ...
            'theTextFontSize', [], ...
            'theFigureTitle', 'cone isomerizations');
    set(ax, 'CLim', [0 1], ...
            'XLim', XLims*theMosaic.micronsPerDegree*1e-6, ...
            'YLim', YLims*theMosaic.micronsPerDegree*1e-6, ...
            'XTick', ticks*theMosaic.micronsPerDegree*1e-6, ...
            'YTick', ticks*theMosaic.micronsPerDegree*1e-6, ...
            'XTickLabel', sprintf('%2.2f\n', ticks), ...
            'YTickLabel', sprintf('%2.2f\n', ticks) ...
            );
        
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    set(ax, 'Color', [0 0 0]);
    colorbarLabel = sprintf('R*/cone/sec');
    colorbarTickLabels = isomerizationsRange(1):isomerizationsTickIncrement:isomerizationsRange(2);
    colorbarTicks = (colorbarTickLabels - isomerizationsRange(1))/( isomerizationsRange(2)- isomerizationsRange(1));
    renderColorBar(ax, colorbarTicks, colorbarTickLabels, colorbarLabel);
    

    ax = subplot('Position', subplotPosVectors(1,4).v);
    noiseFreePhotocurrents = squeeze(noiseFreePhotocurrents(1,:,:));
    [~,idx] = max(noiseFreePhotocurrents(:));
    [peakCone, peakTimeBin] = ind2sub(size(noiseFreePhotocurrents), idx);
    t = theMosaic.timeAxis;
    fprintf('Time bin of max photocurrent response: %d (%2.1f msec)\n', peakTimeBin, t(peakTimeBin)*1000);
    theMosaic.renderActivationMap(ax, noiseFreePhotocurrents(:, peakTimeBin),  ...
        'signalRange', photocurrentsRange, ...
        'visualizedConeAperture', visualizedConeAperture, ...
        'mapType', mapType, ...
        'colorMap', activationColorMap);
    formatFigureForPaper(hFig, ...
            'figureType', 'STIMULUS_OPTICAL_IMAGE_ISOMERIZATIONS_PHOTOCURRENTS', ...
            'theAxes', ax, ...
            'theTextFontSize', [], ...
            'theFigureTitle', 'cone photocurrents');
    set(ax, 'CLim', [0 1], ...
            'XLim', XLims*theMosaic.micronsPerDegree*1e-6, ...
            'YLim', YLims*theMosaic.micronsPerDegree*1e-6, ...
            'XTick', ticks*theMosaic.micronsPerDegree*1e-6, ...
            'YTick', ticks*theMosaic.micronsPerDegree*1e-6, ...
            'XTickLabel', sprintf('%2.2f\n', ticks), ...
            'YTickLabel', sprintf('%2.2f\n', ticks) ...
            );
    set(ax, 'XTickLabel', {}); set(ax, 'YTickLabel', {});
    set(ax, 'Color', [0 0 0]);
    
    colorbarLabel = sprintf('pAmps');
    colorbarTickLabels = photocurrentsRange(1):photocurrentsTickIncrement:photocurrentsRange(2);
    colorbarTicks = (colorbarTickLabels - photocurrentsRange(1))/( photocurrentsRange(2)- photocurrentsRange(1));
    renderColorBar(ax, colorbarTicks, colorbarTickLabels, colorbarLabel);
        
    
    drawnow;
    
    if (~isempty(paramsList))
        % Export to PDF
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('SceneOpticalImageIsomerizationsPhotocurrents');
        rwObject.write(fileName, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
    end
    
end

function renderColorBar(ax, colorbarTicks, colorbarTickLabels,  colorbarLabel)
    
    originalPosition = get(ax, 'position');
    hCbar = colorbar('Ticks', colorbarTicks, 'TickLabels', sprintf('%2.0f\n',colorbarTickLabels));
    hCbar.Orientation = 'horizontal'; 
    hCbar.Location = 'SouthOutside';
    %hCbar.AxisLocation = 'in';
    hCbar.Label.String = colorbarLabel;
    hCbar.FontSize = 16; 
    %hCbar.FontName = 'Menlo'; 
    %hCbar.FontWeight = 'Bold'; 
    hCbar.Color = [0 0 0];
    % The addition changes the figure size, so undo this change
    newPosition = get(gca, 'position');
    set(ax,'position',[newPosition(1) originalPosition(2) originalPosition(3) originalPosition(4)]);
    set(ax,'position',[newPosition(1) originalPosition(2) originalPosition(3) originalPosition(4)]);
    
end

function [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV)
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
end
