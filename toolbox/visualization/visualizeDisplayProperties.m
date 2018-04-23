function hFig = visualizeDisplayProperties(theDisplay, backgroundParams, stimulusScene)

    hFig = figure(); clf; 
    formatFigureForPaper(hFig, 'figureType', 'DISPLAY_PROPERTIES');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',   0.1, ...
           'widthMargin',    0.065, ...
           'leftMargin',     0.06, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.04, ...
           'topMargin',      -0.01);
  
    subplot('Position', subplotPosVectors(1,1).v);
    axesHandle = generateSPDplot(theDisplay, backgroundParams);
    formatFigureForPaper(hFig, ...
        'figureType', 'DISPLAY_PROPERTIES', ...
        'theAxes', axesHandle);
    xtickformat(axesHandle, '%3.0f'); 
    
    subplot('Position', subplotPosVectors(1,2).v); 
    axesHandle = generateCIEplot(theDisplay, backgroundParams);
    formatFigureForPaper(hFig, ...
        'figureType', 'DISPLAY_PROPERTIES', ...
        'theAxes', axesHandle);
    
    
    subplot('Position', subplotPosVectors(1,3).v);
    axesHandle = generateDisplayedStimulus(stimulusScene, backgroundParams);
    formatFigureForPaper(hFig, ...
        'figureType', 'DISPLAY_PROPERTIES', ...
        'theAxes', axesHandle);
    
    subplot('Position', subplotPosVectors(2,1).v);
    axesHandle = generateRadianceMap(stimulusScene, 470, true);
    formatFigureForPaper(hFig, ...
        'figureType', 'DISPLAY_PROPERTIES', ...
        'theAxes', axesHandle);
    
    subplot('Position', subplotPosVectors(2,2).v);
    axesHandle = generateRadianceMap(stimulusScene, 550, false);
    formatFigureForPaper(hFig, ...
        'figureType', 'DISPLAY_PROPERTIES', ...
        'theAxes', axesHandle);
    
    subplot('Position', subplotPosVectors(2,3).v);
    axesHandle = generateRadianceMap(stimulusScene, 630, false);
    formatFigureForPaper(hFig, ...
        'figureType', 'DISPLAY_PROPERTIES', ...
        'theAxes', axesHandle);

end



function axesHandle = generateRadianceMap(stimulusScene, targetWavelength, displayYlabel)
    sceneSpatialSupport = sceneGet(stimulusScene, 'spatial support');
    sceneSpectralSupport = sceneGet(stimulusScene, 'wave');
    sceneFOV = sceneGet(stimulusScene, 'horizontalFOV');
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
    
    scenePhotons = sceneGet(stimulusScene, 'photons');
    maxPhotons = max(scenePhotons(:));
    [~,idx] = min(abs(sceneSpectralSupport-targetWavelength));
    
    imagesc(xSupport, ySupport, squeeze(scenePhotons(:,:,idx)));
    set(gca, 'CLim', [0 maxPhotons],  'XTick', -0.5:0.1:0.5, 'YTick', -0.5:0.1:0.5);
    axis 'square'
    colorbar();
    colorbarLabel = sprintf('photons/sample/sec @ %2.0fnm', sceneSpectralSupport(idx));
    colorbarTicks = (1:2:10)*1.0000e+15;
    colorbarTickLabels = sprintf('%2.0fe15\n', colorbarTicks/1e+15);
    renderColorBar(gca, colorbarTicks, colorbarTickLabels, colorbarLabel);
    colormap(gray(1024));
    
    grid 'on'; box 'on';
    xlabel('x (degs)', 'FontWeight', 'bold');
    if (displayYlabel)
        ylabel('x (degs)', 'FontWeight', 'bold');
    end
    axesHandle = gca;
    
end


function renderColorBar(ax, colorbarTicks, colorbarTickLabels,  colorbarLabel)
    originalPosition = get(ax, 'position');
    hCbar = colorbar('Ticks', colorbarTicks, 'TickLabels', colorbarTickLabels);
    hCbar.Orientation = 'horizontal'; 
    hCbar.Location = 'NorthOutside';
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


function axesHandle = generateDisplayedStimulus(stimulusScene, backgroundParams)

    sceneFOV = sceneGet(stimulusScene, 'horizontalFOV');
    sceneSpatialSupport = sceneGet(stimulusScene, 'spatial support');
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
    
    sceneRGB = sceneGet(stimulusScene, 'rgb image');
    image(xSupport, ySupport, sceneRGB);
    set(gca, 'CLim', [0 1], 'XTick', [-0.5:0.1:0.5], 'YTick', [-0.5:0.1:0.5]);
    axis 'square'
    grid 'on'; box 'on';
    xlabel('x (degs)', 'FontWeight', 'bold');
    ylabel('x (degs)', 'FontWeight', 'bold');
    axesHandle = gca;
end

function axesHandle = generateGammaFunctionPlot(theDisplay)
    gammaTable = displayGet(theDisplay, 'gamma table');
    dacsize = displayGet(theDisplay, 'dacsize')
    vd = displayGet(theDisplay,'viewingdistance')
    dpi = displayGet(theDisplay, 'dpi');
    gammaIn = 0:(size(gammaTable,1)-1);
    gammaInTicks = 0:max(gammaIn)/5:max(gammaIn);
    gammaInTickLabels = sprintf('%0.1f\n', gammaInTicks/2^dacsize);
    [xx,yy] = makeStepSeries(gammaIn, squeeze(gammaTable(:,1)));
    plot(xx(1:end-1),yy(1:end-1), 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', [0 max(gammaIn)], 'YLim', [0 1], 'YTick', 0:0.2:1.0, ...
        'XTick', gammaInTicks, 'XTickLabel', gammaInTickLabels);
    axis 'square'
    grid 'on'; box 'on';
    xlabel('settings value', 'FontWeight', 'bold');
    ylabel('primary value', 'FontWeight', 'bold');
    axesHandle = gca;
end


function axesHandle = generateSPDplot(theDisplay, backgroundParams)
    displayChannelWavelengths = displayGet(theDisplay,'wave');
    displayChannelSpectra = displayGet(theDisplay,'spd') / backgroundParams.lumFactor;
    x = displayChannelWavelengths';
    y = displayChannelSpectra(:,1)*1e3;
    [xx,yy] = makeStepSeries(x,y);
    fill(xx,yy, [1 0.7 0.7], 'EdgeColor', 'r', 'LineWidth', 1.5, 'FaceAlpha', 0.5);
    hold on;
    y = displayChannelSpectra(:,2)*1e3;
    [xx,yy] = makeStepSeries(x,y);
    fill(xx,yy, [0.7 1.0 0.7], 'EdgeColor', [0 0.8 0], 'LineWidth', 1.5, 'FaceAlpha', 0.5);
    y = displayChannelSpectra(:,3)*1e3;
    [xx,yy] = makeStepSeries(x,y);
    fill(xx,yy, [0.7 0.7 1], 'EdgeColor', 'b', 'LineWidth', 1.5, 'FaceAlpha', 0.5);
    set(gca, 'XLim', [displayChannelWavelengths(1) displayChannelWavelengths(end)], ...
        'XTick', [300:50:800], 'YTick', [0:1:5], 'YTickLabel', sprintf('%2.1f\n', [0:1:5]));
    axis 'square'
    grid 'on'; box 'on';
    xlabel('wavelength (nm)', 'FontWeight', 'bold');
    ylabel('energy (mWatts/nm)', 'FontWeight', 'bold');
    axesHandle = gca;
end

function [xx,yy] = makeStepSeries(x,y)
    xx = x(1);
    yy = 0;
    for kk = 1:numel(x)-1
        xx(numel(xx)+1) = x(kk);
        yy(numel(yy)+1) = y(kk);
        xx(numel(xx)+1) = x(kk+1);
        yy(numel(yy)+1) = y(kk);
    end
    xx(numel(xx)+1) = x(kk+1);
    yy(numel(yy)+1) = y(kk+1);
    xx(numel(xx)+1) = x(kk+1);
    yy(numel(yy)+1) = 0;
end
    
function axesHandle = generateCIEplot(theDisplay, backgroundParams)
    load('T_xyzCIEPhys2')
    T_XYZ = 683*T_xyzCIEPhys2;
    S_XYZ = S_xyzCIEPhys2;
 
    displayChannelWavelengths = displayGet(theDisplay,'wave');
    displayChannelS = WlsToS(displayChannelWavelengths);
    displayChannelSpectra = displayGet(theDisplay,'spd') / backgroundParams.lumFactor;
    
    T_XYZforDisplay = SplineCmf(S_XYZ,T_XYZ,displayChannelWavelengths);
    primaryXYZ = displayChannelS(2)* displayChannelSpectra' * T_XYZforDisplay';
    vLambda = squeeze(T_XYZforDisplay(2,:));
    primaryLuminances = displayChannelS(2)*sum(bsxfun(@times, displayChannelSpectra, vLambda'),1)
    
    s = sum(primaryXYZ,2);
    displayXChroma = primaryXYZ(:,1) ./ s;
    displayYChroma = primaryXYZ(:,2) ./ s;
    displayXChroma(end+1) = displayXChroma(1);
    displayYChroma(end+1) = displayYChroma(1);
    
    %chromaticity coordinates of spectrum locus
    wave = SToWls(S_XYZ);
    idx = find(wave > 750);
    s = sum(T_xyzCIEPhys2,1);
    s = s(1:idx(1));
    
    locusX = T_xyzCIEPhys2(1,1:numel(s)) ./ s;
    locusY = T_xyzCIEPhys2(2,1:numel(s)) ./ s;
    hold on
    cieColorPatch(locusX, locusY);
    plot(displayXChroma, displayYChroma, 'w-', 'LineWidth', 1.5);
    plot(displayXChroma(1,:), displayYChroma(1,:), 'wo-',  'MarkerFaceColor', [1 0 0], ...
        'MarkerSize', 12, 'LineWidth', 1.5);
    plot(displayXChroma(2,:), displayYChroma(2,:), 'wo-',  'MarkerFaceColor', [0 1 0], ...
        'MarkerSize', 12, 'LineWidth', 1.5);
    plot(displayXChroma(3,:), displayYChroma(3,:), 'wo-',  'MarkerFaceColor', [0 0 1], ...
        'MarkerSize', 12, 'LineWidth', 1.5);
    
    plot(backgroundParams.backgroundxyY(1), backgroundParams.backgroundxyY(2), 'ks', ...
        'MarkerSize', 14, 'LineWidth', 1.5);
    set(gca, 'XLim', [0 0.85], 'YLim', [0 0.85], ...
        'XTick', 0:0.1:1,  'YTick', 0:0.1:1);
    axis 'square'
    grid 'on'; box 'on';
    xlabel('x-chroma', 'FontWeight', 'bold');
    ylabel('y-chroma', 'FontWeight', 'bold');
    axesHandle = gca;
end


function cieColorPatch(x,y)
    N = length(x);
    i = 1;
    e = 1/3;
    steps = 25;
    xy4rgb = zeros(N*steps*4,5,'double');
    for w = 1:N                                     % wavelength
        w2 = mod(w,N)+1;
        a1 = atan2(y(w)  -e,x(w)  -e);              % start angle
        a2 = atan2(y(w2) -e,x(w2) -e);              % end angle
        r1 = ((x(w)  - e)^2 + (y(w)  - e)^2)^0.5;   % start radius
        r2 = ((x(w2) - e)^2 + (y(w2) - e)^2)^0.5;   % end radius
        for c = 1:steps                               % colourfulness
            % patch polygon
            xyz(1,1) = e+r1*cos(a1)*c/steps;
            xyz(1,2) = e+r1*sin(a1)*c/steps;
            xyz(1,3) = 1 - xyz(1,1) - xyz(1,2);
            xyz(2,1) = e+r1*cos(a1)*(c-1)/steps;
            xyz(2,2) = e+r1*sin(a1)*(c-1)/steps;
            xyz(2,3) = 1 - xyz(2,1) - xyz(2,2);
            xyz(3,1) = e+r2*cos(a2)*(c-1)/steps;
            xyz(3,2) = e+r2*sin(a2)*(c-1)/steps;
            xyz(3,3) = 1 - xyz(3,1) - xyz(3,2);
            xyz(4,1) = e+r2*cos(a2)*c/steps;
            xyz(4,2) = e+r2*sin(a2)*c/steps;
            xyz(4,3) = 1 - xyz(4,1) - xyz(4,2);
            
            % compute sRGB for vertices
            rgb = xyz2srgb(reshape(xyz, [2 2 3]));
            rgb = reshape(rgb, [4 3]);
            % store the results
            xy4rgb(i:i+3,1:2) = xyz(:,1:2);
            xy4rgb(i:i+3,3:5) = rgb;
            i = i + 4;
        end
    end
    [rows cols] = size(xy4rgb);
    f = [1 2 3 4];
    v = zeros(4,3,'double');
    saturation = 1;
    exponent = 1;
    for i = 1:4:rows
        v(:,1:2) = xy4rgb(i:i+3,1:2);
        patch('Vertices',v, 'Faces',f, 'EdgeColor','none', ...
        'FaceVertexCData',((xy4rgb(i:i+3,3:5))*saturation+(1-saturation)*ones(4,3)).^exponent,'FaceColor','interp')
    end
    
    plot(x, y, 'k-');
end



