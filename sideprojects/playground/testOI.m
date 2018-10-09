function testOI
   
    rParams = responseParamsGenerate;
    rParams.backgroundParams.backgroundxyY = [0.31 0.31 10];
    rParams.colorModulationParams.coneContrasts = [0.5 0.5 0.5];
    rParams.spatialParams.row = 256;
    rParams.spatialParams.col = 256;
    rParams.spatialParams.fieldOfViewDegs = 0.3;
    rParams.spatialParams.cyclesPerDegree = 16;
    
    gaborScene = colorSceneCreate(rParams.spatialParams,rParams.backgroundParams,rParams.colorModulationParams,rParams.oiParams);
    gaborOI = oiCreate();
    
    padParam = 'default';
    %padParam = 'padStructLarger';
    %padParam = 'padStructSmaller';
    %padParam = 'padValue';
    %padParam = 'padSize';
    
    switch padParam
        case 'padStructLarger'
            padStruct = struct('sizeDegs', 0.5, 'value', 'zero photons');
            %padStruct = struct('value', 'zero photons');
            %padStruct = struct('sizeDegs', 1);
            %padStruct.value = 'zero photons';
            gaborOI = oiSet(gaborOI, 'pad', padStruct);
        case 'padStructSmaller'
            padStruct = struct('sizeDegs', 0.2, 'value', 'mean photons');
            %padStruct = struct('value', 'zero photons');
            %padStruct = struct('sizeDegs', 1);
            %padStruct.value = 'zero photons';
            gaborOI = oiSet(gaborOI, 'pad', padStruct);
        case 'padValue'
            padValue = 'zero photons';
            gaborOI = oiSet(gaborOI, 'padvalue', padValue);
         case 'padSize'
            padSizeDegs = 0.8;
            gaborOI = oiSet(gaborOI, 'padSizeDegs', padSizeDegs);
    end
    
    
    gaborOI = oiCompute(gaborOI, gaborScene);

    
    visDegs = 0.6;
    figNo = 1;
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 830 380], 'Color', [1 1 1], 'Name', sprintf('%s', padParam));
    visualizeLuminance(gaborScene, nan, 1, visDegs);
    visualizeLuminance(gaborOI, [], 2, visDegs);
    NicePlot.exportFigToPNG(sprintf('%s.png', padParam), hFig, 300);
end

function visualizeLuminance(obj, modulation, subplotID, visDegs)

    if (strcmp(obj.type, 'scene'))
        lumMap = sceneGet(obj, 'luminance');
        if (~isnan(modulation))
            lumMap = lumMap - sceneGet(obj, 'mean luminance');
        end
        spatialSupport = sceneGet(obj, 'spatial support');
        FOV = sceneGet(obj, 'horizontalFOV');
        plotTitle = 'luminance';
    else
        lumMap = oiGet(obj, 'illuminance');
        if (~isnan(modulation))
            lumMap = lumMap - oiGet(obj, 'mean illuminance');
        end
        FOV = oiGet(obj, 'hfov');
        spatialSupport = oiGet(obj, 'spatial support');
        lumMap(1:5,1:5)
        plotTitle = 'illuminance';
    end
    
    [xSupport, ySupport] = getXYspatialSupports(spatialSupport, FOV);
    minLum = min(lumMap(:));
    maxLum = max(lumMap(:));
    
    subplot(1,2,subplotID)
    imagesc(xSupport, ySupport, lumMap);
    hold on;
    plot([-0.8 0.8], [0 0], 'b-');
    plot([0 0], [-0.8 0.8], 'b-');
    axis 'image'
    xlabel('degs');
    
    if (isempty(modulation))
        if (strcmp(obj.type, 'scene'))
            maxVisLum = 14;
        else
            maxVisLum = 0.2;
        end
        set(gca, 'CLim', [0 maxVisLum]);
    elseif (isnan(modulation))
        set(gca, 'CLim', [minLum maxLum]);
    else
        set(gca, 'CLim', modulation*[-1 1]);
    end
    
    set(gca, 'FontSize', 14, 'XLim', 0.5*visDegs*[-1 1], 'YLim', 0.5*visDegs*[-1 1]);
    set(gca, 'XTick', [-0.8:0.1:0.8], 'YTick', [-0.8:0.1:0.8]);
    set(gca, 'XColor', 'b', 'YColor', 'b', 'LineWidth', 1.0);
    title(sprintf('stimulus %s\n(%2.3f - %2.3f)', plotTitle, minLum, maxLum));
    colormap(gray(512));
    c = colorbar;
    c.Label.String = sprintf('visualized %s range',plotTitle);
    drawnow;
end

function [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV)
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
end
