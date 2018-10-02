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
    
    if (1==1)
    padStruct = struct('sizeDegs', 1, 'value', 'mean photons');
    padStruct = struct('value', 'zero photons');
    %padStruct = struct('sizeDegs', 1);
    %padStruct.value = 'zero photons';
    gaborOI = oiSet(gaborOI, 'pad', padStruct);
    end
    
    gaborOI = oiCompute(gaborOI, gaborScene);
    
    
    figNo = 1; 
    visualizeLuminance(gaborScene, 10); figNo = figNo + 1;
    visualizeLuminance(gaborOI, 0.1); figNo = figNo + 1;
    
    figNo = 3; 
    visualizeLuminance(gaborOI, nan); figNo = figNo + 1;
    
end

function visualizeLuminance(obj, modulation)

    if (strcmp(obj.type, 'scene'))
        lumMap = sceneGet(obj, 'luminance');
        if (~isnan(modulation))
            lumMap = lumMap - sceneGet(obj, 'mean luminance');
        end
        spatialSupport = sceneGet(obj, 'spatial support');
        FOV = sceneGet(obj, 'horizontalFOV');
    else
        lumMap = oiGet(obj, 'illuminance');
        if (~isnan(modulation))
            lumMap = lumMap - oiGet(obj, 'mean illuminance');
        end
        FOV = oiGet(obj, 'hfov');
        spatialSupport = oiGet(obj, 'spatial support');
        lumMap(1:5,1:5)
    end
    
    [xSupport, ySupport] = getXYspatialSupports(spatialSupport, FOV);
    
    figure(); clf;
    imagesc(xSupport, ySupport, lumMap);
    axis 'image'
    xlabel('degs');
    if (~isnan(modulation))
        set(gca, 'CLim', modulation*[-1 1]);
    else
        set(gca, 'CLim', [min(lumMap(:)) max(lumMap(:))]);
    end
    colormap(gray(512));
    colorbar;
    drawnow;
end

function [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV)
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
end
