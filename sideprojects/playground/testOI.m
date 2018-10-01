function testOI
   
    rParams = responseParamsGenerate;
    rParams.backgroundParams.backgroundxyY = [0.31 0.31 10];
    rParams.colorModulationParams.coneContrasts = [0.5 0.5 0.5];
    rParams.spatialParams.row = 256;
    rParams.spatialParams.col = 256;
    rParams.spatialParams.fieldOfViewDegs = 0.2;
    rParams.spatialParams.cyclesPerDegree = 16;
    
    gaborScene = colorSceneCreate(rParams.spatialParams,rParams.backgroundParams,rParams.colorModulationParams,rParams.oiParams);
    gaborOI = oiCreate();
    gaborOI = oiCompute(gaborOI, gaborScene);
    
    visualizeLuminance(gaborScene, 10);
    visualizeLuminance(gaborOI, 0.1);
    
end

function visualizeLuminance(obj, modulation)

    if (strcmp(obj.type, 'scene'))
        lumMap = sceneGet(obj, 'luminance');
        lumMap = lumMap- sceneGet(obj, 'mean luminance');
        spatialSupport = sceneGet(obj, 'spatial support');
        FOV = sceneGet(obj, 'horizontalFOV');
    else
        lumMap = oiGet(obj, 'illuminance');
        lumMap = lumMap - oiGet(obj, 'mean illuminance');
        FOV = oiGet(obj, 'hfov');
        spatialSupport = oiGet(obj, 'spatial support');
    end
    
    [xSupport, ySupport] = getXYspatialSupports(spatialSupport, FOV);
    
    figure(); clf;
    imagesc(xSupport, ySupport, lumMap);
    axis 'image'
    xlabel('degs');
    set(gca, 'CLim', modulation*[-1 1])
    colormap(gray(512))
    drawnow;
end

function [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV)
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
end
