function realizedScene = realizeSceneInDisplay(scene, display)
        
    % Extract the scene's XYZ components
    sceneXYZ = sceneGet(scene, 'xyz');
    
    % Extract the display's RGB-to-XYZ transformation matrix
    displayRGBtoXYZ = displayGet(display, 'rgb2xyz');
    
    % Generate RGB primaries for generating the scene's XYZ components
    rgbPrimariesImage = imageLinearTransform(sceneXYZ, inv(displayRGBtoXYZ));

    % Extract inverse gamma table 
    inverseGammaTable = displayGet(display, 'inverse gamma');
    
    % Pass linear RGB primaries via inverse gamma to generate settings values
    rgbSettingsImage = ieLUTLinear(rgbPrimariesImage, inverseGammaTable/max(inverseGammaTable(:)));
    
    if (any(rgbSettingsImage(:)>1.0))
        fprintf(2,'Image is out of gamut > 1). Clipping to gamut.\n');
        rgbSettingsImage(rgbSettingsImage>1.0) = 1.0;
    end
    
    if (any(rgbSettingsImage(:)<0.0))
        fprintf(2,'Image is out of gamut < 0). Clipping to gamut.\n');
        rgbSettingsImage(rgbSettingsImage<0.0) = 0.0;
    end
 
    % Generate a scene based on these RGB settings
    meanLuminance = [];
    realizedScene = sceneFromFile(rgbSettingsImage, 'rgb', meanLuminance, display);

    % Set the realized scene size and view distance to match those of the scene
    realizedScene = sceneSet(realizedScene, 'wangular', sceneGet(scene,'wangular'));
    realizedScene = sceneSet(realizedScene, 'distance', sceneGet(scene,'distance'));
end
