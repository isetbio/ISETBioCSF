function realizedScene = realizeSceneInDisplay(scene, display)
% Realize the desired scene into the provided display.
%
% Syntax:
%   realizedScene = realizeSceneInDisplay(scene, display)
%
% Description:
%    Transform the provided scene according to the stipulations provided by
%    the specified display.
%
% Inputs:
%    scene         - Struct. A scene structure to display.
%    display       - Struct. A structure containing information about the
%                    display to be used.
%
% Outputs:
%    realizedScene - Struct. A scene structure tailored for the display.
%
% Optional key/value pairs:
%    None.
%

    % Extract the scene's XYZ components
    sceneXYZ = sceneGet(scene, 'xyz');

    % Extract the display's RGB-to-XYZ transformation matrix
    displayRGBtoXYZ = displayGet(display, 'rgb2xyz');

    % Generate RGB primaries for generating the scene's XYZ components
    rgbPrimariesImage = ...
        imageLinearTransform(sceneXYZ, inv(displayRGBtoXYZ));

    % Extract inverse gamma table 
    inverseGammaTable = displayGet(display, 'inverse gamma');

    % Pass linear RGB primaries via an inverse gamma to generate the
    % settings values.
    rgbSettingsImage = ieLUTLinear(rgbPrimariesImage, ...
        inverseGammaTable / max(inverseGammaTable(:)));

    if (any(rgbSettingsImage(:) > 1.0))
        fprintf(2, 'Image is out of gamut > 1). Clipping to gamut.\n');
        rgbSettingsImage(rgbSettingsImage > 1.0) = 1.0;
    end

    if (any(rgbSettingsImage(:) < 0.0))
        fprintf(2, 'Image is out of gamut < 0). Clipping to gamut.\n');
        rgbSettingsImage(rgbSettingsImage < 0.0) = 0.0;
    end

    % Generate a scene based on these RGB settings
    meanLuminance = [];
    realizedScene = sceneFromFile(rgbSettingsImage, 'rgb', ...
        meanLuminance, display);

    % Set the realized scene size and view distance to match those of the
    % provided scene.
    realizedScene = sceneSet(realizedScene, ...
        'wangular', sceneGet(scene, 'wangular'));
    realizedScene = sceneSet(realizedScene, ...
        'distance', sceneGet(scene, 'distance'));

end
