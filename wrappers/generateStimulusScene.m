function scene = generateStimulusScene(stimParams, presentationDisplay)
% Generate a stimulus scene
%
% Syntax:
%   scene = generateStimulusScene(stimParams, presentationDisplay)
%
% Description:
%    Generate a stimulus scene for the provided display, using supplied
%    stimulus parameters.
%
% Inputs:
%    stimParams          - Struct. A structure of stimulus parameters,
%                          containing numeric values for spatial frequency
%                          cycles per degree, orientation degrees, width
%                          degrees, contrast, and mean luminance cd per m^2
%    presentationDisplay - Struct. A display structure.
%
% Outputs:
%    scene               - Struct. A scene structure.
%
% Optional key/value pairs:
%    None.
%

    %% Compute the display's pixel size
    % Do this in degrees to figure out the scene pixels.
    displayPixelSizeMeters = ...
        displayGet(presentationDisplay, 'sample spacing');
    viewingDistanceMeters = ...
        displayGet(presentationDisplay, 'viewing distance');
    displayPixelSizeDegrees = 2 * ...
        atand(0.5 * displayPixelSizeMeters / viewingDistanceMeters);
    scenePixelsNum = ...
        round(stimParams.widthDegs / displayPixelSizeDegrees(1));

    imageHarmonicParams = struct('freq', ...
        stimParams.widthDegs * stimParams.spatialFrequencyCyclesPerDeg, ...
        'ang', stimParams.orientationDegs / 180 * pi, ...
        'contrast', stimParams.contrast, ...
        'GaborFlag', 0.2, 'row', scenePixelsNum , 'col', scenePixelsNum);

    %% Step 3. Generate a scene representing a 10 c/deg Gabor stimulus
    scene = sceneCreate('harmonic', imageHarmonicParams);
    scene = sceneSet(scene, 'distance', viewingDistanceMeters);
    scene = sceneSet(scene, 'wangular', stimParams.widthDegs);

    % Adjust scene luminance to desired level
    scene = sceneAdjustLuminance(scene, stimParams.meanLuminanceCdPerM2);

end