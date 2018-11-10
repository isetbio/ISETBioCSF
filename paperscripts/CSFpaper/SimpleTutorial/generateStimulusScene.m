function scene = generateStimulusScene(stimParams, presentationDisplay)
    %% Compute the display's pixel size in degrees to figure out the scene pixels
    displayPixelSizeMeters = displayGet(presentationDisplay, 'sample spacing');
    viewingDistanceMeters = displayGet(presentationDisplay, 'viewing distance');
    displayPixelSizeDegrees = 2 * atand(0.5*displayPixelSizeMeters/viewingDistanceMeters);
    scenePixelsNum = round(stimParams.widthDegs/displayPixelSizeDegrees(1));
    
    imageHarmonicParams = struct(...
        'freq', stimParams.widthDegs * stimParams.spatialFrequencyCyclesPerDeg, ...
        'ang', stimParams.orientationDegs/180*pi, ...
        'contrast', stimParams.contrast, ...
        'GaborFlag', 0.2, ...
        'row', scenePixelsNum , 'col', scenePixelsNum);
    
    %% Step 3. Generate a scene representing a 10 c/deg Gabor stimulus
    scene = sceneCreate('harmonic', imageHarmonicParams);
    scene = sceneSet(scene, 'distance', viewingDistanceMeters);
    scene = sceneSet(scene, 'wangular', stimParams.widthDegs);
    
    % Adjust scene luminance to desired level
    scene = sceneAdjustLuminance(scene, stimParams.meanLuminanceCdPerM2);
end