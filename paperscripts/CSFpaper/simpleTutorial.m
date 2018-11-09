function simpleTutorial
    %% Set the viewing distance
    viewingDistanceMeters = 0.57;
    
    %% Step 1. Generate a display for presenting stimuli and place it at the desired viewing distance
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', viewingDistanceMeters);
    
    %% Step 2. Specify stimulus params
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', 10, ... 
        'orientationDegs', 45, ...
        'widthDegs', 0.4, ...
        'contrast', 0.8, ...
        'meanLuminanceCdPerM2', 40);
    
    %% Step 3. Generate a scene describing the stimulus
    scene = generateStimulusScene(stimParams, presentationDisplay);
    
    %% Step 4. Realize this scene into a particular display
    realizedScene = realizeSceneInDisplay(scene, presentationDisplay);
    
    %% Step 5. Generate human optics
    opticalImage = oiCreate('wvf human');
    
    %% Step 6. Compute retinal image
    opticalImage = oiCompute(opticalImage, realizedScene);
    
    %% Step 7. Generate hexagonal cone mosaic
    coneMosaic = coneMosaicHex(7, ...
        'fovDegs', stimParams.widthDegs, ...
        'eccBasedConeDensity', true, ...
        'eccBasedConeQuantalEfficiency', true, ...
        'maxGridAdjustmentIterations', 50);
    
    %% Compute using a 10 ms integration time
    coneMosaic.integrationTime = 10/1000;
    
    %% Step 8. Compute isomerizations
    coneMosaic.compute(opticalImage, 'emPath', zeros(3, 1,2));
    
    
    %% Step 9. Display components
    displayedWavelengths = [450:25:725];
    scenePhotons = sceneGet(scene, 'photons');
    realizedScenePhotons = sceneGet(realizedScene, 'photons');
    maxPhotons = max([ max(scenePhotons(:)) max(realizedScenePhotons(:))]);
    
    %% Visualize scene 
    figNo = 1;
    visualizeScene(scene, maxPhotons, displayedWavelengths , figNo, 'scene');
    
    %% Visualize display
    figNo = figNo + 1;
    visualizeDisplay(presentationDisplay, figNo);

    %% Visualize scene rendered on display
    figNo = 3;
    visualizeScene(realizedScene, maxPhotons, displayedWavelengths, figNo, 'realized scene');
    
    %% Compare the 2 scenes
    figNo = 4;
    visualizeSpectralSlices(scene, realizedScene, maxPhotons, displayedWavelengths, figNo);
    
    
    %% Visualize the PSF
    figNo = 5;
    visualizePSF(opticalImage, displayedWavelengths, figNo);
    

    %% Visualize the optical image
    figNo = 6;
    visualizeOpticalImage(opticalImage, displayedWavelengths, figNo);
    

    size(coneMosaic.absorptions)
    %% Display results
    figure(7); clf;
    axHandle = subplot(1,2,1);
    coneMosaic.visualizeGrid('axesHandle', axHandle, 'displayVisualDegs', true);
    
    axHandle = subplot(1,2,2);
    coneMosaic.renderActivationMap(axHandle, coneMosaic.absorptions, 'mapType', 'modulated disks');
    
end
    
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

    if (1==2)
        dImage = displayCompute(display, rgbSettingsImage);
        size(dImage)
        imagescRGB(dImage);
        pause
    end
    
    % Generate a scene based on these RGB settings
    meanLuminance = [];
    realizedScene = sceneFromFile(rgbSettingsImage, 'rgb', meanLuminance, display);

    % Set the realized scene size and view distance to match those of the scene
    realizedScene = sceneSet(realizedScene, 'wangular', sceneGet(scene,'wangular'));
    realizedScene = sceneSet(realizedScene, 'distance', sceneGet(scene,'distance'));
end


% ===== VISUALIZATION HELPER METHODS ====

function visualizePSF(opticalImage, displayedWavelengths, figNo)

    % Get the optics data
    optics = oiGet(opticalImage, 'optics');
    
    % Get the psf 
    psf = opticsGet(optics, 'psf data');
    
    % Extract wavevelength support
    psfWavelengths = opticsGet(optics, 'otf wave');
    
    % Compute spatial support
    psfSampleSpacing = opticsGet(optics, 'psf spacing', 'um');
    psfSupport = (1:size(psf,1))*psfSampleSpacing;
    psfSupportMicrons = (psfSupport - mean(psfSupport));
    
    figure(figNo); clf;
    %% Display the PSFs
    for k = 1:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-psfWavelengths));
        imagesc(psfSupportMicrons, psfSupportMicrons, squeeze(psf(:,:,idx))); hold on;
        plot([0 0], [-20 20], 'g-'); plot([-20 20], [0 0],'g-'); axis 'xy';  axis 'image';
        set(gca, 'XLim', 20*[-1 1], 'YLim', 20*[-1 1], 'XTick', [-20:10:20], 'YTick', [-20:10:20]);
        if (k == 1)
            ylabel('retinal space (microns)');
        else
            set(gca, 'XTick', [], 'YTick', []);
        end
        title(sprintf('%d nm', psfWavelengths(idx)));
    end
    colormap(gray);
    drawnow;
end

function visualizeOpticalImage(opticalImage, displayedWavelengths, figNo)
    
    spatialSupport = oiGet(opticalImage, 'spatial support', 'mm');
    oiWavelengths = oiGet(opticalImage, 'wave');
    photons = oiGet(opticalImage, 'photons');
    rgbImage = oiGet(opticalImage, 'rgb image');
    figure(figNo); clf;
    
    subplot(3,4,1);
    image(spatialSupport(1,:,1),spatialSupport(:,1,2), rgbImage); axis 'xy';  axis 'image'
    xlabel('retinal space (mm)'); ylabel('retinal space (mm)');
    % Display more info
    title(sprintf('FOV: %2.2f degs', oiGet(opticalImage, 'wAngular')));
    
    photonRange = [0 max(photons(:))];
    for k = 2:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-oiWavelengths));
        photonsAtWavelengthBand = squeeze(photons(:,:,idx));
        imagesc(spatialSupport(1,:,1),spatialSupport(:,1,2), squeeze(photons(:,:,idx)), photonRange); 
        axis 'xy'; axis 'image'
        title(sprintf('%d nm, %2.3f E15 photons/sec', oiWavelengths(idx), mean(photonsAtWavelengthBand(:))*1E-15));
    end
    colormap(gray);
    drawnow;
end


function visualizeSpectralSlices(scene, realizedScene, maxPhotons, displayedWavelengths, figNo)

    luminanceMap = sceneGet(scene, 'luminance');
    [~, idx] = max(luminanceMap(:));
    [row,col] = ind2sub(size(luminanceMap), idx);
    
    scenePhotons = sceneGet(scene, 'photons')/1e15;
    realizedScenePhotons = sceneGet(realizedScene, 'photons')/1e15;
    
    figure(figNo); clf;
    sceneRadiance = squeeze(scenePhotons(row,col,:));
    realizedSceneRadiance = squeeze(realizedScenePhotons(row,col,:));
    plot(sceneGet(scene, 'wave'), sceneRadiance, 'ko-'); hold on;
    plot(sceneGet(realizedScene, 'wave'), realizedSceneRadiance, 'ro-');
    set(gca, 'YLim', [0 maxPhotons]*1e-15, 'XTick', displayedWavelengths, 'XLim', [370 790]);
    xlabel('wave');
    ylabel('photon rate (1e15 photons/pixel/sec)');
end


function visualizeScene(scene, maxPhotons, displayedWavelengths, figNo, sceneName)
    % Extract spatial support
    spatialSupportMilliMeters = sceneGet(scene, 'spatial support', 'mm');
    % Extract wavevelength support
    wavelengths = sceneGet(scene, 'wave');
    % Extact field of view (degs) and viewing distance (meters)
    fovDegrees = sceneGet(scene, 'wAngular');
    viewingDistanceMeters = sceneGet(scene, 'distance');
    % Exctract luminance map
    luminanceMap = sceneGet(scene, 'luminance');
    % Exctract chromaticity map
    XYZmap = sceneGet(scene, 'xyz');
    
    xMap = squeeze(XYZmap(:,:,1))./sum(XYZmap,3);
    yMap = squeeze(XYZmap(:,:,2))./sum(XYZmap,3);
    % Compute mean luminance and mean chromaticity
    meanLuminance = mean(luminanceMap(:));
    meanChromaticity = [mean(xMap(:)) mean(yMap(:))];
    
    % Extract scene radiance (photon rate)
    photons = sceneGet(scene, 'photons');
    
    % Extract scene as linear RGB primary values
    rgbImage = sceneGet(scene, 'rgb image');
    
    hFig = figure(figNo); clf;
    set(hFig, 'Name', sceneName);
    subplot(3,4,1);
    image(spatialSupportMilliMeters(1,:,1),spatialSupportMilliMeters(:,1,2), rgbImage); 
    axis 'xy'; axis 'image'
    xlabel('visual space (mm)'); ylabel('visual space (mm)');
    % Display more info
    title(sprintf('FOV: %2.2f degs, VD: %2.2f meters\n mean lum.: %2.1f cd/m2, mean xy: (%0.3f,%0.3f)', ...
        fovDegrees, viewingDistanceMeters, meanLuminance, meanChromaticity(1), meanChromaticity(2)));
    
  
    photonRange = [0 maxPhotons];  
    for k = 2:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-wavelengths));
        photonsAtWavelengthBand = squeeze(photons(:,:,idx));
        imagesc(spatialSupportMilliMeters(1,:,1),spatialSupportMilliMeters(:,1,2), ...
            photonsAtWavelengthBand, photonRange); 
        axis 'xy';  axis 'image';
        title(sprintf('%d nm, %2.3f E15 photons/sec', wavelengths(idx), mean(photonsAtWavelengthBand(:))*1E-15 ));
        set(gca, 'XTick', [], 'YTick', []);
    end
     colormap(gray);
end

function visualizeDisplay(display, figNo)
    % Retrieve wavelength support and spectral power distribution of the
    % primaries
    wave = displayGet(display, 'wave');
    spds = displayGet(display, 'spd primaries');
    % Retrieve resolution in dots/inch
    resolutionDPI = displayGet(display, 'dpi');
    figure(figNo); clf;
    subplot(1,2,1);
    plot(wave, spds(:,1), 'ro-'); hold on;
    plot(wave, spds(:,2), 'go-'); plot(wave, spds(:,3), 'bo-'); 
    xlabel('wavelength (nm)'); ylabel('energy');
    title(sprintf('resolution: %2.2f dots/inch', resolutionDPI));
    subplot(1,2,2);
    gammaTable = displayGet(display, 'gamma table');
    plot(1:size(gammaTable,1), gammaTable(:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(1:size(gammaTable,1), gammaTable(:,2), 'g-', 'LineWidth', 1.5); hold on;
    plot(1:size(gammaTable,1), gammaTable(:,3), 'b-', 'LineWidth', 1.5); hold on;
    drawnow;
end

