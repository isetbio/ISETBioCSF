function simpleTutorial

    %% Generate a scene representing a 10 c/deg Gabor stimulus
    imagePixels = 300; viewingDistanceInMeters = 1.0;
    stimFrequencyCyclerPerDegree = 16;
    stimOrientationDegs = 60;
    stimSizeDegs = 0.4;
    stimParams = struct(...
        'freq', stimSizeDegs*stimFrequencyCyclerPerDegree, ...
        'ang', stimOrientationDegs/180*pi, ...
        'contrast', 0.8, ...
        'GaborFlag', 0.2, ...
        'row', imagePixels , 'col', imagePixels);
    scene = sceneCreate('harmonic', stimParams);
    scene = sceneSet(scene, 'distance', viewingDistanceInMeters);
    scene = sceneSet(scene, 'wangular', stimSizeDegs);
    
    %% Display scene as RGB
    % Obtain support and rgb data 
    spatialSupport = sceneGet(scene, 'spatial support', 'mm');
    rgbImage = sceneGet(scene, 'rgb image');
    figure(1); clf;
    subplot(3,4,1);
    image(spatialSupport(1,:,1),spatialSupport(:,1,2), rgbImage); axis 'xy'; axis 'image'
    xlabel('visual space (mm)'); ylabel('visual space (mm)');
    % Display more info
    title(sprintf('FOV: %2.2f degs\nVD: %2.2f meters', sceneGet(scene, 'wAngular'), sceneGet(scene, 'distance')));
    
    %% Generate a display for realizing this scene
    renderingDisplay = displayCreate('LCD-Apple.mat');
    sceneWidthInMeters = sceneGet(scene, 'width');
    scenePixelSizeInMeters = sceneWidthInMeters/imagePixels;
    sceneDotsPerMeter = 1/scenePixelSizeInMeters;
    sceneDotsPerInch = sceneDotsPerMeter / 39.37;
    renderingDisplay = displaySet(renderingDisplay, 'dpi', sceneDotsPerInch);
    renderingDisplay = displaySet(renderingDisplay, 'viewing distance', viewingDistanceInMeters);
    wave = displayGet(renderingDisplay, 'wave');
    spds = displayGet(renderingDisplay, 'spd primaries');
    figure(2); clf;
    plot(wave, spds(:,1), 'ro-'); hold on;
    plot(wave, spds(:,2), 'go-'); plot(wave, spds(:,3), 'bo-'); 
    xlabel('wavelength (nm)'); ylabel('energy');
    drawnow;
    
    %% Realize this scene into a particular display
    scene = sceneFromFile(rgbImage, 'rgb', [], renderingDisplay);

%     spatialSupport2 = sceneGet(scene, 'spatial support', 'mm');
%     rgbImage = sceneGet(scene, 'rgb image');
%     figure(11); clf;
%     subplot(3,4,1);
%     image(spatialSupport2(1,:,1),spatialSupport2(:,1,2), rgbImage); axis 'image'
%     xlabel('space (mm)'); ylabel('space (mm)');
    
    
    %% Display the emitted radiance for this scene
    photons = sceneGet(scene, 'photons');
    wavelengths = sceneGet(scene, 'wave');
    displayedWavelengths = [450:25:725];
    photonRange = [min(photons(:)) max(photons(:))];
    figure(1);
    for k = 2:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-wavelengths));
        imagesc(spatialSupport(1,:,1),spatialSupport(:,1,2), squeeze(photons(:,:,idx)), photonRange); axis 'xy';  axis 'image';
        title(sprintf('%d nm', wavelengths(idx)));
        set(gca, 'XTick', [], 'YTick', []);
    end
    colormap(gray);
    drawnow;
    
    %% Generate human optics
    opticalImage = oiCreate('wvf human');
    
    %% Compute retinal image
    opticalImage = oiCompute(opticalImage, scene);
    
    %% Extract useful information
    % Wavelenghts
    oiWavelengths = oiGet(opticalImage, 'wave');
    % Get the psf  data
    optics = oiGet(opticalImage, 'optics');
    psf = opticsGet(optics, 'psf data');
    psfWavelengths = opticsGet(optics, 'otf wave');
    psfSampleSpacing = opticsGet(optics, 'psf spacing', 'um');
    psfSupport = (1:size(psf,1))*psfSampleSpacing;
    psfSupportMicrons = (psfSupport - mean(psfSupport));
    
    figure(3); clf;
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
    
    %% Display the optical image
    spatialSupport = oiGet(opticalImage, 'spatial support', 'mm');
    rgbImage = oiGet(opticalImage, 'rgb image');
    figure(4); clf;
    subplot(3,4,1);
    image(spatialSupport(1,:,1),spatialSupport(:,1,2), rgbImage); axis 'xy';  axis 'image'
    xlabel('retinal space (mm)'); ylabel('retinal space (mm)');
    % Display more info
    title(sprintf('FOV: %2.2f degs', oiGet(opticalImage, 'wAngular')));
    
    for k = 2:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-oiWavelengths));
        imagesc(spatialSupport(1,:,1),spatialSupport(:,1,2), squeeze(photons(:,:,idx)), photonRange); axis 'image';
        set(gca, 'XTick', [], 'YTick', []);
        title(sprintf('%d nm', oiWavelengths(idx)));
    end
    colormap(gray);
    drawnow;
    
    %% Generate hexagonal cone mosaic
    coneMosaic = coneMosaicHex(7, ...
        'fovDegs', stimSizeDegs, ...
        'eccBasedConeDensity', true, ...
        'maxGridAdjustmentIterations', 50);
    
    %% Compute isomerizations
    isomerizations = coneMosaic.compute(opticalImage);
    
    %% Display results
    figure(5); clf;
    axHandle = subplot(1,2,1);
    coneMosaic.visualizeGrid('axesHandle', axHandle, 'displayVisualDegs', true);
    
    axHandle = subplot(1,2,2);
    coneMosaic.renderActivationMap(axHandle, isomerizations);
    
end