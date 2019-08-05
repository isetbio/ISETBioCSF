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
    set(hFig, 'Name', sceneName, 'Position', [10 10 1300 1024], 'Color', [1 1 1]);
    subplot(3,4,1);
    image(spatialSupportMilliMeters(1,:,1),spatialSupportMilliMeters(:,1,2), rgbImage); 
    axis 'xy'; axis 'image'
    set(gca, 'XTick', [-10:1:10], 'YTick', [-10:1:10]);
    xlabel('visual space (mm)'); ylabel('visual space (mm)');
    xtickformat('%0.2f'); ytickformat('%0.2f');
    set(gca, 'FontSize', 14);
    % Display more info
    title(sprintf('FOV: %2.2f degs\nVD: %2.2f meters\n mean lum.: %2.1f cd/m2\nmean xy: (%0.3f,%0.3f)', ...
        fovDegrees, viewingDistanceMeters, meanLuminance, meanChromaticity(1), meanChromaticity(2)));
    
  
    photonRange = [0 maxPhotons];  
    for k = 2:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-wavelengths));
        photonsAtWavelengthBand = squeeze(photons(:,:,idx));
        imagesc(spatialSupportMilliMeters(1,:,1),spatialSupportMilliMeters(:,1,2), ...
            photonsAtWavelengthBand, photonRange); 
        axis 'xy';  axis 'image';
        set(gca, 'XTick', [], 'YTick', [], 'FontSize', 14);
        title(sprintf('%d nm\n%2.3fE15 photons/sample/sec', wavelengths(idx), mean(photonsAtWavelengthBand(:))*1E-15 ));
    end
     colormap(gray);
end