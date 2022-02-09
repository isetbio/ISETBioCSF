function visualizeOpticalImageAtWavelengths(opticalImage, displayedWavelengths, figNo)
    
    spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
    
    % Convert spatial support in degrees
    optics = oiGet(opticalImage, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportXDegs = spatialSupportDegs(1,:,1);
    spatialSupportYDegs = spatialSupportDegs(:,1,2);
    
    oiWavelengths = oiGet(opticalImage, 'wave');
    photons = oiGet(opticalImage, 'photons');
    rgbImage = oiGet(opticalImage, 'rgb image');
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1300 1024], 'Color', [1 1 1], 'Name', 'Retinal Image');
    subplot(3,4,1);
    image(spatialSupportXDegs,spatialSupportYDegs, rgbImage); axis 'xy';  axis 'image'
    set(gca, 'XTick', [-1:0.1:11], 'YTick', [-1:0.1:1]);
    xtickformat('%0.1f'); ytickformat('%0.1f');
    xlabel('retinal space (degs)'); ylabel('retinal space (degs)');
    set(gca, 'FontSize', 14);
    % Display more info
    title(sprintf('FOV: %2.2f degs', oiGet(opticalImage, 'wAngular')));
    
    photonRange = [0 max(photons(:))];
    for k = 2:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-oiWavelengths));
        photonsAtWavelengthBand = squeeze(photons(:,:,idx));
        imagesc(spatialSupportXDegs,spatialSupportYDegs, squeeze(photons(:,:,idx)), photonRange); 
        set(gca, 'XTick', [], 'YTick', [], 'FontSize', 14);
        axis 'xy'; axis 'image'
        title(sprintf('%d nm\n%2.3fE15 photons/sample/sec', oiWavelengths(idx), mean(photonsAtWavelengthBand(:))*1E-15));
    end
    colormap(gray);
    drawnow;
end