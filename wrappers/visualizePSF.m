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
    
    % Convert spatial support in minutes of arc
    focalLength = opticsGet(optics, 'focal length');
    micronsPerDegree = focalLength*tand(1)*1e6;
    psfSupportArcMin = psfSupportMicrons/micronsPerDegree*60;
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1000 900], 'Color', [1 1 1]);
    psfRangeArcMin = 10*[-1 1];
    %% Display the PSFs
    for k = 1:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-psfWavelengths));
        imagesc(psfSupportArcMin, psfSupportArcMin, squeeze(psf(:,:,idx))); hold on;
        plot([0 0], psfRangeArcMin, 'g-'); plot(psfRangeArcMin, [0 0],'g-'); axis 'xy';  axis 'image';
        set(gca, 'XLim', psfRangeArcMin, 'YLim', psfRangeArcMin, 'XTick', [-20:5:20], 'YTick', [-20:5:20], 'FontSize', 14);
        if (k == 1)
            xlabel('retinal space (arc min.)');
            ylabel('retinal space (arc min.)');
        else
            set(gca, 'XTick', [], 'YTick', []);
        end
        title(sprintf('%d nm', psfWavelengths(idx)));
    end
    colormap(gray);
    drawnow;
end