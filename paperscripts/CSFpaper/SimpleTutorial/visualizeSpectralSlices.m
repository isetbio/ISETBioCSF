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