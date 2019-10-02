function visualizeSpectralSlices(scene, realizedScene, maxPhotons, ...
    displayedWavelengths, figNo)
% Visualize spectral slices
%
% Syntax:
%   visualizeSpectralSlices(scene, realizedScene, maxPhotons, ...
%       displayedWavelengths, figNo)
%
% Description:
%    A function to visualize the spectral scenes.
%
% Inputs:
%    scene         - Struct. A scene structure.
%    realizedScene - Struct. The realized scene's structure.
%    maxPhotons    - Numeric. The maximum number of photons.
%    displayedWavelengths
%                  - Vector. A numeric vector of the displayed wavelengths.
%    figNo         - Numeric. The figure number.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    luminanceMap = sceneGet(scene, 'luminance');
    [~, idx] = max(luminanceMap(:));
    [row,col] = ind2sub(size(luminanceMap), idx);
    
    scenePhotons = sceneGet(scene, 'photons') / 1e15;
    realizedScenePhotons = sceneGet(realizedScene, 'photons')/1e15;
    
    hFig = figure(figNo);
    clf;
    set(hFig, 'Position', [10 10 780 580], 'Color', [1 1 1]);
    sceneRadiance = squeeze(scenePhotons(row,col,:));
    realizedSceneRadiance = squeeze(realizedScenePhotons(row,col,:));
    plot(sceneGet(scene, 'wave'), sceneRadiance, 'ko-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], ...
        'LineWidth', 1.5);
    hold on;
    plot(sceneGet(realizedScene, 'wave'), realizedSceneRadiance, 'ro-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], ...
        'LineWidth', 1.5);
    legend({'scene pixel (0,0)','realized scene pixel (0,0)'});
    set(gca, 'YLim', [0 maxPhotons] * 1e-15, ...
        'XTick', displayedWavelengths, 'XLim', [370 790], 'FontSize', 14);
    xlabel('wavelength (nm)');
    ylabel('photon rate (1e15 photons/sample/sec)');
end