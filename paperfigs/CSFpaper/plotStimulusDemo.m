function plotStimulusDemo

    spatialParams = spatialParamsGenerate();
    backgroundParams = backgroundParamsGenerate();
    backgroundParams.backgroundxyY = [0.33 0.33 45];
    colorModulationParams = colorModulationParamsGenerate();
    colorModulationParams.coneContrasts = [0.75 0.75 0.75];
    
    oiParams = oiParamsGenerate();
    [theScene, gamutScaleFactor, varargout] = ...
        colorSceneCreate(spatialParams,backgroundParams,colorModulationParams,oiParams,false)

    rgb = xyz2rgb(0.012*sceneGet(theScene, 'xyz'));
    %rgb = sum(rgb,3);
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 600 600]);
    subplot('Position', [0.02 0.02 0.96 0.96]);
    
    image(rgb);
    colormap(gray(1024));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    box on;
    
    localDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    NicePlot.exportFigToPDF(sprintf('%s/TestStimulus.pdf', localDir), hFig, 300);

    
end

