function visualizeSceneLuminanceProfiles(nullScene, lowFrequencyScenes, highFrequencyScenes, contrastLevels, noiseInstances)

    nContrasts = numel(contrastLevels);
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:noiseInstances 
            realizedScene = lowFrequencyScenes{theContrastLevel, theInstance};
            lumMap1(theInstance, theContrastLevel,:,:) = sceneGet(realizedScene, 'luminance');

            realizedScene = highFrequencyScenes{theContrastLevel, theInstance};
            lumMap2(theInstance, theContrastLevel,:,:) = sceneGet(realizedScene, 'luminance');
        end
    end


    fovDegs = sceneGet(realizedScene, 'horizontal fov');
    cols = sceneGet(realizedScene, 'cols');

    sampleSpacing = fovDegs / cols;
    spatialSupportDegs = (1:cols)*sampleSpacing;
    spatialSupportDegs = spatialSupportDegs - mean(spatialSupportDegs);
        
    nContrasts = size(lumMap1,2);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', nContrasts, ...
       'colsNum', 2, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.04);
   
    
    N = round(size(lumMap1,3)/2);
    CLim = [0 max([max(lumMap1(:)) max(lumMap2(:))])];
    

    figure(444); clf;
    colormap(gray);
    lumMap  = sceneGet(nullScene, 'luminance');
    imagesc(spatialSupportDegs, spatialSupportDegs, lumMap);
    set(gca, 'CLim', CLim);
    axis 'image';
    title('null scene');
    
    figure(555);
    colormap(gray);
    for theContrastLevel = 1:nContrasts 
        subplot('Position', subplotPosVectors(theContrastLevel,1).v);
        imagesc(spatialSupportDegs, spatialSupportDegs, squeeze(lumMap1(1,theContrastLevel,:,:)));
        set(gca, 'CLim', CLim);
        axis 'image';
        
        
        subplot('Position', subplotPosVectors(theContrastLevel,2).v);
        meanSpatialProfile = mean(squeeze(lumMap1(:, theContrastLevel, N, :)),1);
        plot(spatialSupportDegs, meanSpatialProfile, 'rs-');
        axis 'square'
        set(gca, 'YLim', CLim);
        drawnow
    end
    
    figure(666);
    colormap(gray);
    for theContrastLevel = 1:nContrasts 
        subplot('Position', subplotPosVectors(theContrastLevel,1).v);
        imagesc(spatialSupportDegs, spatialSupportDegs, squeeze(lumMap2(1,theContrastLevel,:,:)));
        set(gca, 'CLim', CLim);
        axis 'image';
        
        subplot('Position', subplotPosVectors(theContrastLevel,2).v);
        meanSpatialProfile = mean(squeeze(lumMap2(:, theContrastLevel, N, :)),1);
        plot(spatialSupportDegs, meanSpatialProfile, 'rs-');
        axis 'square'
        set(gca, 'YLim', CLim);
        drawnow
    end
    
end