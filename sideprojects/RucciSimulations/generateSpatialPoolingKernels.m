function spatialPoolingKernels = generateSpatialPoolingKernels(theMosaic, ...
            lowFrequencyTemplate, lowFrequencyTemplateOrtho, ...
            highFrequencyTemplate, highFrequencyTemplateOrtho, ...
            lowFrequencyScene, ...
            lowFrequencySceneOrtho, ...
            highFrequencyScene, ...
            highFrequencySceneOrtho, ...
            spatialSupportDegs)
    
    coneLocsMeters = theMosaic.coneLocsHexGrid;
    coneLocsDegs(:,1) = coneLocsMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsDegs(:,2) = coneLocsMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    
    [X,Y] = meshgrid(spatialSupportDegs, spatialSupportDegs);
    rfCoordsDegs = [X(:) Y(:)];
    
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
    
    
    spatialPoolingKernels.lowFrequencyPoolingWeightsLinear           = lowFrequencyTemplate.linear(idx);
    spatialPoolingKernels.lowFrequencyPoolingWeightsQuadrature       = lowFrequencyTemplate.quadrature(idx);

    spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsLinear      = lowFrequencyTemplateOrtho.linear(idx);
    spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsQuadrature  = lowFrequencyTemplateOrtho.quadrature(idx);
    
    spatialPoolingKernels.highFrequencyPoolingWeightsLinear          = highFrequencyTemplate.linear(idx);
    spatialPoolingKernels.highFrequencyPoolingWeightsQuadrature      = highFrequencyTemplate.quadrature(idx);

    spatialPoolingKernels.highFrequencyOrthoPoolingWeightsLinear     = highFrequencyTemplateOrtho.linear(idx);
    spatialPoolingKernels.highFrequencyOrthoPoolingWeightsQuadrature = highFrequencyTemplateOrtho.quadrature(idx);
    
    
    coneRadiusMicrons = 1.4;
    coneRadiusDegs = coneRadiusMicrons/theMosaic.micronsPerDegree;
    spatialPoolingKernels.coneAperture(:,1) = coneRadiusDegs * cosd(0:60:360);
    spatialPoolingKernels.coneAperture(:,2) = coneRadiusDegs * sind(0:60:360);
    spatialPoolingKernels.coneLocsDegs = coneLocsDegs;



    hFig = figure(333); clf;
    set(hFig, 'Position', [10 10 1640 910], 'Color', [1 1 1]);


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 3, ...
        'colsNum', 4, ...
        'heightMargin',  0.03, ...
        'widthMargin',    0.03, ...
        'leftMargin',     0.03, ...
        'rightMargin',    0.03, ...
        'bottomMargin',   0.05, ...
        'topMargin',      0.05);

    % The low frequency stimulus and templates
    ax1 = subplot('Position', subplotPosVectors(1,1).v);
    ax2 = subplot('Position', subplotPosVectors(2,1).v);
    ax3 = subplot('Position', subplotPosVectors(3,1).v);
    sceneRGB = sceneGet(lowFrequencyScene, 'rgbimage');
    visualizeSpatialPoolingScheme(ax1, ax2, ax3, spatialPoolingKernels.coneLocsDegs, spatialPoolingKernels.coneAperture, ...
        spatialPoolingKernels.lowFrequencyPoolingWeightsLinear, ...
        spatialPoolingKernels.lowFrequencyPoolingWeightsQuadrature, ...
        spatialSupportDegs, sceneRGB, ...
        'low frequency stimulus (0 deg)');

    % The orthogonal low frequency stimulus and templates
    ax1 = subplot('Position', subplotPosVectors(1,2).v);
    ax2 = subplot('Position', subplotPosVectors(2,2).v);
    ax3 = subplot('Position', subplotPosVectors(3,2).v);
    sceneRGB = sceneGet(lowFrequencySceneOrtho, 'rgbimage');
    visualizeSpatialPoolingScheme(ax1, ax2, ax3, spatialPoolingKernels.coneLocsDegs, spatialPoolingKernels.coneAperture, ...
        spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsLinear, ...
        spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsQuadrature, ...
        spatialSupportDegs, sceneRGB, ...
        'low frequency stimulus (90 deg)');


    % The high frequency stimulus and templates
    ax1 = subplot('Position', subplotPosVectors(1,3).v);
    ax2 = subplot('Position', subplotPosVectors(2,3).v);
    ax3 = subplot('Position', subplotPosVectors(3,3).v);
    sceneRGB = sceneGet(highFrequencyScene, 'rgbimage');
    visualizeSpatialPoolingScheme(ax1, ax2, ax3, spatialPoolingKernels.coneLocsDegs, spatialPoolingKernels.coneAperture, ...
        spatialPoolingKernels.highFrequencyPoolingWeightsLinear, ...
        spatialPoolingKernels.highFrequencyPoolingWeightsQuadrature, ...
        spatialSupportDegs, sceneRGB, ...
        'high frequency stimulus (0 deg)');

    % The orthogonal low frequency stimulus and templates
    ax1 = subplot('Position', subplotPosVectors(1,4).v);
    ax2 = subplot('Position', subplotPosVectors(2,4).v);
    ax3 = subplot('Position', subplotPosVectors(3,4).v);
    sceneRGB = sceneGet(highFrequencySceneOrtho, 'rgbimage');
    visualizeSpatialPoolingScheme(ax1, ax2, ax3, spatialPoolingKernels.coneLocsDegs, spatialPoolingKernels.coneAperture, ...
        spatialPoolingKernels.highFrequencyOrthoPoolingWeightsLinear, ...
        spatialPoolingKernels.highFrequencyOrthoPoolingWeightsQuadrature, ...
        spatialSupportDegs, sceneRGB, ...
        'high frequency stimulus (90 deg)');

end

function visualizeSpatialPoolingScheme(ax1, ax2, ax3, coneLocsDegs, coneAperture, ...
    poolingWeights, poolingWeightsQuadr, spatialSupportDegs, stimulusRGB, stimulusTitle)

    % Plot stimulus
    imagesc(ax1, spatialSupportDegs, spatialSupportDegs, stimulusRGB);
    set(ax1, 'CLim', [-1 1], 'XLim', max(abs(spatialSupportDegs))*[-1 1], 'YLim', max(abs(spatialSupportDegs))*[-1 1], 'FontSize', 14);
    ylabel(ax1,'spatial position (degs)');
    axis(ax1, 'image');
    axis(ax1, 'xy');
    set(ax1, 'XTickLabel', {});
    ticks = max(abs(spatialSupportDegs))*[-1 0 1];
    tickLabels = sprintf('%2.1f\n', ticks);
    set(ax1, 'YTick', ticks, 'YTickLabel', tickLabels);
    colormap(ax1,gray(1024));
    title(ax1, stimulusTitle);
    drawnow;

    visualizePoolingKernel(ax2, coneLocsDegs, coneAperture, poolingWeights, spatialSupportDegs, false);
    visualizePoolingKernel(ax3, coneLocsDegs, coneAperture, poolingWeightsQuadr, spatialSupportDegs, true);
end

function  visualizePoolingKernel(ax, coneLocsDegs, coneAperture, poolingWeights, spatialSupportDegs, showXLabel)
    hold(ax, 'on');
    
    quantizationLevels = 1024;
    cmap = brewermap(quantizationLevels, '*RdBu');
    faceColorsNormalizedValues = round(quantizationLevels * (1+poolingWeights));

    edgeColor = [0 0 0];
    lineWidth = 0.1;
    renderWeights(ax, coneAperture, coneLocsDegs(:,1), coneLocsDegs(:,2), ...
        faceColorsNormalizedValues, edgeColor, lineWidth);

    ticks = max(abs(spatialSupportDegs))*[-1 0 1];
    tickLabels = sprintf('%2.1f\n', ticks);
    set(ax, 'YTick', ticks, 'YTickLabel', tickLabels);
    
    set(ax, 'XLim', max(abs(spatialSupportDegs))*[-1 1], 'YLim', max(abs(spatialSupportDegs))*[-1 1], 'FontSize', 14);
    if (showXLabel)
        xlabel(ax,'spatial position (degs)');
        set(ax, 'XTick', ticks, 'XTickLabel', tickLabels);
    else
        set(ax, 'XTickLabel', {});
    end
    set(ax, 'YTick', ticks, 'YTickLabel', tickLabels);
    ylabel(ax,'spatial position (degs)');
    axis(ax, 'square');
    box(ax, 'on');
    axis(ax, 'xy');
    colormap(ax, cmap);
    drawnow;
    
end

function renderWeights(axesHandle, coneAperture, xCoords, yCoords, ...
    faceColorsNormalizedValues, edgeColor, lineWidth)

    verticesPerCone = size(coneAperture,1);
    
    verticesList = zeros(verticesPerCone * numel(xCoords), 2);
    facesList = [];
    colors = [];
    
    for coneIndex = 1:numel(xCoords)
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        
        verticesList(idx, 1) = coneAperture(:,1) + xCoords(coneIndex);
        verticesList(idx, 2) = coneAperture(:,2) + yCoords(coneIndex);
        facesList = cat(1, facesList, idx);
        colors = cat(1, colors, ...
            repmat(faceColorsNormalizedValues(coneIndex), ...
            [verticesPerCone 1]));
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end