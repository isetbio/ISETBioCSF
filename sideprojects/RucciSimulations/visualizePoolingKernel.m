function  visualizePoolingKernel(ax, coneLocsDegs, coneAperture, poolingWeights, spatialSupportDegs, showXLabel)
    
    poolingWeights = poolingWeights / max(abs(poolingWeights));
    
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