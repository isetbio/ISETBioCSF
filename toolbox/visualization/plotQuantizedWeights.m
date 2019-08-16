function plotQuantizedWeights(axesHandle, quantizedWeights, quantizationLevels, coneLocsInDegs, pixelOutline)

    quantizedWeights(quantizedWeights >  1) = 1;
    quantizedWeights(quantizedWeights < -1) = -1;

    faceColorsNormalizedValues = 0.5*(1+quantizedWeights);
    xCoords = coneLocsInDegs(:,1);
    yCoords = coneLocsInDegs(:,2);
    lineWidth = 0.2;
    edgeColor = [0 0 0];
    renderWeights(axesHandle, pixelOutline, xCoords, yCoords, ...
        faceColorsNormalizedValues, edgeColor, lineWidth);
    
end

function renderWeights(axesHandle, pixelOutline, xCoords, yCoords, ...
    faceColorsNormalizedValues, edgeColor, lineWidth)

    
    facesList = [];
    colors = [];
    if (size(pixelOutline.x,1) == 1)
        pixelOutline.x = repmat(pixelOutline.x, [numel(xCoords) 1]);
        pixelOutline.y = repmat(pixelOutline.y, [numel(xCoords) 1]);
    end

    for coneIndex = 1:numel(xCoords)
        coneAperture(:,1) = pixelOutline.x(coneIndex,:);
        coneAperture(:,2) = pixelOutline.y(coneIndex,:);
        if (coneIndex == 1)
            verticesPerCone = size(coneAperture,1);
            verticesList = zeros(verticesPerCone * numel(xCoords), 2);
        end
        
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