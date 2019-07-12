function plotQuantizedWeights(axesHandle, quantizedWeights, quantizationLevels, coneLocsInDegs, coneX, coneY)

    quantizedWeights(quantizedWeights >  1) = 1;
    quantizedWeights(quantizedWeights < -1) = -1;

    faceColorsNormalizedValues = (1+quantizedWeights);
    coneAperture(:,1) = coneX;
    coneAperture(:,2) = coneY;
    xCoords = coneLocsInDegs(:,1);
    yCoords = coneLocsInDegs(:,2);
    lineWidth = 0.2;
    edgeColor = [0 0 0];
    renderWeights(axesHandle, coneAperture, xCoords, yCoords, ...
        faceColorsNormalizedValues, edgeColor, lineWidth);
    
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