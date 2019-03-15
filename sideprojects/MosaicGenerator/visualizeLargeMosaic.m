function visualizeLargeMosaic()
    load('ConeMosaic_15.0Degs_Best.mat', 'theMosaic');

    center = [0 0];
    radius = 1.5;
    
    visualizeSubregion(theMosaic, center, radius);
end

function visualizeSubregion(theMosaic, center, radius)

    hexCoords = theMosaic.coneLocsHexGrid;
    
    metersToDegs = 1e6/300;
    hexCoordsDegs = hexCoords * metersToDegs;
    
    
    rr = bsxfun(@minus, hexCoordsDegs, center);
    rr = sqrt(sum(rr.^2,2));
    
    idx = find(rr < radius);
    numel(idx)
    
    meshEdgeColor = [0.5 0.5 0.5];
    meshFaceColor = [1.0 0.9 0.9];
    meshFaceAlpha = 0.5;
    meshEdgeAlpha = 0.8;
    lineStyle = '-';
    figure()
    axesHandle = gca;
    theMosaic.renderHexMesh(axesHandle, hexCoordsDegs(idx,1), hexCoordsDegs(idx,2), ...
        meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    hold(axesHandle, 'on');
    edgeColor = [0 0 0];
    faceColor = [0.6 0.6 0.6];
    lineWidth = 0.5;
    
    dxInner = diameterForCircularApertureFromWidthForSquareAperture(...
        theMosaic.pigment.pdWidth);
    dxOuter = diameterForCircularApertureFromWidthForSquareAperture(...
         theMosaic.pigment.width);
    iTheta = (0:10:360) / 180 * pi;
    outerApertureOutline.x = dxOuter / 2.0 * cos(iTheta);
    outerApertureOutline.y = dxOuter / 2.0 * sin(iTheta);
    innerApertureOutline.x = dxInner / 2.0 * cos(iTheta);
    innerApertureOutline.y = dxInner / 2.0 * sin(iTheta);
    
    [innerApertureOutlineVarying, outerApertureOutlineVarying] = theMosaic.computeApertureSizes(...
            dxInner, dxOuter, ...
            innerApertureOutline, outerApertureOutline, ...
            hexCoords(idx,1), hexCoords(idx,2) ...
        );
    
    innerApertureOutlineVarying.x = innerApertureOutlineVarying.x * metersToDegs;
    innerApertureOutlineVarying.y = innerApertureOutlineVarying.y * metersToDegs;
    
    theMosaic.renderPatchArray(axesHandle, innerApertureOutlineVarying, ...
        (hexCoordsDegs(idx,1))', hexCoordsDegs(idx,2), edgeColor, ...
        faceColor, lineStyle, lineWidth);
    axis 'equal';
    
end

