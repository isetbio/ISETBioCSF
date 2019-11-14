function generateRGCmosaicFromConeMosaic()

    % Instantiate a WatsonRGCModel object. Set the 'generateAllFigures' flag 
    % to true to generate several figures of the Watson 2014 paper
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    meridians = {'superior', 'inferior', 'nasal', 'temporal'};
    meridianColors = [0 0 1; 0 0 0; 0 0.8 0; 1 0 0];
    
    horizontalMeridians = {'nasal', 'temporal'};
    verticalMeridians = {'superior', 'inferior'};
    
    eccDegs = 0.2:0.1:8;
    midgetRGCRFtoConeRatioHorizontal = zeros(1,numel(eccDegs));
    midgetRGCRFtoConeRatioVertical = zeros(1,numel(eccDegs));
    for k = 1:numel(meridians)
        meridianName = sprintf('%s meridian', meridians{k});
        midgetRGCRFDensity = WatsonRGCCalc.midgetRGCRFDensity(eccDegs, meridianName, 'RFs per deg2');
        [coneSpacingDegs, coneDensity] = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, meridianName, 'Cones per deg2');
        midgetRGCRFtoConeRatio(k,:) = (0.5*midgetRGCRFDensity) ./ coneDensity;
        if (ismember(meridians{k}, horizontalMeridians))
            midgetRGCRFtoConeRatioHorizontal = midgetRGCRFtoConeRatioHorizontal + midgetRGCRFtoConeRatio(k,:);
        else
            midgetRGCRFtoConeRatioVertical = midgetRGCRFtoConeRatioVertical + midgetRGCRFtoConeRatio(k,:);
        end
    end
    
    figure(2); clf; 
    subplot(1,2,1);
    hold on;
    for k = 1:numel(meridians)
        conesPerMidgetRGC = 1./midgetRGCRFtoConeRatio(k,:);
        plot(eccDegs, conesPerMidgetRGC, '-', 'LineWidth', 1.5, 'Color', meridianColors(k,:)); hold on;
    end
    legend(meridians, 'Location', 'NorthWest');
    set(gca, 'FontSize', 14, 'XLim', [0 8], 'YLim', [1 5], 'YTick', 1:0.5:10, 'XTick', [0:1:10]);
    grid on; box on
    xlabel('ecc (degs)');
    ylabel('cones per midget RGC RF center');
    
    subplot(1,2,2);
    hold on;
    conesPerMidgetRGCHorizontalMeridian = 1./(0.5*midgetRGCRFtoConeRatioHorizontal);
    conesPerMidgetRGCHorizontalVertical = 1./(0.5*midgetRGCRFtoConeRatioVertical);
    plot(eccDegs, conesPerMidgetRGCHorizontalMeridian, 'k-', 'LineWidth', 1.5); hold on;
    plot(eccDegs, conesPerMidgetRGCHorizontalVertical, 'k--', 'LineWidth', 1.5); hold on;
    legend({'horizontal', 'vertical'}, 'Location', 'NorthWest');
    set(gca, 'FontSize', 14, 'XLim', [0 8], 'YLim', [1 5], 'YTick', 1:0.5:10, 'XTick', [0:1:10]);
    grid on; box on
    xlabel('ecc (degs)');
    ylabel('cones per midget RGC RF center');
    
    pause
    
    % Load cone mosaic
    load('HexConeMosaic15DegsBest.mat', 'theMosaic');
    
    % Visualize only part of the mosaic X(degs): [0 6], Y(degs): [-0.25 0.25]
    visualizedFOV = struct('xo', 3.0, 'yo', 0.0, 'width', 6, 'height', 0.5);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2800 600]);
    axesHandle  = subplot('Position', [0.02 0.5 0.97 0.45]);
    theMosaic.visualizeGrid('axesHandle', axesHandle, 'ticksInVisualDegs', true, 'visualizedFOV', visualizedFOV);
end


function visualizeSubregion(theMosaic, coneTypesHexGrid, roi)
    metersToDegs = 1e6/theMosaic.micronsPerDegree;
    hexCoordsMeters = theMosaic.coneLocsHexGrid;
    hexCoordsDegs = hexCoordsMeters * metersToDegs;
    rr = abs(bsxfun(@minus, hexCoordsDegs, roi.centerDegs));
    marginDegs = 5/theMosaic.micronsPerDegree;
    
    % Determine distance from ROI center
    selectedConeIndices = find((rr(:,1) <= roi.sizeDegs(1)+marginDegs) & (rr(:,2) <= roi.sizeDegs(2)+marginDegs));
    xLimDegs = (roi.centerDegs(1) + 0.5*roi.sizeDegs(1)*[-1 1]);
    yLimDegs = (roi.centerDegs(2) + 0.5*roi.sizeDegs(2)*[-1 1]);
    
    fprintf('Number of cones in the selected region: %d (total mosaic contains: %d cones)', numel(selectedConeIndices), size(hexCoordsDegs,1));
    
    meshEdgeColor = [0.5 0.5 0.5];
    meshFaceColor = [1.0 0.9 0.9];
    meshFaceAlpha = 0.5;
    meshEdgeAlpha = 0.8;
    lineStyle = '-';
    hFig = figure(1);clf;
    set(hFig, 'Position', [10 10 2800 600]);
    axesHandle  = subplot('Position', [0.02 0.5 0.97 0.45]);
    % First render the mesh
    theMosaic.renderHexMesh(axesHandle, hexCoordsDegs(selectedConeIndices,1), hexCoordsDegs(selectedConeIndices,2), ...
        meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    hold(axesHandle, 'on');
    
    
    lineWidth = 1;
    
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
            hexCoordsMeters(selectedConeIndices,1), hexCoordsMeters(selectedConeIndices,2) ...
        );
    
    innerApertureOutlineVarying.x = innerApertureOutlineVarying.x * metersToDegs;
    innerApertureOutlineVarying.y = innerApertureOutlineVarying.y * metersToDegs;
    
    % Now render the L-cones
    lConeIndices = find(coneTypesHexGrid(selectedConeIndices) == 2);
    if (~isempty(lConeIndices))
        faceColor = [1 0.5 0.5];
        edgeColor = [1 0 0];
        coneApertures.x = innerApertureOutlineVarying.x(lConeIndices,:);
        coneApertures.y = innerApertureOutlineVarying.y(lConeIndices,:);
        theMosaic.renderPatchArray(axesHandle, coneApertures, ...
            (hexCoordsDegs(selectedConeIndices(lConeIndices),1))', hexCoordsDegs(selectedConeIndices(lConeIndices),2), edgeColor, ...
            faceColor, lineStyle, lineWidth);
    end
    
    % Now render the M-cones
    mConeIndices = find(coneTypesHexGrid(selectedConeIndices) == 3);
    if (~isempty(mConeIndices))
        faceColor = [0.5 1 0.5];
        edgeColor = [0 1 0];
        coneApertures.x = innerApertureOutlineVarying.x(mConeIndices,:);
        coneApertures.y = innerApertureOutlineVarying.y(mConeIndices,:);
        theMosaic.renderPatchArray(axesHandle, coneApertures, ...
            (hexCoordsDegs(selectedConeIndices(mConeIndices),1))', hexCoordsDegs(selectedConeIndices(mConeIndices),2), edgeColor, ...
            faceColor, lineStyle, lineWidth);
    end
    
    
    % Now render the S-cones
    sConeIndices = find(coneTypesHexGrid(selectedConeIndices) == 4);
    if (~isempty(sConeIndices))
        faceColor = [0.5 0.5 1];
        edgeColor = [0 0 1];
        coneApertures.x = innerApertureOutlineVarying.x(sConeIndices,:);
        coneApertures.y = innerApertureOutlineVarying.y(sConeIndices,:);
        theMosaic.renderPatchArray(axesHandle, coneApertures, ...
            (hexCoordsDegs(selectedConeIndices(sConeIndices),1))', hexCoordsDegs(selectedConeIndices(sConeIndices),2), edgeColor, ...
            faceColor, lineStyle, lineWidth);
    end
    
    axis(axesHandle, 'equal')
    set(axesHandle, 'XLim', xLimDegs, 'YLim', yLimDegs, 'FontSize', 12);
    xlabel(axesHandle,'spatial position, x (degs)');
    ylabel(axesHandle, 'spatial position, y (degs)');
end

