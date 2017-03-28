function hFig = visualizeSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            spatialPoolingKernelParams, spatialPoolingFilter, coneLocsInDegs, mosaicFOVDegs, stimulusFOVDegs)
        
        
    zLevels = [0.025:0.05:1.0];
    zLevels = [-fliplr(zLevels) zLevels];
        
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 1, ...
           'heightMargin',   0.00, ...
           'widthMargin',    0.00, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.01, ...
           'topMargin',      0.005);
       
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 850 740], 'Color', [1 1 1]);
    
    coneRadius = 1.0/300;
    coneX = coneRadius * cos(2*pi*(0:30:360)/360);
    coneY = coneRadius * sin(2*pi*(0:30:360)/360);
    coneMarkerSize = 7;
    quantizationLevels = 1024;
    if (strcmp(spatialPoolingKernelParams.type, 'GaussianRF'))
        subplot('Position', subplotPosVectors(1,1).v);
        quantizedWeights = round(spatialPoolingFilter.poolingWeights*quantizationLevels);
        quantizedWeights(quantizedWeights > quantizationLevels) = quantizationLevels;
        title('spatial pooling weights (cone mosaic view)');
        imagesc(xaxis, yaxis, 0.5 + 0.2*spatialModulation);
        hold on;
        for k = 1:size(coneLocsInDegs,1)
            fill(squeeze(coneLocsInDegs(k,1))+coneX, squeeze(coneLocsInDegs(k,2))+coneY, [0.5 0.5 0.5]);
        end
        for iLevel = 1:quantizationLevels
            idx = find(quantizedWeights == iLevel);
            if (~isempty(idx))
            for k = 1:numel(idx)
            c = [0.5 0.5 0.5] + (iLevel/max(abs(quantizedWeights)))*[0.5 -0.4 -0.4];
            fill(squeeze(coneLocsInDegs(idx(k),1))+coneX, squeeze(coneLocsInDegs(idx(k),2))+coneY,  c);
            end
            end
            idx = find(quantizedWeights == -iLevel);
            if (~isempty(idx))
                for k = 1:numel(idx)
            c = [0.5 0.5 0.5] + (iLevel/max(abs(quantizedWeights)))*[-0.4 -0.4 0.5];
            fill(squeeze(coneLocsInDegs(idx(k),1))+coneX, squeeze(coneLocsInDegs(idx(k),2))+coneY,  c);
                end
            end
        end

        hold off;
        axis 'xy'; axis 'image'
        set(gca, 'CLim', [0 1], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
        set(gca, 'FontSize', 20);  
    else
        error('UNknown spatialPooling filter: %s\n', spatialPoolingKernelParams.type);
    end
      
    colormap(gray(quantizationLevels));
    drawnow;
end

function renderPatchArray(pixelOutline, xCoords, yCoords, faceColorsNormalizedValues,  edgeColor, lineWidth)
    verticesPerCone = numel(pixelOutline.x);
    verticesList = zeros(verticesPerCone * numel(xCoords), 2);
    facesList = [];
    colors = [];
    for coneIdx = 1:numel(xCoords)
        idx = (coneIdx-1)*verticesPerCone + (1:verticesPerCone);
        verticesList(idx,1) = pixelOutline.x(:) + xCoords(coneIdx);
        verticesList(idx,2) = pixelOutline.y(:) + yCoords(coneIdx);
        facesList = cat(1, facesList, idx);
        colors = cat(1, colors, repmat(faceColorsNormalizedValues(coneIdx), [verticesPerCone 1]));
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.LineWidth = lineWidth;
    patch(S);
end


function hFig = visualizeSpatialPoolingSchemeOLD(xaxis, yaxis, spatialModulation, ...
            spatialPoolingFilter, coneLocsInDegs, mosaicFOVDegs, stimulusFOVDegs)

    zLevels = [0.05:0.05:1.0]; % 0.05:0.2:1.0;
    zLevels = [-fliplr(zLevels) zLevels];
        
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.015, ...
           'topMargin',      0.005);
       
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2050 1230], 'Color', [1 1 1]);
        
    subplot('Position', subplotPosVectors(1,1).v);
    imagesc(xaxis, yaxis, spatialModulation);
    hold on;
    % Spatial RF in cyan
    contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, zLevels(zLevels>0), 'Color', 'c', 'LineWidth', 1.0, 'LineStyle', '-');
    % outline mosaic extent in green
    x = mosaicFOVDegs * [-0.5 0.5 0.5 -0.5 -0.5];
    y = mosaicFOVDegs * [-0.5 -0.5 0.5 0.5 -0.5];
       
    plot(x,y, 'g-', 'LineWidth', 1.5);
    hold off
    axis 'xy'; axis 'image'
    set(gca, 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]);
    set(gca, 'FontSize', 14, 'XTickLabels', {});
    ylabel('degrees');
    title('stimulus, cone mosaic, and spatial pooling profile (stimulus view)');
         
    subplot('Position', subplotPosVectors(2,1).v);
    imagesc(xaxis, yaxis, spatialModulation);
    hold on;
    plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'ro', 'MarkerSize', 6);
    hold off;
    axis 'xy'; axis 'image'
    set(gca, 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02);
    set(gca, 'FontSize', 14, 'XTickLabels', {});
    ylabel('degrees');
    title('stimulus and cone mosaic (cone mosaic view)', 'FontSize', 16);

    if (isfield(spatialPoolingFilter, 'RFprofile'))
        subplot('Position', subplotPosVectors(1,2).v);
        imagesc(xaxis, yaxis, 0.5 + 0.5*spatialModulation);
        hold on;
        % outline mosaic extent in green
        x = mosaicFOVDegs * [-0.5 0.5 0.5 -0.5 -0.5];
        y = mosaicFOVDegs * [-0.5 -0.5 0.5 0.5 -0.5];
        plot(x,y, 'g-', 'LineWidth', 1.5);
        % Gaussian pooling
        contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, zLevels(zLevels>0), 'Color', 'r', 'LineWidth', 1.0, 'LineStyle', '-');
        contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, zLevels(zLevels<0), 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '-');
        plot(x,y, 'g-', 'LineWidth', 1.5);
        hold off
        axis 'xy'; axis 'image'
        set(gca, 'CLim', [0 1], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'XTickLabels', {});
        set(gca, 'FontSize', 16);
        title('stimulus, cone mosaic and spatial pooling profile (stimulus view)');
        
    elseif (isfield(spatialPoolingFilter, 'cosPhasePoolingProfile'))
        
        for k = 2:3
            subplot('Position', subplotPosVectors(1,k).v);
            imagesc(xaxis, yaxis, spatialModulation);
            axis 'xy';
            axis 'image';
            hold on;
            % outline mosaic extent in green
            x = mosaicFOVDegs * [-0.5 0.5 0.5 -0.5 -0.5];
            y = mosaicFOVDegs * [-0.5 -0.5 0.5 0.5 -0.5];
            plot(x,y, 'g-', 'LineWidth', 1.5);
            if (k == 2)
                % V1 cos-phase filter
                contour(xaxis, yaxis, spatialPoolingFilter.cosPhasePoolingProfile, zLevels(zLevels>0), 'Color', 'r', 'LineWidth', 1.0, 'LineStyle', '-');
                contour(xaxis, yaxis, spatialPoolingFilter.cosPhasePoolingProfile, zLevels(zLevels<0), 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '-');
                title('stimulus, cone mosaic and cos-phase V1 filter profile (stimulus view)');
            else
                % V1 sin-phase filter
                contour(xaxis, yaxis, spatialPoolingFilter.sinPhasePoolingProfile, zLevels(zLevels>0), 'Color', 'r', 'LineWidth', 1.0, 'LineStyle', '-');
                contour(xaxis, yaxis, spatialPoolingFilter.sinPhasePoolingProfile, zLevels(zLevels<0), 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '-');
                title('stimulus, cone mosaic and sin-phase V1 filter profile (stimulus view)');
            end
            plot(x,y, 'g-', 'LineWidth', 1.5);
            hold off
            set(gca, 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'XTickLabels', {});
            set(gca, 'FontSize', 16);
        end % k
    end
    
    % Show cone pooling weights
    coneMarkerSize = 10;
    quantizationLevels = 1024;
    if (isfield(spatialPoolingFilter, 'RFprofile')) 
        subplot('Position', subplotPosVectors(2,2).v);
        quantizedWeights = round(spatialPoolingFilter.poolingWeights*quantizationLevels);
        quantizedWeights(quantizedWeights > quantizationLevels) = quantizationLevels;
        title('spatial pooling weights (cone mosaic view)');
        imagesc(xaxis, yaxis, 0.5 + 0.2*spatialModulation);
        hold on;
        plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'ro', 'MarkerSize', coneMarkerSize, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3]);
        for iLevel = 1:quantizationLevels
            idx = find(quantizedWeights == iLevel);
            c = [0.5 0.5 0.5] + (iLevel/max(abs(quantizedWeights)))*[0.5 -0.5 -0.5];
            plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', coneMarkerSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', [0.3 0.3 0.3]);
            idx = find(quantizedWeights == -iLevel);
            c = [0.5 0.5 0.5] + (iLevel/max(abs(quantizedWeights)))*[-0.5 -0.5 0.5];
            plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', coneMarkerSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', [0.3 0.3 0.3]);
        end
        hold off;
        axis 'xy'; axis 'image'
        set(gca, 'CLim', [0 1], 'Color', [0 0 0], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
        set(gca, 'FontSize', 20);  
    else
        for k = 2:3
            subplot('Position', subplotPosVectors(2,k).v); 
            if (k == 2)
                quantizedWeights = round(V1filterBank.cosPhasePoolingWeights*quantizationLevels);
                title('cos-phase V1 filter pooling weights (cone mosaic view)');
            else
                quantizedWeights = round(V1filterBank.sinPhasePoolingWeights*quantizationLevels);
                title('sin-phase V1 filter pooling weights (cone mosaic view)');
            end
            imagesc(xaxis, yaxis, 0.5 + 0.2*spatialModulation);
            hold on;
            plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'ro', 'MarkerSize', coneMarkerSize, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3]);
            for iLevel = 1:quantizationLevels
                idx = find(quantizedWeights == iLevel);
                c = [0.5 0.5 0.5] + (iLevel/max(abs(quantizedWeights)))*[0.5 -0.5 -0.5];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', coneMarkerSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', [0.3 0.3 0.3]);
                idx = find(quantizedWeights == -iLevel);
                c = [0.5 0.5 0.5] + (iLevel/max(abs(quantizedWeights)))*[-0.5 -0.5 0.5];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', coneMarkerSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', [0.3 0.3 0.3]);
            end
            hold off;
            axis 'xy'; axis 'image'
            set(gca, 'CLim', [0 1], 'Color', [0 0 0], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
            set(gca, 'FontSize', 20);
        end % for k 
    end
            
    colormap(gray(1024));
    drawnow;
end

