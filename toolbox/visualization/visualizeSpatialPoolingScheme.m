function hFig = visualizeSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            spatialPoolingFilter, coneLocsInDegs, mosaicFOVDegs, stimulusFOVDegs)

    zLevels = 0.05:0.2:1.0;
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
    set(hFig, 'Position', [10 10 2050 1230]);
        
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
    plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
    hold off;
    axis 'xy'; axis 'image'
    set(gca, 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02);
    set(gca, 'FontSize', 14, 'XTickLabels', {});
    ylabel('degrees');
    title('stimulus and cone mosaic (cone mosaic view)', 'FontSize', 16);

    if (isfield(spatialPoolingFilter, 'RFprofile'))
        subplot('Position', subplotPosVectors(1,2).v);
        imagesc(xaxis, yaxis, spatialModulation);
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
        set(gca, 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1], 'XTickLabels', {});
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
    quantizationLevels = 1024;
    if (isfield(spatialPoolingFilter, 'RFprofile')) 
        subplot('Position', subplotPosVectors(2,2).v);
        quantizedWeights = round(spatialPoolingFilter.poolingWeights*quantizationLevels);
        title('spatial pooling weights (cone mosaic view)');
        hold on;
        for iLevel = 1:quantizationLevels
            idx = find(quantizedWeights == iLevel);
            c = (iLevel/max(abs(quantizedWeights)))*[1 0.5 0.5];
            plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            idx = find(quantizedWeights == -iLevel);
            c = (iLevel/max(abs(quantizedWeights)))*[0.5 0.5 1.0];
            plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', 2, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
        end
        hold off;
        axis 'xy'; axis 'image'
        set(gca, 'Color', [0 0 0], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
        set(gca, 'FontSize', 14);  
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
            hold on;
            for iLevel = 1:quantizationLevels
                idx = find(quantizedWeights == iLevel);
                c = (iLevel/max(abs(quantizedWeights)))*[1 0.5 0.5];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                idx = find(quantizedWeights == -iLevel);
                c = (iLevel/max(abs(quantizedWeights)))*[0.5 0.5 1.0];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            end
            hold off;
            axis 'xy'; axis 'image'
            set(gca, 'Color', [0 0 0], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
            set(gca, 'FontSize', 14);
        end % for k 
    end
            
    colormap(gray(1024));
    drawnow;
end

