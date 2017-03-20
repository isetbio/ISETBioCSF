function spatialPoolingFilter = generateSpatialPoolingKernel(spatialParams, mosaicParams, topLevelDirParams, visualizeSpatialScheme, spatialPoolingKernelType)

    % Load the mosaic
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    
    % Get the cone locs in degrees
    if (strcmp(mosaicParams.conePacking,'rect'))
        coneLocsInMeters = theMosaic.coneLocs;
    else
        coneLocsInMeters = theMosaic.coneLocsHexGrid;
    end
    coneLocsInDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsInDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    
    % Find the density map around each cone
    eccInMeters = sqrt(sum(coneLocsInMeters.^2, 2));
    ang = atan2(squeeze(coneLocsInMeters(:,2)), squeeze(coneLocsInMeters(:,1)))/pi*180;
    [~, ~, coneDensity] = coneSize(eccInMeters(:),ang(:));
    
    if (visualizeSpatialScheme)
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
    end % visualizeSpatialScheme
    
    switch(spatialParams.spatialType)
        case 'pedestalDisk' 
            spatialModulation = pedestalModulationDisk(spatialParams, 1.0);
        otherwise
            error('Unknown spatial type specified: ''%s''.', spatialParams.spatialType);
    end
    
    xaxis = (0:(size(spatialModulation,2)-1))/size(spatialModulation,2) * spatialParams.fieldOfViewDegs;
    xaxis = xaxis - mean(xaxis);
    yaxis = (0:(size(spatialModulation,1)-1))/size(spatialModulation,1) * spatialParams.fieldOfViewDegs;
    yaxis = yaxis - mean(yaxis);
    
    % Generate the spatial pooling filter      
    spatialPoolingFilter = makeSpatialPoolingFilter(spatialParams, coneLocsInDegs, xaxis, yaxis, coneDensity, spatialPoolingKernelType);
    
    if (visualizeSpatialScheme)
        zLevels = 0.05:0.2:1.0;
        zLevels = [-fliplr(zLevels) zLevels];
        
        subplot('Position', subplotPosVectors(1,1).v);
            imagesc(xaxis, yaxis, spatialModulation);
            hold on;
            % Spatial RF in cyan
            contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, zLevels(zLevels>0), 'Color', 'c', 'LineWidth', 1.0, 'LineStyle', '-');
            % outline mosaic extent in green
            x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
            y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
            plot(x,y, 'g-', 'LineWidth', 1.5);
            hold off
            axis 'xy'; axis 'image'
            set(gca, 'XLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1], 'YLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1]);
            set(gca, 'FontSize', 14, 'XTickLabels', {});
            ylabel('degrees');
            title('stimulus, cone mosaic, and spatial pooling profile (stimulus view)');
         
        subplot('Position', subplotPosVectors(2,1).v);
            imagesc(xaxis, yaxis, spatialModulation);
            hold on;
            plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
            hold off;
            axis 'xy'; axis 'image'
            set(gca, 'XLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1]*1.02);
            set(gca, 'FontSize', 14, 'XTickLabels', {});
            ylabel('degrees');
            title('stimulus and cone mosaic (cone mosaic view)', 'FontSize', 16);

        subplot('Position', subplotPosVectors(1,2).v);
            imagesc(xaxis, yaxis, spatialModulation);
            hold on;
            % outline mosaic extent in green
            x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
            y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
            plot(x,y, 'g-', 'LineWidth', 1.5);
            % V1 cos-phase filter
            contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, zLevels(zLevels>0), 'Color', 'r', 'LineWidth', 1.0, 'LineStyle', '-');
            contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, zLevels(zLevels<0), 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '-');
            plot(x,y, 'g-', 'LineWidth', 1.5);
            hold off
            axis 'xy'; axis 'image'
            set(gca, 'XLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1], 'YLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1], 'XTickLabels', {});
            set(gca, 'FontSize', 16);
            title('stimulus, cone mosaic and spatial pooling profile (stimulus view)');
            
        quantizationLevels = 1024;
        subplot('Position', subplotPosVectors(2,2).v); 
            quantizedWeights = round(spatialPoolingFilter.poolingWeights*quantizationLevels);
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
            set(gca, 'Color', [0 0 0], 'XLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicParams.fieldOfViewDegs spatialParams.fieldOfViewDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
            set(gca, 'FontSize', 14);
            title('spatial pooling weights (cone mosaic view)');
            
        colormap(gray(1024));
        drawnow;
    end % (visualizeSpatialScheme)

end

function spatialPoolingFilter = makeSpatialPoolingFilter(spatialParams, coneLocsDegs, xaxis, yaxis, coneDensity, spatialPoolingKernelType)

    spatialPoolingFilter = [];
  
    [X,Y] = meshgrid(xaxis, yaxis);
    switch spatialPoolingKernelType
        case 'svmGaussianRF'
            sigma = spatialParams.testDiameterDegs/1.5;
            RFprofile = exp(-0.5*(X/sigma).^2) .* exp(-0.5*(Y/sigma).^2);
        otherwise
            error('Unknown spatialPoolingKernelType: ''%s''.', spatialPoolingKernelType);
    end % switch
    
    spatialPoolingFilter.RFprofile = RFprofile / max(abs(RFprofile(:)));
    
    rfCoordsDegs = [X(:) Y(:)];
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
    spatialPoolingFilter.poolingWeights = spatialPoolingFilter.RFprofile(idx);
    
    % Adjust weights by the inverse of the coneDensity
    maxWeight = max(spatialPoolingFilter.poolingWeights(:));
    spatialPoolingFilter.poolingWeights = spatialPoolingFilter.poolingWeights ./ coneDensity;
    spatialPoolingFilter.poolingWeights = spatialPoolingFilter.poolingWeights / max(abs(spatialPoolingFilter.poolingWeights(:))) * maxWeight;
    
    % Normalize with respect to total energy over space
    netWeight = 1/sqrt(2.0) * sqrt(sum(spatialPoolingFilter.poolingWeights(:).^2));
    spatialPoolingFilter.poolingWeights = spatialPoolingFilter.poolingWeights / netWeight;
end

