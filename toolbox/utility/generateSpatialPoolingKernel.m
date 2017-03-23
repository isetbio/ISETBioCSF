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
    
    switch(spatialParams.spatialType)
        case 'pedestalDisk' 
            spatialModulation = pedestalModulationDisk(spatialParams, 1.0);
        otherwise
            error('Unknown spatial type specified: ''%s''.', spatialParams.spatialType);
    end % switch
    
    xaxis = (0:(size(spatialModulation,2)-1))/size(spatialModulation,2) * spatialParams.fieldOfViewDegs;
    xaxis = xaxis - mean(xaxis);
    yaxis = (0:(size(spatialModulation,1)-1))/size(spatialModulation,1) * spatialParams.fieldOfViewDegs;
    yaxis = yaxis - mean(yaxis);
    
    % Generate the spatial pooling filter      
    spatialPoolingFilter = makeSpatialPoolingFilter(spatialParams, coneLocsInDegs, xaxis, yaxis, coneDensity, spatialPoolingKernelType);
    
    if (visualizeSpatialScheme)
        hFig = visualizeSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            spatialPoolingFilter, coneLocsInDegs, mosaicParams.fieldOfViewDegs, spatialParams.fieldOfViewDegs);
    end % visualizeSpatialScheme

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

