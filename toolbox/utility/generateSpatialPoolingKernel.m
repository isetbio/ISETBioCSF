function [spatialPoolingFilter, hFig] = generateSpatialPoolingKernel(spatialParams, mosaicParams, topLevelDirParams, visualizeSpatialScheme, thresholdParams, paramsList)

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
    
    
    if (strcmp(mosaicParams.conePacking, 'hex'))
        % Find the density map around each cone
        eccInMeters = sqrt(sum(coneLocsInMeters.^2, 2));
        ang = atan2(squeeze(coneLocsInMeters(:,2)), squeeze(coneLocsInMeters(:,1)))/pi*180;
        [~, ~, coneDensity] = coneSizeReadData('eccentricity',eccInMeters(:),'angle',ang(:));
    else
        coneDensity = 1;
    end
    
    switch(spatialParams.spatialType)
        case 'pedestalDisk' 
            spatialModulation = pedestalModulationDisk(spatialParams, 1.0);
        case 'spot' 
            spatialModulation = drawSpot(spatialParams);
        otherwise
            error('Unknown spatial type specified: ''%s''.', spatialParams.spatialType);
    end % switch
    
    xaxis = (0:(size(spatialModulation,2)-1))/size(spatialModulation,2) * spatialParams.fieldOfViewDegs;
    xaxis = xaxis - mean(xaxis);
    yaxis = (0:(size(spatialModulation,1)-1))/size(spatialModulation,1) * spatialParams.fieldOfViewDegs;
    yaxis = yaxis - mean(yaxis);
    
    % Generate the spatial pooling filter      
    spatialPoolingFilter = makeSpatialPoolingFilter(spatialParams, coneLocsInDegs, xaxis, yaxis, coneDensity, thresholdParams.spatialPoolingKernelParams);
    
    if (visualizeSpatialScheme)
        hFig = visualizeSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            thresholdParams.spatialPoolingKernelParams, spatialPoolingFilter, coneLocsInDegs, mosaicParams.fieldOfViewDegs, spatialParams.fieldOfViewDegs);
    
        % Save figure
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('SpatialPooling');
        paramsList{numel(paramsList)+1} = thresholdParams;
        rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
       
    end % visualizeSpatialScheme

end

function spatialPoolingFilter = makeSpatialPoolingFilter(spatialParams, coneLocsDegs, xaxis, yaxis, coneDensity, spatialPoolingKernelParams)

    spatialPoolingFilter = [];
    [X,Y] = meshgrid(xaxis, yaxis);
    
    switch spatialParams.spatialType
        case 'Gabor'
            testDiameterDegs = spatialParams.testDiameterDegs;
        case 'spot'
            testDiameterDegs = spatialParams.spotSizeDegs;
        case 'pedestalDisk'
            testDiameterDegs = spatialParams.testDiameterDegs;
        otherwise
            error('Do not know how to get testDiameter from ''%s''.', spatialParams.spatialType)
    end
    
    switch spatialPoolingKernelParams.type
        case 'GaussianRF'
            if (spatialPoolingKernelParams.shrinkageFactor > 0)
                sigma = spatialPoolingKernelParams.shrinkageFactor * testDiameterDegs/2.0;
            else
                sigma = spatialPoolingKernelParams.shrinkageFactor;
            end
            RFprofile = exp(-0.5*(X/sigma).^2) .* exp(-0.5*(Y/sigma).^2);
        otherwise
            error('Unknown spatialPoolingKernelType: ''%s''.', spatialPoolingKernelParams.type);
    end % switch
    
    spatialPoolingFilter.RFprofile = RFprofile / max(abs(RFprofile(:)));
    
    rfCoordsDegs = [X(:) Y(:)];
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
    spatialPoolingFilter.poolingWeights = spatialPoolingFilter.RFprofile(idx);
    maxWeight = max(spatialPoolingFilter.poolingWeights(:));
     
    % Adjust weights by the inverse of the coneDensity
    if (spatialPoolingKernelParams.adjustForConeDensity)
        spatialPoolingFilter.poolingWeights = spatialPoolingFilter.poolingWeights ./ coneDensity;
    end
    spatialPoolingFilter.poolingWeights = spatialPoolingFilter.poolingWeights / max(abs(spatialPoolingFilter.poolingWeights(:))) * maxWeight;
end

