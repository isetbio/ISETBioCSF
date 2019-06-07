function spatialPoolingKernels = generateSpatialPoolingKernels(theMosaic, ...
            lowFrequencyTemplate, lowFrequencyTemplateOrtho, ...
            highFrequencyTemplate, highFrequencyTemplateOrtho, spatialSupportDegs)
    
    coneLocsMeters = theMosaic.coneLocsHexGrid;
    coneLocsDegs(:,1) = coneLocsMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsDegs(:,2) = coneLocsMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    
    [X,Y] = meshgrid(spatialSupportDegs, spatialSupportDegs);
    rfCoordsDegs = [X(:) Y(:)];
    
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
    
    size(lowFrequencyTemplate.linear)
    figure(234234);
    subplot(2,2,1)
    imagesc(lowFrequencyTemplate.linear)
    subplot(2,2,2);
    imagesc(lowFrequencyTemplateOrtho.linear);
    subplot(2,2,3)
    imagesc(lowFrequencyTemplate.quadrature)
    subplot(2,2,4);
    imagesc(lowFrequencyTemplateOrtho.quadrature);
    pause

    max(abs(lowFrequencyTemplate.quadrature(:)))
    max(abs(lowFrequencyTemplateOrtho.linear(:)))
    pause
    
    spatialPoolingKernels.lowFrequencyPoolingWeightsLinear           = lowFrequencyTemplate.linear(idx);
    spatialPoolingKernels.lowFrequencyPoolingWeightsQuadrature       = lowFrequencyTemplate.quadrature(idx);

    spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsLinear      = lowFrequencyTemplateOrtho.linear(idx);
    spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsQuadrature  = lowFrequencyTemplateOrtho.quadrature(idx);
    
    spatialPoolingKernels.highFrequencyPoolingWeightsLinear          = highFrequencyTemplate.linear(idx);
    spatialPoolingKernels.highFrequencyPoolingWeightsQuadrature      = highFrequencyTemplate.quadrature(idx);

    spatialPoolingKernels.highFrequencyOrthoPoolingWeightsLinear     = highFrequencyTemplateOrtho.linear(idx);
    spatialPoolingKernels.highFrequencyOrthoPoolingWeightsQuadrature = highFrequencyTemplateOrtho.quadrature(idx);
    
    coneRadiusMicrons = 1.4;
    coneRadiusDegs = coneRadiusMicrons/theMosaic.micronsPerDegree;
    spatialPoolingKernels.coneAperture(:,1) = coneRadiusDegs * cosd(0:60:360);
    spatialPoolingKernels.coneAperture(:,2) = coneRadiusDegs * sind(0:60:360);
    spatialPoolingKernels.coneLocsDegs = coneLocsDegs;

end