function pedestalModulationDisk = pedestalModulationDisk(spatialParams, contrast)

    border = min([min([spatialParams.row spatialParams.col])-1 max([0 0.5*(spatialParams.fieldOfViewDegs - spatialParams.pedestalDiameterDegs)/spatialParams.fieldOfViewDegs*spatialParams.row])]);
    i1 = border;
    i2 = spatialParams.row - (border);
    j1 = border;
    j2 = spatialParams.col - (border);
    
    pedestalModulationDisk = zeros(spatialParams.row, spatialParams.col);
    pedestalModulationDisk(i1:i2, j1:j2) = contrast;
end

