function pedestalModulationDisk = pedestalModulationDisk(spatialParams, contrast)

    x = (0:(spatialParams.col-1))/((spatialParams.col-1))*spatialParams.fieldOfViewDegs;
    x = x - mean(x); y = x;
    [X,Y] = meshgrid(x,y);
   
    % Make a hard edge disk
    pedestalModulationDisk = zeros(size(X));
    pedestalModulationDisk(sqrt(X.^2 + Y.^2) <= spatialParams.testDiameterDegs/2) = 1;
    
    % Make the edges softer
    sigma = round(spatialParams.col/50);
    pedestalModulationDisk = conv2(pedestalModulationDisk, ones(sigma,sigma), 'same');
    
    % Adjust to peak contrast
    pedestalModulationDisk = pedestalModulationDisk / max(pedestalModulationDisk(:)) * contrast;
end


    
    