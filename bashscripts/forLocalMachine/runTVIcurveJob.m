function runTVIcurveJob()

    close all;
    
    % Set 
    nTrainingSamples = 256;
    lowContrast =  3e-4;
    highContrast = 0.3;
    fieldOfViewDegs = 1.0;
    testDiameterDegs = 0.5;
    nPedestalLuminanceLevels = 10;
    
    computeMosaic = false;
    osNoise = 'random';
    computeResponses = ~true;
    ramPercentageEmployed = 1.0;
    pedestalLuminanceListIndicesUsed = 1:10;
    
    visualizeResponses = ~true;
    visualizeOuterSegmentFilters = ~true;
    findPerformance = true;
    visualizePerformance = ~true;
    visualizeSpatialScheme = ~true;
    
    if (computeResponses) || (visualizeResponses) || (visualizeOuterSegmentFilters)
        c_TVIcurve(...
            'pupilDiamMm', 3.0, ...
            'computeMosaic',computeMosaic, ...
            'osNoise', osNoise, ...
            'computeResponses', computeResponses, ...
            'ramPercentageEmployed', ramPercentageEmployed, ...
            'nTrainingSamples', nTrainingSamples, ...
            'lowContrast', lowContrast, ...
            'highContrast', highContrast, ...
            'nPedestalLuminanceLevels', nPedestalLuminanceLevels, ...
            'pedestalLuminanceListIndicesUsed', pedestalLuminanceListIndicesUsed, ...
            'fieldOfViewDegs', fieldOfViewDegs, ...
            'testDiameterDegs', testDiameterDegs, ...
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'visualizeResponses', visualizeResponses, ...
            'visualizeOuterSegmentFilters', visualizeOuterSegmentFilters, ...
            'findPerformance', false, ...
            'visualizePerformance', false ...
            );
    end
    
    if (findPerformance) || (visualizePerformance) || (visualizeSpatialScheme)
        c_TVIcurve( ...
            'pupilDiamMm', 3.0, ...
            'osNoise', osNoise, ...
            'nTrainingSamples', nTrainingSamples, ...
            'lowContrast', lowContrast, ...
            'highContrast', highContrast, ...
            'nPedestalLuminanceLevels', nPedestalLuminanceLevels, ...
            'pedestalLuminanceListIndicesUsed', pedestalLuminanceListIndicesUsed, ...
            'fieldOfViewDegs', fieldOfViewDegs, ...
            'testDiameterDegs', testDiameterDegs, ...
            'computeResponses', false, ...
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'visualizeTransformedSignals', true, ...
            'findPerformance', findPerformance, ...
            'visualizePerformance', visualizePerformance, ...
            'performanceSignal', 'photocurrents', ... % 'isomerizations', ... % 'photocurrents', ...
            'performanceClassifier', 'svmGaussianRF', ... % 'mlpt', 'svmGaussianRF', ...
            'visualizeResponses', false);
    end
    
end


