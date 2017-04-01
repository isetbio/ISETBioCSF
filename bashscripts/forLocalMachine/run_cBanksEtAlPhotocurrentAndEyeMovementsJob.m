function run_cBanksEtAlPhotocurrentAndEyeMovementsJob()

 
    % 'originalBanks'; 'defaultIsetbio';  'fullIsetbioNoScones'; 'fullIsetbioWithScones'
    mosaicType = 'fullIsetbioWithScones';
    
    % 'singleExposure'; 'timeSeriesNoPhotocurrents'; 'timeSeriesPhotocurrents'
    temporalAnalysis = 'timeSeriesPhotocurrents';
    
    % 'random'; 'frozen0';
    emPathType = 'frozen0'; %'random';     
    centeredEMPaths = false;
    
    % 'isomerizations', 'photocurrents'
    performanceSignal = 'isomerizations';
    
    % 'mlpt', 'svm', 'svmV1FilterBank'
    performanceClassifier =  svmV1FilterBank';
    
    spatialPoolingKernelParams.type = 'V1QuadraturePair';
    spatialPoolingKernelParams.activationFunction = 'energy'; %'fullWaveRectifier'
    spatialPoolingKernelParams.adjustForConeDensity = true;
    
    % What to do ?
    nTrainingSamples = 1024;
    computationInstance = 1;
    
    if (computationIntance == 0)
        % All conditions in 1 MATLAB session
        ramPercentageEmployed = 1.0;  % use all the RAM
        cyclesPerDegreeExamined =  [10 20 40];
    elseif (computationIntance == 1)
        % First half of the conditions in session 1 of 2 parallel MATLAB sessions
        ramPercentageEmployed = 1.0;  % use all of the RAM
        cyclesPerDegreeExamined =  [10];
    elseif (computationIntance == 2)
        % Second half of the conditions in session 2 of 2 parallel MATLAB sessions
        ramPercentageEmployed = 0.5;  % use 1/2 the RAM
        cyclesPerDegreeExamined =  [20];
    elseif (computationIntance == 3)
        % Second half of the conditions in session 2 of 2 parallel MATLAB sessions
        ramPercentageEmployed = 0.5;  % use 1/2 the RAM
        cyclesPerDegreeExamined =  [40];
    end
    
    
    
    computeMosaic = true;
    computeResponses = true;
    visualizeResponses = true;
    visualizeSpatialScheme = true;
    findPerformance = true;
    visualizePerformance = true;
    
    
    
    switch mosaicType
        case 'originalBanks'
            % 1. Original Banks et al mosaic params
            coneSpacingMicrons = 3.0;
            innerSegmentDiameter = 3.0;    % for a circular sensor
            conePacking = 'hexReg';
            
        case 'defaultIsetbio'
            % 2. Default isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hexReg';
    
        case 'fullIsetbioNoScones'
            % 3. spatially-varying density isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hex';
            LMSRatio = [0.67 0.33 0];
            
        case 'fullIsetbioWithScones'
            % 3. spatially-varying density isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hex';
            LMSRatio = [0.60 0.30 0.10];
    end
    
    
    switch temporalAnalysis 
        case 'singleExposure'
            % For 100 ms simulation
            responseStabilizationMilliseconds = 0;
            responseExtinctionMilliseconds = 100;
            integrationTimeMilliseconds =  100;
            lowContrast = 0.001;
            highContrast = 0.3;
            nContrastsPerDirection =  12;    
            
        case 'timeSeriesNoPhotocurrents'
            % For 5 milliseconds simulation
            responseStabilizationMilliseconds = 10;
            responseExtinctionMilliseconds = 100;
            integrationTimeMilliseconds =  7.5;
            lowContrast = 0.001;
            highContrast = 0.3;
            nContrastsPerDirection =  12;    
    
        case 'timeSeriesPhotocurrents'
            % For 5 milliseconds simulation
            responseStabilizationMilliseconds = 100;
            responseExtinctionMilliseconds = 400;   % use 400 for photocurrents computations
            integrationTimeMilliseconds =  5;
            lowContrast = 0.001;
            highContrast = 1.0;
            nContrastsPerDirection =  18;    
    end
    
    
    
    freezeNoise = ~true;
    luminancesExamined =  [34];
    
    
    if (computeResponses) || (visualizeResponses) || (visualizeSpatialScheme) 
        c_BanksEtAlPhotocurrentAndEyeMovements(...
            'cyclesPerDegree', cyclesPerDegreeExamined, ...
            'luminances', luminancesExamined, ...
            'nTrainingSamples', nTrainingSamples, ...
            'lowContrast', lowContrast, ...
            'highContrast', highContrast, ...
            'nContrastsPerDirection', nContrastsPerDirection, ...
            'ramPercentageEmployed', ramPercentageEmployed, ...
            'emPathType', emPathType, ...
            'centeredEMPaths', centeredEMPaths, ...
            'responseStabilizationMilliseconds', responseStabilizationMilliseconds, ...
            'responseExtinctionMilliseconds', responseExtinctionMilliseconds, ...
            'freezeNoise', freezeNoise, ...
            'integrationTime', integrationTimeMilliseconds/1000, ...
            'coneSpacingMicrons', coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter), ...
            'conePacking', conePacking, ...
            'LMSRatio', LMSRatio, ...
            'computeMosaic',computeMosaic, ...
            'computeResponses', computeResponses, ...
            'visualizeResponses', visualizeResponses, ...
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'findPerformance', false, ...
            'visualizePerformance', false, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier ...
        );
    end
    
    if (findPerformance) || (visualizePerformance) 
        c_BanksEtAlPhotocurrentAndEyeMovements(...
            'cyclesPerDegree', cyclesPerDegreeExamined, ...
            'luminances', luminancesExamined, ...
            'nTrainingSamples', nTrainingSamples, ...
            'nContrastsPerDirection', nContrastsPerDirection, ...
            'lowContrast', lowContrast, ...
            'highContrast', highContrast, ...
            'ramPercentageEmployed', ramPercentageEmployed, ...
            'emPathType', emPathType, ...
            'centeredEMPaths', centeredEMPaths, ...
            'responseStabilizationMilliseconds', responseStabilizationMilliseconds, ...
            'responseExtinctionMilliseconds', responseExtinctionMilliseconds, ...
            'freezeNoise', freezeNoise, ...
            'integrationTime', integrationTimeMilliseconds/1000, ...
            'coneSpacingMicrons', coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter), ...
            'conePacking', conePacking, ...
            'LMSRatio', LMSRatio, ...
            'computeMosaic', false, ...
            'computeResponses', false, ...
            'visualizeResponses', false, ...
            'visualizeSpatialScheme', visualizeSpatialScheme, ...
            'findPerformance', findPerformance, ...
            'visualizePerformance', visualizePerformance, ...
            'visualizeTransformedSignals', true, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier, ...
            'spatialPoolingKernelParams', spatialPoolingKernelParams ...
        );
    end
    
    
end

