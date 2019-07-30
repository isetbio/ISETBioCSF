function varargout = run_BanksPhotocurrentEyeMovementConditions(params)
 
varargout = {};
varargout{1} = [];  % the mosaic
varargout{2} = [];  % psychometric functions
varargout{3} = [];  % the figData

    % Assign parfor workers num based on computer name
    [~, computerNetworkName] = system('uname -n');
    if (contains(computerNetworkName, 'ionean'))
        params.parforWorkersNum = 10;
    elseif (contains(computerNetworkName, 'leviathan'))
        params.parforWorkersNum = 20;
    else
        params.parforWorkersNum = 4;
    end
    fprintf('\nWill use %d workers\n', params.parforWorkersNum);

    if (params.deleteResponseInstances)
        c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', params.opticsModel, ...
            'wavefrontSpatialSamples', params.wavefrontSpatialSamples, ...              
            'opticalImagePadSizeDegs', params.opticalImagePadSizeDegs, ...     
            'pupilDiamMm', params.pupilDiamMm, ...
            'blur', params.blur, ...
            'apertureBlur', params.apertureBlur, ...
            'imagePixels', params.imagePixels, ...
            'cyclesPerDegree', params.cyclesPerDegreeExamined, ...
            'stimulusOrientationDegs', params.stimulusOrientationDegs, ...
            'luminances', params.luminancesExamined, ...
            'coneContrastDirection', params.coneContrastDirection, ...
            'nTrainingSamples', params.nTrainingSamples, ...
            'lowContrast', params.lowContrast, ...
            'highContrast', params.highContrast, ...
            'nContrastsPerDirection', params.nContrastsPerDirection, ...
            'frameRate', params.frameRate, ...
            'ramPercentageEmployed', params.ramPercentageEmployed, ...
            'emPathType', params.emPathType, ...
            'centeredEMPaths', params.centeredEMPaths, ...
            'stimulusDurationInSeconds', params.stimulusDurationInSeconds, ...
            'responseStabilizationMilliseconds', params.responseStabilizationMilliseconds, ...
            'responseExtinctionMilliseconds', params.responseExtinctionMilliseconds, ...
            'secondsToInclude', params.secondsToInclude, ...
            'secondsToIncludeOffset', params.secondsToIncludeOffset, ...
            'singleBinTemporalWindowIfPossible', params.singleBinTemporalWindowIfPossible, ...
            'freezeNoise', params.freezeNoise, ...
            'minimumMosaicFOVdegs', params.minimumMosaicFOVdegs, ...
            'integrationTime', params.integrationTimeMilliseconds/1000, ...
            'coneSpacingMicrons', params.coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', params.innerSegmentDiameter, ...
            'conePacking', params.conePacking, ...
            'eccBasedConeQuantalEfficiency', params.eccBasedConeQuantalEfficiency, ...
            'eccBasedMacularPigment', params.eccBasedMacularPigment, ...
            'LMSRatio', params.LMSRatio, ...
            'mosaicRotationDegs', params.mosaicRotationDegs, ...
            'sConeMinDistanceFactor', params.sConeMinDistanceFactor, ...
            'sConeFreeRadiusMicrons', params.sConeFreeRadiusMicrons, ...
            'latticeAdjustmentPositionalToleranceF', params.latticeAdjustmentPositionalToleranceF, ...
            'latticeAdjustmentDelaunayToleranceF', params.latticeAdjustmentDelaunayToleranceF, ...
            'maxGridAdjustmentIterations', params.maxGridAdjustmentIterations, ...   
            'marginF', params.marginF, ...
            'resamplingFactor', params.resamplingFactor, ...                           
            'computeMosaic', false, ...
            'visualizeMosaic', false, ...
            'computeResponses', false, ...
            'computePhotocurrentResponseInstances', false, ...
            'visualizeResponses', false, ...
            'visualizeOuterSegmentFilters', false, ...
            'visualizeSpatialScheme', false, ...
            'findPerformance', false, ...
            'visualizePerformance', false, ...
            'deleteResponseInstances', true);
        return;
    end
    
    if (params.computeResponses) || (params.visualizeResponsesWithSpatialPoolingSchemeInVideo) || (params.visualizeMosaicWithFirstEMpath) || (params.visualizeResponses) || (params.visualizeMosaic) || (params.computeMosaic)

        [~, ~,~,~, theMosaic] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', params.opticsModel, ...
            'wavefrontSpatialSamples', params.wavefrontSpatialSamples, ...               
            'opticalImagePadSizeDegs', params.opticalImagePadSizeDegs, ...  
            'pupilDiamMm', params.pupilDiamMm, ...
            'blur', params.blur, ...
            'apertureBlur', params.apertureBlur, ...
            'imagePixels', params.imagePixels, ...
            'cyclesPerDegree', params.cyclesPerDegreeExamined, ...
            'stimulusOrientationDegs', params.stimulusOrientationDegs, ...
            'luminances', params.luminancesExamined, ...
            'coneContrastDirection', params.coneContrastDirection, ...
            'nTrainingSamples', params.nTrainingSamples, ...
            'lowContrast', params.lowContrast, ...
            'highContrast', params.highContrast, ...
            'nContrastsPerDirection', params.nContrastsPerDirection, ...
            'frameRate', params.frameRate, ...
            'ramPercentageEmployed', params.ramPercentageEmployed, ...
            'parforWorkersNum', params.parforWorkersNum, ...
            'emPathType', params.emPathType, ...
            'centeredEMPaths', params.centeredEMPaths, ...
            'stimulusDurationInSeconds', params.stimulusDurationInSeconds, ...
            'responseStabilizationMilliseconds', params.responseStabilizationMilliseconds, ...
            'responseExtinctionMilliseconds', params.responseExtinctionMilliseconds, ...
            'secondsToInclude', params.secondsToInclude, ...
            'secondsToIncludeOffset', params.secondsToIncludeOffset, ...
            'singleBinTemporalWindowIfPossible', params.singleBinTemporalWindowIfPossible, ...
            'freezeNoise', params.freezeNoise, ...
            'minimumMosaicFOVdegs', params.minimumMosaicFOVdegs, ...
            'integrationTime', params.integrationTimeMilliseconds/1000, ...
            'coneSpacingMicrons', params.coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', params.innerSegmentDiameter, ...
            'conePacking', params.conePacking, ...
            'eccBasedConeQuantalEfficiency', params.eccBasedConeQuantalEfficiency, ...
            'eccBasedMacularPigment', params.eccBasedMacularPigment, ...
            'LMSRatio', params.LMSRatio, ...
            'mosaicRotationDegs', params.mosaicRotationDegs, ...
            'sConeMinDistanceFactor', params.sConeMinDistanceFactor, ...
            'sConeFreeRadiusMicrons', params.sConeFreeRadiusMicrons, ...
            'latticeAdjustmentPositionalToleranceF', params.latticeAdjustmentPositionalToleranceF, ...
            'latticeAdjustmentDelaunayToleranceF', params.latticeAdjustmentDelaunayToleranceF, ...
            'maxGridAdjustmentIterations', params.maxGridAdjustmentIterations, ...   
            'marginF', params.marginF, ...
            'resamplingFactor', params.resamplingFactor, ...                             
            'computeMosaic', params.computeMosaic, ...
            'visualizeMosaic', params.visualizeMosaic, ...
            'computeResponses', params.computeResponses, ...
            'computePhotocurrentResponseInstances', params.computePhotocurrentResponseInstances, ...
            'visualizeDisplay', params.visualizeDisplay, ...
            'visualizeOptics', params.visualizeOptics, ...
            'visualizeStimulusAndOpticalImage', params.visualizeStimulusAndOpticalImage, ...
            'visualizeOIsequence', params.visualizeOIsequence, ...
            'visualizeResponses', params.visualizeResponses, ...
            'visualizeResponsesWithSpatialPoolingSchemeInVideo', params.visualizeResponsesWithSpatialPoolingSchemeInVideo, ...
            'visualizeOuterSegmentFilters', params.visualizeOuterSegmentFilters, ...
            'visualizeMosaicWithFirstEMpath', params.visualizeMosaicWithFirstEMpath, ...
            'visualizeSpatialScheme', params.visualizeSpatialScheme, ...
            'findPerformance', false, ...
            'visualizePerformance', false, ...
            'performanceSignal' , params.performanceSignal, ...
            'performanceClassifier', params.performanceClassifier, ...
            'spatialPoolingKernelParams', params.spatialPoolingKernelParams ...
        );
        varargout{1} = theMosaic;
    end
    
    if (params.findPerformance) || (params.visualizePerformance)
            [perfData, ~, ~, ~, theMosaic, thePsychometricFunctions, theFigData] = c_BanksEtAlPhotocurrentAndEyeMovements(...
                'opticsModel', params.opticsModel, ...
                'wavefrontSpatialSamples', params.wavefrontSpatialSamples, ...               
                'opticalImagePadSizeDegs', params.opticalImagePadSizeDegs, ...      
                'pupilDiamMm', params.pupilDiamMm, ...
                'blur', params.blur, ...
                'apertureBlur', params.apertureBlur, ...
                'imagePixels', params.imagePixels, ...
                'cyclesPerDegree', params.cyclesPerDegreeExamined, ...
                'stimulusOrientationDegs', params.stimulusOrientationDegs, ...
                'luminances', params.luminancesExamined, ...
                'coneContrastDirection', params.coneContrastDirection, ...
                'nTrainingSamples', params.nTrainingSamples, ... 
                'nContrastsPerDirection', params.nContrastsPerDirection, ...
                'frameRate', params.frameRate, ...
                'lowContrast', params.lowContrast, ...
                'highContrast', params.highContrast, ...
                'ramPercentageEmployed', params.ramPercentageEmployed, ...
                'parforWorkersNumForClassification', params.parforWorkersNumForClassification, ...
                'emPathType', params.emPathType, ...
                'centeredEMPaths', params.centeredEMPaths, ...
                'stimulusDurationInSeconds', params.stimulusDurationInSeconds, ...
                'responseStabilizationMilliseconds', params.responseStabilizationMilliseconds, ...
                'responseExtinctionMilliseconds', params.responseExtinctionMilliseconds, ...
                'secondsToInclude', params.secondsToInclude, ...
                'secondsToIncludeOffset', params.secondsToIncludeOffset, ...
                'singleBinTemporalWindowIfPossible', params.singleBinTemporalWindowIfPossible, ...
                'freezeNoise', params.freezeNoise, ...
                'minimumMosaicFOVdegs', params.minimumMosaicFOVdegs, ...
                'integrationTime', params.integrationTimeMilliseconds/1000, ...
                'coneSpacingMicrons', params.coneSpacingMicrons, ...
                'innerSegmentSizeMicrons', params.innerSegmentDiameter, ...
                'conePacking', params.conePacking, ...
                'eccBasedConeQuantalEfficiency', params.eccBasedConeQuantalEfficiency, ...
                'eccBasedMacularPigment', params.eccBasedMacularPigment, ...
                'LMSRatio', params.LMSRatio, ...
                'mosaicRotationDegs', params.mosaicRotationDegs, ...
                'sConeMinDistanceFactor', params.sConeMinDistanceFactor, ...
                'sConeFreeRadiusMicrons', params.sConeFreeRadiusMicrons, ...
                'latticeAdjustmentPositionalToleranceF', params.latticeAdjustmentPositionalToleranceF, ...
                'latticeAdjustmentDelaunayToleranceF', params.latticeAdjustmentDelaunayToleranceF, ...
                'maxGridAdjustmentIterations', params.maxGridAdjustmentIterations, ...   
                'marginF', params.marginF, ...
                'resamplingFactor', params.resamplingFactor, ...                            
                'computeMosaic', false, ...
                'visualizeMosaic', params.visualizeMosaic, ...
                'computeResponses', false, ...
                'visualizeResponses', false, ...
                'visualizeResponsesWithSpatialPoolingSchemeInVideo', false, ...
                'visualizeOuterSegmentFilters',false, ...
                'visualizeSpatialPoolingScheme', params.visualizeSpatialPoolingScheme, ...
                'findPerformance', params.findPerformance, ...
                'visualizePerformance', params.visualizePerformance, ...
                'visualizeKernelTransformedSignals', params.visualizeKernelTransformedSignals, ...
                'performanceSignal' , params.performanceSignal, ...
                'performanceClassifier', params.performanceClassifier, ...
                'performanceTrialsUsed', params.performanceTrialsUsed, ...
                'spatialPoolingKernelParams', params.spatialPoolingKernelParams ...
                );
           varargout{1} = theMosaic;
           varargout{2} = thePsychometricFunctions;
           varargout{3} = theFigData;
    end
end

