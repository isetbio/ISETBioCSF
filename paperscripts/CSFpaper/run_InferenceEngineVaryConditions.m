function run_InferenceEngineVaryConditions
% This is the script used to assess the impact of different inference engines on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Mosaic to use
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'; 
    
    % Optics to use
    opticsName = 'ThibosBestPSFSubject3MMPupil';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    % Adjust any params we want to change from their default values
    params.opticsModel = opticsName;
    
    % Chromatic direction params
    params.coneContrastDirection = 'L+M+S';
    params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];

    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;

    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;
    
    params.emPathType = 'random';
    params.centeredEMpaths = true;
    
    examinedInferenceEngines = {...
        'mlpt' ...                  % ideal discriminator
        'svm' ...                   % SVM - PCA
        'svmV1FilterBank' ...       % stimulus-matched (linear)
        'svmV1FilterBank' ...       % stimulus-matched (non-linear, HW)
        'svmV1FilterBank' ...       % stimulus-matched (non-linear, FW)
        'svmV1FilterBank' ...       % quadrature energy
        'svmV1FilterEnsemble' ...
        'svmV1FilterEnsemble' ...
        'svmV1FilterEnsemble' ...
        'svmV1FilterEnsemble' ...
        'svmV1FilterEnsemble' ...
        'svmV1FilterEnsemble' ...
    };
    examinedInferenceEngineLegends = {...
        'ideal observer' ...
        'SVM-PCA' ...
        'SVM-Template-L' ...
        'SVM-Template-NL (HW)' ...
        'SVM-Template-NL (FW)' ...
        'SVM-Template-Q' ...
        'SVM (QPhE) population (1)' ...
        'SVM (QPhE) population (2)' ...
        'SVM (QPhE) population (3)' ...
        'SVM (QPhE) population (4)' ...
        'SVM (QPhE) population (5)' ...
        'SVM (QPhE) population (6)' ...
    };

    idx = 1:6;
    visualizedConditions = idx;
    examinedInferenceEngines = {examinedInferenceEngines{visualizedConditions}};
    examinedInferenceEngineLegends = {examinedInferenceEngineLegends{visualizedConditions}};

    
    ensembleFilterParamsStructs{1} = struct(...
        'spatialPositionsNum',  9, ...
        'cyclesPerRFs', 1.0, ...
        'orientations', 0);
    
    ensembleFilterParamsStructs{2} = struct(...
        'spatialPositionsNum',  9, ...
        'cyclesPerRFs', 1.5, ...
        'orientations', 0);
    
    ensembleFilterParamsStructs{3} = struct(...
        'spatialPositionsNum',  9, ...
        'cyclesPerRFs', 2.0, ...
        'orientations', 0);
    
    ensembleFilterParamsStructs{4} = struct(...
        'spatialPositionsNum',  9, ...
        'cyclesPerRFs', 2.5, ...
        'orientations', 0);
    
    ensembleFilterParamsStructs{5} = struct(...
        'spatialPositionsNum',  9, ...
        'cyclesPerRFs', [1.0 1.5 2.0 2.5], ...
        'orientations', 0);
    
    ensembleFilterParamsStructs{6} = struct(...
        'spatialPositionsNum',  9, ...
        'cyclesPerRFs', [1.0 1.5 2.0 2.5], ...
        'orientations', 7.5*[-1 0 1]);
    
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = true;
    params.deleteResponseInstances = ~true;
    
    % Go
    for engineIndex = 1:numel(examinedInferenceEngines)
        params.performanceClassifier = examinedInferenceEngines{engineIndex};
        
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            if strcmp(examinedInferenceEngineLegends{engineIndex},  'SVM-Template-NL (FW)')
                params.spatialPoolingKernelParams.type = 'V1CosUnit';
                params.spatialPoolingKernelParams.activationFunction = 'fullWaveRectifier';
                
            elseif strcmp(examinedInferenceEngineLegends{engineIndex},  'SVM-Template-NL (HW)')
                params.spatialPoolingKernelParams.type = 'V1CosUnit';
                params.spatialPoolingKernelParams.activationFunction = 'halfWaveRectifier';
                
            elseif strcmp(examinedInferenceEngineLegends{engineIndex},  'SVM-Template-L')
                params.spatialPoolingKernelParams.type = 'V1CosUnit';
                params.spatialPoolingKernelParams.activationFunction = 'linear';
                
            elseif strcmp(examinedInferenceEngineLegends{engineIndex},  'SVM-Template-Q')
                params.spatialPoolingKernelParams.type = 'V1QuadraturePair';
                params.spatialPoolingKernelParams.activationFunction = 'energy';
            end 
            
        elseif strcmp(params.performanceClassifier, 'svmV1FilterEnsemble')
            if strcmp(examinedInferenceEngineLegends{engineIndex}, 'SVM (QPhE) population (1)')
                ensembleFilterParams = ensembleFilterParamsStructs{1};
            elseif strcmp(examinedInferenceEngineLegends{engineIndex}, 'SVM (QPhE) population (2)')
                ensembleFilterParams = ensembleFilterParamsStructs{2};
            elseif strcmp(examinedInferenceEngineLegends{engineIndex}, 'SVM (QPhE) population (3)')
                ensembleFilterParams = ensembleFilterParamsStructs{3};
            elseif strcmp(examinedInferenceEngineLegends{engineIndex}, 'SVM (QPhE) population (4)')
                ensembleFilterParams = ensembleFilterParamsStructs{4};
            elseif strcmp(examinedInferenceEngineLegends{engineIndex}, 'SVM (QPhE) population (5)')
                ensembleFilterParams = ensembleFilterParamsStructs{5};
            elseif strcmp(examinedInferenceEngineLegends{engineIndex}, 'SVM (QPhE) population (6)')
                ensembleFilterParams = ensembleFilterParamsStructs{6};
            else
                error('Not handling case: ''%s''', examinedInferenceEngineLegends{engineIndex});
            end
            
            fNames = fieldnames(ensembleFilterParams);
            for fNameIndex = 1:numel(fNames)
                fName = fNames{fNameIndex};
                params.spatialPoolingKernelParams.(fName) = ensembleFilterParams.(fName);
            end

        end
        [~,~, theFigData{engineIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'InferenceEngine';
        theRatioLims = [0.05 1.0];
        theRatioTicks = [0.05  0.1 0.2 0.5 1.0];
        generateFigureForPaper(theFigData, examinedInferenceEngineLegends, variedParamName, sprintf('%s_%s',mosaicName, opticsName), ...
            'figureType', 'CSF', ...
            'inGraphText', '', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
end

