function run_paper2InferenceEngine
% Compute and contrast performance via different inference engines.
%
% Syntax:
%   run_paper2InferenceEngine
%
% Description:
%    Compute and contrast performance via different inference engines.
%
%    The computation is done via the ecc-based cone efficiency & macular pigment
%    mosaic and the default Thibos subject. We use a 5 msec integration
%    time but we will examine the effect of shorter integration times. 
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History
%    01/31/19  npc  Wrote it.

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    
    % Whether to compute responses
    computeMosaic = ~true;
    computeResponses = ~true;
    visualizeResponses = ~true;
    findPerformance = true;
    visualizePerformance = true;
    
    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use: Here set to 5.0 ms, but 2.5 ms may be better 
    % for capturing the dynamics of fixationalEM
    integrationTimeMilliseconds = 5.0;
    
    % How long the stimulus is. We might be changing the duration. 100 ms is the default
    stimulusDurationInSeconds = 100/1000;
    % Frame rate in Hz. 10 Hz, so each frame is 100 msec long
    % Will need to change this to study shorter stimulus durations.
    frameRate = 10; 
    
    % Compute photocurrent responses
    computePhotocurrents = true;
    
    performanceSignal = 'isomerizations'; %'isomerizations'; % 'photocurrents';
    emPathType = 'randomNoSaccades';
    centeredEMPaths =  true; % 'atStimulusModulationMidPoint';
    nTrainingSamples = 1030;
     
    
    % Assemble conditions list to be examined
    % Init condition index
    condIndex = 0;
    

%     condIndex = condIndex+1;
%     examinedCond(condIndex).label = 'SVM-Template-Linear, noEM';
%     examinedCond(condIndex).minimumMosaicFOVdegs = [];  % no minimum mosaic size, so spatial pooling is matched to stimulus
%     examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
%     examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
%     examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
%     examinedCond(condIndex).performanceSignal = performanceSignal;
%     examinedCond(condIndex).emPathType = 'frozen0';
%     examinedCond(condIndex).centeredEMPaths = true;
% %     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).label = 'SVM-Template-Energy-Matched, drift';
%     examinedCond(condIndex).minimumMosaicFOVdegs = [];  % no minimum mosaic size, so spatial pooling is matched to stimulus
%     examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
%     examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
%     examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
%     examinedCond(condIndex).performanceSignal = performanceSignal;
%     examinedCond(condIndex).emPathType = emPathType;
%     examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
% %     
% %     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).label = 'SVM-Template-Energy-0.16degs, drift';
%     examinedCond(condIndex).minimumMosaicFOVdegs = 0.158;   % stimuli smaller than this, will use spatial pooling based on this mosaic size
%     examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
%     examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
%     examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
%     examinedCond(condIndex).performanceSignal = performanceSignal;
%     examinedCond(condIndex).emPathType = emPathType;
%     examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
%     
%     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).label = 'SVM-Template-Energy-0.33degs, drift';
%     examinedCond(condIndex).minimumMosaicFOVdegs = 0.328;   % stimuli smaller than this, will use spatial pooling based on this mosaic size
%     examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
%     examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
%     examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
%     examinedCond(condIndex).performanceSignal = performanceSignal;
%     examinedCond(condIndex).emPathType = emPathType;
%     examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
%     
%     condIndex = condIndex+1;
%     examinedCond(condIndex).label = 'SVM-Template-Energy-0.50degs, drift';
%     examinedCond(condIndex).minimumMosaicFOVdegs = 0.492;   % stimuli smaller than this, will use spatial pooling based on this mosaic size
%     examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
%     examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
%     examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
%     examinedCond(condIndex).performanceSignal = performanceSignal;
%     examinedCond(condIndex).emPathType = emPathType;
%     examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'SVM-Template-Energy-Ensemble-0.33 degs-7, drift';
    examinedCond(condIndex).performanceClassifier = 'svmV1FilterEnsemble';
    examinedCond(condIndex).minimumMosaicFOVdegs = -0.328;   % nagative sign means that stimuli smaller than this, will use spatial ensemble pooling based on this mosaic size
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedCond(condIndex).performanceSignal = performanceSignal;
    examinedCond(condIndex).emPathType = 'randomNoSaccades';
    examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
    examinedCond(condIndex).ensembleFilterParams = struct(...
                        'spatialPositionsNum',  1, ...   % 1 results in a 3x3 grid of spatial pooling templates
                        'cyclesPerRFs', 7, ...           % each template contains 7 cycles of the stimulus
                        'orientations', 0);
                    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'SVM-Template-Energy-Ensemble-0.33 degs-5, drift';
    examinedCond(condIndex).performanceClassifier = 'svmV1FilterEnsemble';
    examinedCond(condIndex).minimumMosaicFOVdegs = -0.328;   % nagative sign means that stimuli smaller than this, will use spatial ensemble pooling based on this mosaic size
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1QuadraturePair';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'energy';
    examinedCond(condIndex).performanceSignal = performanceSignal;
    examinedCond(condIndex).emPathType = 'randomNoSaccades';
    examinedCond(condIndex).centeredEMPaths = centeredEMPaths;
    examinedCond(condIndex).ensembleFilterParams = struct(...
                        'spatialPositionsNum',  1, ...   % 1 results in a 3x3 grid of spatial pooling templates
                        'cyclesPerRFs', 5, ...           % each template contains 5 cycles of the stimulus
                        'orientations', 0);                

    % Go
    examinedLegends = {};

    for condIndex = 1:numel(examinedCond)
        % Get default params
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, frameRate, stimulusDurationInSeconds, computationInstance);
        params.cyclesPerDegreeExamined = [24];

        
        % Update params
        cond = examinedCond(condIndex);
        params.performanceClassifier = cond.performanceClassifier;
        params.performanceSignal = cond.performanceSignal;
        params.emPathType = cond.emPathType;
        params.centeredEMPaths = cond.centeredEMPaths;
        params.nTrainingSamples = nTrainingSamples;
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
    
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
            params.minimumMosaicFOVdegs = cond.minimumMosaicFOVdegs;
        elseif (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
            params.minimumMosaicFOVdegs = cond.minimumMosaicFOVdegs;
        elseif (strcmp(params.performanceClassifier, 'svmV1FilterEnsemble'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
            params.minimumMosaicFOVdegs = cond.minimumMosaicFOVdegs;
            fNames = fieldnames(cond.ensembleFilterParams);
            for fNameIndex = 1:numel(fNames)
                fName = fNames{fNameIndex};
                params.spatialPoolingKernelParams.(fName) = cond.ensembleFilterParams.(fName);
            end
        end
        
        % Update params
        params = getRemainingDefaultParams(params, ...
            computePhotocurrents, computeResponses, computeMosaic, ...
            visualizeResponses, findPerformance, visualizePerformance); 
        
        % Go !
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
   
    if (makeSummaryFigure)
        variedParamName = sprintf('EyeMovements_%s', performanceSignal);
        theRatioLims = [0.07 1.15];
        theRatioTicks = [0.05 0.1 0.2 0.5 1.0];
        formatLabel = 'ComparedToBanksSubjects';
        generateFigureForPaper(theFigData, examinedLegends, variedParamName, formatLabel, ...
            'figureType', 'CSF', ...
            'showSubjectData', ~true, ...
            'showSubjectMeanData', ~true, ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks, ... 
            'theLegendPosition', [0.23,0.885,0.75,0.08], ...    % custom legend position and size
            'paperDir', 'CSFpaper2', ...                        % sub-directory where figure will be exported
            'figureHasFinalSize', true ...                      % publication-ready size
            );
    end
end

function params = getRemainingDefaultParams(params, computePhotocurrents, computeResponses, computeMosaic, visualizeResponses, findPerformance, visualizePerformance)
                         
    
    % Simulation steps to perform
    params.computeMosaic = computeMosaic; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
    params.computePhotocurrentResponseInstances = computePhotocurrents && computeResponses;
    params.visualizeResponses = visualizeResponses;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;
    
    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = findPerformance;
    params.visualizePerformance = visualizePerformance;
    params.deleteResponseInstances = ~true;
end
