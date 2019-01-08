function run_paper2IsomerizationsVsPhotocurrents

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all but
    % the 2 largest), or some specific spatial frequency, like 16
    computationInstance = 16;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = false;
    
    % Whether to compute responses
    computeResponses = true;

    % Pupil diameter to be used
    pupilDiamMm = 3.0;
    
    % Integration time to use, 2.5 is better for capturing fixationalEM dynamics
    integrationTimeMilliseconds = 5.0;
    
    % Init condition index
    condIndex = 0;
    
    if (1==2)
        condIndex = condIndex+1;
        examinedCond(condIndex).label = 'Ideal observer, isomerizations';
        examinedCond(condIndex).performanceClassifier = 'mlpt';
        examinedCond(condIndex).performanceSignal = 'isomerizations';

        condIndex = condIndex+1;
        examinedCond(condIndex).label = 'SVM-Template-Linear, isomerizations';
        examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
        examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
        examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
        examinedCond(condIndex).performanceSignal = 'isomerizations';
    end
    
    condIndex = condIndex+1;
    examinedCond(condIndex).label = 'SVM-Template-Linear, photocurrents';
    examinedCond(condIndex).performanceClassifier = 'svmV1FilterBank';
    examinedCond(condIndex).spatialPoolingKernelParams.type = 'V1CosUnit';
    examinedCond(condIndex).spatialPoolingKernelParams.activationFunction = 'linear';
    examinedCond(condIndex).performanceSignal = 'photocurrents';
    
    % Go
    examinedLegends = {};
    
    for condIndex = 1:numel(examinedCond)
        params = getCSFPaper2DefaultParams(pupilDiamMm, integrationTimeMilliseconds, computationInstance);
        
        cond = examinedCond(condIndex);
        
        examinedLegends{numel(examinedLegends) + 1} = cond.label;
        params.performanceClassifier = cond.performanceClassifier;
        params.performanceSignal = cond.performanceSignal;
        
        if (strcmp(params.performanceClassifier, 'svmV1FilterBank'))
            params.spatialPoolingKernelParams.type = cond.spatialPoolingKernelParams.type;
            params.spatialPoolingKernelParams.activationFunction = cond.spatialPoolingKernelParams.activationFunction;
        end
        
        % Update params
        if (strcmp(cond.label, 'Ideal observer, isomerizations'))
            computePhotocurrents = false;
        else
            computePhotocurrents = true;
        end
        
        params = getRemainingDefaultParams(params, computeResponses, computePhotocurrents);  
        [~,~, theFigData{condIndex}] = run_BanksPhotocurrentEyeMovementConditions(params);
    end % condIndex
    
    
    if (makeSummaryFigure)
    end
    

    
end

function params = getRemainingDefaultParams(params, computeResponses, computePhotocurrents)
                         
    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;
    
    params.computeResponses = computeResponses;
    params.computePhotocurrentResponseInstances = computePhotocurrents && computeResponses;
    params.visualizeResponses = true;
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
end
