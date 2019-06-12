function estimatePerformanceForStimulus(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, contrastLevels, resourcesDir, figNo)
        
    % Load energy mechanism responses to the standard orientation stimulus
    fName = energyResponsesDataFileName(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    load(fName, 'energyConeExcitationResponseOutput', 'energyConeExcitationResponseOrthoOutput', ...
                'energyPhotoCurrentResponseOutput', 'energyPhotoCurrentResponseOrthoOutput', 'timeAxis');
    
    coneExcitationResponseStandardOriStimulusOutput = energyConeExcitationResponseOutput;
    coneExcitationResponseStandardOriStimulusOrthoOutput = energyConeExcitationResponseOrthoOutput;
    photoCurrentResponseStandardOriStimulusOutput = energyPhotoCurrentResponseOutput;
    photoCurrentResponseStandardOriStimulusOrthoOutput = energyPhotoCurrentResponseOrthoOutput;
    
    % Load energy mechanism responses to the orthogonal orientation stimulus
    fName = energyResponsesDataFileName(sprintf('%sOrtho',stimDescriptor), analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    load(fName, 'energyConeExcitationResponseOutput', 'energyConeExcitationResponseOrthoOutput', ...
                'energyPhotoCurrentResponseOutput', 'energyPhotoCurrentResponseOrthoOutput', 'timeAxis');
    coneExcitationResponseOrthogonalOriStimulusOutput = energyConeExcitationResponseOutput;
    coneExcitationResponseOrthogonalOriStimulusOrthoOutput = energyConeExcitationResponseOrthoOutput;
    photoCurrentResponseOrthogonalOriStimulusOutput = energyPhotoCurrentResponseOutput;
    photoCurrentResponseOrthogonalOriStimulusOrthoOutput = energyPhotoCurrentResponseOrthoOutput;
    
    % Find discriminability contrast threshold at the level of cone excitations
    performanceAtConeExcitations = computeDiscriminabilityContrastThreshold(stimDescriptor, 'cone excitations', ...
        coneExcitationResponseStandardOriStimulusOutput, coneExcitationResponseStandardOriStimulusOrthoOutput, ...
        coneExcitationResponseOrthogonalOriStimulusOutput, coneExcitationResponseOrthogonalOriStimulusOrthoOutput, ...
        timeAxis, contrastLevels, figNo);
    
    % Find discriminability contrast threshold at the level of photocurrents
    performanceAtPhotocurrents = computeDiscriminabilityContrastThreshold(stimDescriptor, 'photocurrents', ...
        photoCurrentResponseStandardOriStimulusOutput , photoCurrentResponseStandardOriStimulusOrthoOutput, ...
        photoCurrentResponseOrthogonalOriStimulusOutput, photoCurrentResponseOrthogonalOriStimulusOrthoOutput, ...
        timeAxis, contrastLevels, figNo+100);
    
    fName = performanceDataFileName(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    save(fName, 'performanceAtConeExcitations', 'performanceAtPhotocurrents');
    fprintf('Saved performance for %s to %s\n', stimDescriptor, fName);
    
    figNo = 4000;
    plotPerformance(performanceAtConeExcitations, stimDescriptor, eyePosition, 'cone excitations', figNo);
end

function performance = computeDiscriminabilityContrastThreshold(stimDescriptor, signalName, ...
    standardOriStimulusResponse, standardOriStimulusOrthoResponse, ...
    orthogonalOriStimulusResponse, orthogonalOriStimulusOrthoResponse, timeAxis, contrastLevels, figNo)
    
    visualizeEnergyResponses(stimDescriptor, signalName, ...
        standardOriStimulusResponse, standardOriStimulusOrthoResponse, ...
        orthogonalOriStimulusResponse, orthogonalOriStimulusOrthoResponse, ...
        contrastLevels, timeAxis, figNo);
        
    nContrasts = size(standardOriStimulusResponse,1);
    nTrials = size(standardOriStimulusResponse,2);
    nTimeBins = size(standardOriStimulusResponse,3);
    
    % Use all time bins for now
    timeBinsIncludedInClassification = 1:nTimeBins;
    
    % The psychometric function
    rawPsychometricFunction.contrast = contrastLevels;
    rawPsychometricFunction.performance = zeros(1, nContrasts);
    
    taskIntervals = 2;
    kFold = 10;
    for theContrastLevel = 1:nContrasts
        % The standard and orthogonal tuned mechanism responses to the standard stimulus
        r11 = squeeze(standardOriStimulusResponse(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        r12 = squeeze(standardOriStimulusOrthoResponse(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        % Concatenate responses from 2 mechanisms in 1 long response
        r1 = [r11 r12];

        % The standard and orthogonal tuned mechanism responses to the orthogonal stimulus
        r21 = squeeze(orthogonalOriStimulusResponse(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        r22 = squeeze(orthogonalOriStimulusOrthoResponse(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        % Concatenate responses from 2 mechanisms in 1 long response
        r2 = [r21 r22];

        % Make classification matrix
        [classificationMatrix, classLabels] = assembleBinaryClassificationMatrix(taskIntervals, r1, r2);
        
        % Train a binary SVM classifier 
        svm = fitcsvm(classificationMatrix,classLabels);

        % Perform a 10-fold cross-validation on the trained SVM model
        CVSVM = crossval(svm,'KFold',kFold);

        % Compute classification loss for the in-sample responses using a model trained on out-of-sample responses
        fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
        
        % Average percent correct across all folds 
        rawPsychometricFunction.performance(theContrastLevel) = mean(fractionCorrect)*100;
    end % theContrastLevel
    

    nTrialsForPsychometricFunction = nTrials/kFold;
    performanceCriterion = 0.75;
    [contrastThreshold, smoothPsychometricFunction] = fitWeibulToPsychometricFunction(rawPsychometricFunction.contrast, rawPsychometricFunction.performance, performanceCriterion, nTrialsForPsychometricFunction);

%     % Plot the raw and fitted psychometric function
%     hFig = figure(figNo+1); clf;
%     set(hFig, 'Color', [1 1 1]);
%     
%     % The smooth (fitted) psychometric function
%     plot(smoothPsychometricFunction.contrast*100, smoothPsychometricFunction.performance, 'r-', 'LineWidth', 1.5); hold on;
%     % The raw (measured) psychometric function
%     plot(rawPsychometricFunction.contrast*100, rawPsychometricFunction.performance, 'ko', 'MarkerSize', 12, ...
%         'MarkerFaceColor', [0.8 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'LineWidth', 1.0);
%     % The contrast threshold
%     plot(contrastThreshold*[1 1]*100, [0 performanceThreshold]*100, 'b-', 'LineWidth', 1.5);
%     plot([0.0001 contrastThreshold]*100, performanceThreshold*[1 1]*100, 'b-', 'LineWidth', 1.5);
%     contrastLims = [min(rawPsychometricFunction.contrast) max(rawPsychometricFunction.contrast)]*100;
%     set(gca, 'XScale', 'log', 'XLim', contrastLims, 'YLim', [0.4 1.0]*100, 'XScale', 'log', 'FontSize', 14);
%     xlabel('contrast (%)');
%     ylabel('classification accuracy');
%     title(sprintf('%s (%s)', stimDescriptor, signalName));
    
    performance = struct(...
        'performanceCriterion', performanceCriterion, ...
        'contrastThreshold', contrastThreshold, ...
        'smoothPsychometricFunction', smoothPsychometricFunction, ...
        'rawPsychometricFunction', rawPsychometricFunction);
end


function [contrastThreshold, smoothPsychometricFunction] = fitWeibulToPsychometricFunction(contrasts, rawPsychometricFunction, performanceThreshold, nTrials)
    % Set up psychometric function model. Here we use a cumulative Weibull function
    psychometricFunctionModel = @PAL_Weibull;
    rawPsychometricFunction = rawPsychometricFunction/100;
    % Set up search grid
    gridLevels = 100;
    searchGridParams.alpha = logspace(log10(min(contrasts)),log10(max(contrasts)),gridLevels);
    searchGridParams.beta = 10.^linspace(-4,4,gridLevels);
    searchGridParams.gamma = 0.5;
    searchGridParams.lambda = 0.0;

    % Optimization settings for the fit
    optionsParams             = optimset('fminsearch');
    optionsParams.TolFun      = 1e-09;
    optionsParams.MaxFunEvals = 1000;
    optionsParams.MaxIter     = 1000;
    optionsParams.Display     = 'off';

    % Parameters for the curve fitting
    % Parameters that are allowed to vary
    % The parameters are: threshold, slope, guess-rate, lapse-rate
    paramsFree = [1 1 0 0];
    trialsNumCorrectPerContrastLevel = round(nTrials*rawPsychometricFunction);
    trialsNumPerContrastLevel = repmat(nTrials,1,length(rawPsychometricFunction));
    
    % Fit the data and get the best fit params
    paramsValues = PAL_PFML_Fit(contrasts(:), trialsNumCorrectPerContrastLevel(:), trialsNumPerContrastLevel(:), ...
            searchGridParams, paramsFree, psychometricFunctionModel, 'SearchOptions', optionsParams);
        
    % Obtain the threshold at which performance cross a threshold performance
    contrastThreshold = psychometricFunctionModel(paramsValues, performanceThreshold, 'inverse');
    
    % Obtain a high resolution version of the fitted function
    smoothPsychometricFunction.contrast = searchGridParams.alpha;
    smoothPsychometricFunction.performance = PAL_Weibull(paramsValues, searchGridParams.alpha)*100;

end


function [classificationMatrix, classLabels] = assembleBinaryClassificationMatrix(taskIntervals, responseStandardOriStimulus, responseOrthogonalOriStimulus)

    [nTrials, responseSize] = size(responseStandardOriStimulus);
    
    if (taskIntervals == 1)
        % In the one interval task, the standard and orthogonal orientation response instances are labelled as the 2 classes.
        % Allocate matrices
        classificationMatrix = nan(2*nTrials, responseSize);
        classLabels = nan(2*nTrials, 1);
        % Class 1
        classificationMatrix(1:nTrials,:) = responseOrthogonalOriStimulus;
        classLabels((1:nTrials)) = 0;
        % Class 2
        classificationMatrix(nTrials+(1:nTrials),:) = responseStandardOriStimulus;
        classLabels(nTrials+(1:nTrials)) = 1;
    elseif (taskIntervals == 2)
        % In the two inteval task, we concatenate 
        % [responseStandardOriStimulus responseOrthogonalOriStimulus] as one class and 
        % [responseOrthogonalOriStimulus] as the other. 
        % Allocate matrices
        classificationMatrix = nan(nTrials, 2*responseSize);
        classLabels = nan(nTrials, 1);
        halfTrials = floor(nTrials/2);
        % Class 1
        classificationMatrix(1:halfTrials,:) = [...
            responseOrthogonalOriStimulus(1:halfTrials,:) ...
            responseStandardOriStimulus(1:halfTrials,:)];
        classLabels((1:halfTrials)) = 0;
        % Class 2
        idx = halfTrials+(1:halfTrials);
        classificationMatrix(idx,:) = [...
            responseStandardOriStimulus(idx,:) ...
            responseOrthogonalOriStimulus(idx,:)];
        classLabels(idx) = 1;
    else
        error('Task can have 1 or 2 intervals only.')
    end
end
