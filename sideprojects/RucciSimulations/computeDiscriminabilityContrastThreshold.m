function computeDiscriminabilityContrastThreshold(stimDescriptor, signalName, responseStandardOriStimulus, responseOrthogonalOriStimulus, timeAxis, contrastLevels, figNo)
    
    %visualizeEnergyResponses(stimDescriptor, signalName, responseStandardOriStimulus, responseOrthogonalOriStimulus, contrastLevels, timeAxis, figNo);
        
    
    nContrasts = size(responseStandardOriStimulus.output,1);
    nTrials = size(responseStandardOriStimulus.output,2);
    nTimeBins = size(responseStandardOriStimulus.output,3);
    
    timeBinsIncludedInClassification = 1:nTimeBins;
    taskIntervals = 2;
    
    for theContrastLevel = 1:nContrasts
        % The standard and orthogonal tuned mechanism responses to the standard stimulus
        r11 = squeeze(responseStandardOriStimulus.output(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        r12 = squeeze(responseStandardOriStimulus.orthoOutput(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        % Concatenate responses from 2 mechanisms in 1 long response
        r1 = [r11 r12];

        % The standard and orthogonal tuned mechanism responses to the orthogonal stimulus
        r21 = squeeze(responseOrthogonalOriStimulus.output(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        r22 = squeeze(responseOrthogonalOriStimulus.orthoOutput(theContrastLevel,1:nTrials,timeBinsIncludedInClassification));
        % Concatenate responses from 2 mechanisms in 1 long response
        r2 = [r21 r22];

        % Make classification matrix
        [classificationMatrix, classLabels] = assembleBinaryClassificationMatrix(taskIntervals, r1, r2);
        
        size(classificationMatrix)
        size(classLabels)
        
        % Train a binary SVM classifier 
        svm = fitcsvm(classificationMatrix,classLabels);

        % Perform a 10-fold cross-validation on the trained SVM model
        kFold = 10;
        CVSVM = crossval(svm,'KFold',kFold);

        % Compute classification loss for the in-sample responses using a model trained on out-of-sample responses
        fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
        
        % Average percent correct across all folds 
        percentCorrect(theContrastLevel) = mean(fractionCorrect)*100;
    end % theContrastLevel
    
    figure(figNo+1); clf;
    plot(contrastLevels, percentCorrect);
    xlabel('contrast');
    ylabel('classification accuracy');
    title(stimDescriptor);
    
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
