% Examine SVM performance as a function  of  trialsNum
function testSVM

    clear; close all;

    % Stimulus spatial  frequency
    testSF = 2.0;
    
    % Stimulus size (inversely  proportional  to SF)
    sizeDegs = 2.0/testSF;
    
    theMosaic = coneMosaicTreeShrewCreate(75, ...% theOI.optics.micronsPerDegree, ...
        'fovDegs', sizeDegs, ...        % match mosaic width to stimulus size
        'integrationTimeSeconds', 10/1000);
    
    trialsNumList = [100 200 400 800 1600 3200 6400 12800];
    testContrasts  = [0.01 0.015 0.02 0.03];
    
    for k = 1:numel(testContrasts)
        [testTrialsNumList,percentCorrectMean,percentCorrectSEM] = ...
            doAnalysisForContrast(theMosaic, testContrasts(k), testSF, sizeDegs, trialsNumList);
        figure(1)
        subplot(2,2,k);
        plotData(testTrialsNumList,percentCorrectMean,percentCorrectSEM, testContrasts(k));
    end

end


function plotData(trialsNumList,percentCorrectMean,percentCorrectSEM, testContrast  )
    errorbar(trialsNumList,percentCorrectMean,percentCorrectSEM, 's-', ...
        'MarkerSize', 12, 'LineWidth',  1.5);
    set(gca, 'XLim', [0.5*trialsNumList(1) 1.1*trialsNumList(end)], 'XScale', 'log', ...
        'XTick',  [10 30 100 300 1000], 'YLim', [40 100], 'YTick', 0:20:100, ...
        'FontSize', 12);
    grid on;
    title(sprintf('contrast = %2.1f%%',  testContrast*100));
    xlabel('testTrialsNum')
    ylabel('% Correct')
    
end

function [testTrialsNumList,percentCorrectMean,percentCorrectSEM] = ...
    doAnalysisForContrast(theMosaic, testContrast, testSF, sizeDegs, trialsNumList)
    % Create presentation display and place it 5 cm in front of the eye
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 20/100);

    % Stimulus  contrast
    contrast = testContrast;

    % parameter struct for a low spatial frequency Gabor stimulus
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', testSF, ... % changing cycles/deg
        'orientationDegs', 0, ...               % 0 degrees
        'phaseDegs', 0, ...                    % spatial phase degrees, 0 = cos, 90 = sin
        'sizeDegs', sizeDegs, ...                     % 14 x 14 size
        'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
        'contrast', contrast,...                   % 0.005 Michelson contrast
        'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
        'pixelsAlongWidthDim', [], ...          % pixels- width dimension
        'pixelsAlongHeightDim', [] ...          % pixel- height dimension
        );

    %
    % Generate a scene representing the Gabor stimulus with the above params as
    % realized on the presentationDisplay
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);

    % zero contrast for the null stimulus
    stimParams.contrast = 0.0;

    nullScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);

    % Generate wavefront-aberration derived ts optics
    theOI = oiTreeShrewCreate();

    % Compute the retinal image of the test stimulus
    theTestOI = oiCompute(theOI, testScene);

    % Compute the retinal image of the null stimulus
    theNullOI = oiCompute(theOI, nullScene);

    % Obtain the indices of the grid nodes that contain cones
    [~,~,~, nonNullConeIndices] = theMosaic.indicesForCones;
        
    % Examine performance as a function of trialsNum
    percentCorrectMean = zeros(1, numel(trialsNumList));
    percentCorrectSEM = zeros(1, numel(trialsNumList));
    
    for iTrials = 1:numel(trialsNumList)
        nTrialsNum = trialsNumList(iTrials);
        fprintf('\n Starting computation of %d trials ...', nTrialsNum);
        
        emPathLength = 1;
        emPath = zeros(nTrialsNum, emPathLength, 2);

        % Compute mosaic excitation responses to the test stimulus
        coneExcitationsTest = theMosaic.compute(theTestOI, 'emPath', emPath);
        
        % Compute mosaic excitation responses to the null stimulus
        coneExcitationsNull = theMosaic.compute(theNullOI, 'emPath', emPath);

        taskIntervals = 2;
        kFold = 10;
        [percentCorrectMean(iTrials), percentCorrectSEM(iTrials)] = ...
            computePerformance(coneExcitationsTest, coneExcitationsNull, nonNullConeIndices, taskIntervals, kFold); 
        
        fprintf('Done!');
    end
    
    % We are doing a k-Fold cross-validation, so SVM is trained on
    % (kFold-1)/kFold of the trials
    % and performance is tested on the remaining 1/kFold of the trials
    testTrialsNumList = trialsNumList/kFold;
    
end

function [percentCorrectMean, percentCorrectSEM] = computePerformance(coneExcitationsTest, coneExcitationsNull, nonNullConeIndices, ...
    taskIntervals, kFold)

    % Extract the response vectors for nodes containing cones
    [nTrials, nRows, mCols, nTimeBins] = size(coneExcitationsTest);
    coneExcitationsTestReshaped = reshape(coneExcitationsTest, [nTrials nRows*mCols nTimeBins]);
    coneExcitationsNullReshaped = reshape(coneExcitationsNull, [nTrials nRows*mCols nTimeBins]);
    testResponses = coneExcitationsTestReshaped(:, nonNullConeIndices, :);
    nullResponses = coneExcitationsNullReshaped(:, nonNullConeIndices, :);

    % Collapse response vectors across space and time
    responseSize = numel(nonNullConeIndices)*nTimeBins;
    testResponses = reshape(testResponses, [nTrials responseSize]);
    nullResponses = reshape(nullResponses, [nTrials responseSize]);

    %
    if (taskIntervals == 1)
        % In the one interval task, the null and test response instances are labelled as the 2 classes.
        % Allocate matrices
        classificationMatrix = nan(2*nTrials, responseSize);
        classes = nan(2*nTrials, 1);
        % Class 1
        classificationMatrix(1:nTrials,:) = nullResponses;
        classes((1:nTrials)) = 0;
        % Class 2
        classificationMatrix(nTrials+(1:nTrials),:) = testResponses;
        classes(nTrials+(1:nTrials)) = 1;
    elseif (taskIntervals == 2)
        % In the two inteval task, we concatenate [null test] as one class and [test null] as the other.
        % Allocate matrices
        classificationMatrix = nan(nTrials, 2*responseSize);
        classes = nan(nTrials, 1);
        halfTrials = floor(nTrials/2);
        % Class 1
        classificationMatrix(1:halfTrials,:) = [...
            nullResponses(1:halfTrials,:) ...
            testResponses(1:halfTrials,:)];
        classes((1:halfTrials)) = 0;
        % Class 2
        idx = halfTrials+(1:halfTrials);
        classificationMatrix(idx,:) = [...
            testResponses(idx,:) ...
            nullResponses(idx,:)];
        classes(idx) = 1;
    else
        error('Task can have 1 or 2 intervals only.')
    end

    % Find principal components of the responses
    [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

    % Project the responses onto the space formed by the first 2 PC vectors
    pcComponentsNum = 2;
    classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNum);

    % Train a binary SVM classifier
    svm = fitcsvm(classificationMatrixProjection,classes);

    % Perform a 10-fold cross-validation on the trained SVM model
    CVSVM = crossval(svm,'KFold',kFold);

    % Compute classification loss for the in-sample responses using a model
    % trained on out-of-sample responses
    %
    fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
    
    %
    % Mean percent correct across all folds
    percentCorrectMean = mean(fractionCorrect)*100;
    
    %  Standard  error  of the mean across all folds
    percentCorrectSEM = std(fractionCorrect)/sqrt(numel(fractionCorrect))*100;
end