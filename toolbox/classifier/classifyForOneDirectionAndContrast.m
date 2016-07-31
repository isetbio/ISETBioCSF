function [usePercentCorrect,useStdErr] = classifyForOneDirectionAndContrast(stimData,classificationData,classes,thresholdParams)
% [usePercentCorrect,useStdErr] = classifyForOneDirectionAndContrast(stimData,data,classes,thresholdParams)
%
% Do classification for one stimulus color direction and contrast.
%
% 7/14/16  dhb, xd  Pulled this out.  Hoping it will work.

% Insert stimulus responses into data for classifier
[classificationData,classes] = classificationDataStimDataInsert(classificationData,classes,stimData,thresholdParams);

% Do PCA?
if (thresholdParams.PCAComponents > 0)
    fprintf('\tDoing PCA ... ');
    theData = transformDataWithPCA(classificationData,thresholdParams.PCAComponents);
    fprintf('done\n');
else
    theData = classificationData;
end

% Perform SVM classification for this stimulus vs the zero contrast stimulus
fprintf('\tRunning SVM ...');
[usePercentCorrect, useStdErr] = classifyWithSVM(theData,classes,thresholdParams.kFold);
fprintf(' correct: %2.2f%%\n', usePercentCorrect*100);
end