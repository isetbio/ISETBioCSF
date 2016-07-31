function [usePercentCorrect,useStdErr] = ClassifyForOneDirection(stimData,data,classes,thresholdParams)
% [usePercentCorrect,useStdErr] = ClassifyForOneDirection(stimData,data,classes,thresholdParams)
%
% Called within a parfor loop in t_colorGaborDetectFindPerformance.  We had
% to put this into a function to avoid some parfor errors.
%
% 7/14/16  dhb, xd  Pulled this out.  Hoping it will work.

% For each direction/contrast, insert stimulus data into the array that already has
% the blank data in it.
%
% We can do this to simulate a one interval or a two interval task.  In the
% one interval task, the blanks and modulation instances are labelled as the
% two classes.  In the two inteval task, we concatenate [blank modulation]
% as one class and [modulation blank] as the other.  The same amount of
% data is used in each case, but the number of training instances is half
% for the two interval case, but with effective response vectors that are
% twice as long.
nTrials = numel(stimData.responseInstanceArray);
responseSize = numel(stimData.responseInstanceArray(1).theMosaicPhotoCurrents(:));
if (thresholdParam.nIntervals == 1)
    for iTrial = 1:nTrials
        % Put data into the right form for SVM.
        % This loop overwrites the stimlus data each time through, a
        % little risky, coding wise, but will work unless someone
        % modifies the data generation tutorial to produce a different
        % number of noisy instances for different test directions or
        % contrasts.
        if (strcmp(thresholdParam.signalSource,'photocurrents'))
            data(nTrials+iTrial,:) = stimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
        else
            data(nTrials+iTrial,:) = stimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
        end
    end
elseif (thresholdParam.nIntervals == 2)
    for iTrial = 1:nTrials/2
        if (strcmp(thresholdParam.signalSource,'photocurrents'))
            data(iTrial,responseSize+1:end) = stimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
            data(nTrials/2+iTrial,1:responseSize) = stimData.responseInstanceArray(nTrials/2+iTrial).theMosaicPhotoCurrents(:);
        else
            data(iTrial,responseSize+1:end) = stimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
            data(nTrials/2+iTrial,1:responseSize) = stimData.responseInstanceArray(nTrials/2+iTrial).theMosaicIsomerizations(:);
        end
    end
end

% Do PCA?
if (thresholdParam.PCAComponents > 0)
    fprintf('\tDoing PCA ... ');
    theData = transformDataWithPCA(data,thresholdParam.PCAComponents);
    fprintf('done\n');
else
    theData = data;
end

% Perform SVM classification for this stimulus vs the zero contrast stimulus
fprintf('\tRunning SVM ...');
[usePercentCorrect, useStdErr] = classifyWithSVM(theData,classes,thresholdParam.kFold);
fprintf(' correct: %2.2f%%\n', usePercentCorrect*100);
end