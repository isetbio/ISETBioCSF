function [classificationData,classes] = classificationDataStimDataInsert(classificationData,classes,stimData,thresholdParams)
% [classificationData,classes] = classificationDataStimDataInsert(classificationData,classes,stimData,thresholdParams)
%
% Insert the stimulus data into the classificationData and classes arrays
% that were set up by
%   classificationDataNoStimDataInitialize
%
% See also
%   classificationDataNoStimDataInitialize

% For each direction/contrast, insert stimulus data into the array that already has
% the blank data in it.
nTrials = numel(stimData.responseInstanceArray);
responseSize = numel(stimData.responseInstanceArray(1).theMosaicIsomerizations(:));
if (thresholdParams.nIntervals == 1)
    for iTrial = 1:nTrials
        % Put data into the right form for SVM.
        % This loop overwrites the stimlus data each time through, a
        % little risky, coding wise, but will work unless someone
        % modifies the data generation tutorial to produce a different
        % number of noisy instances for different test directions or
        % contrasts.
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(nTrials+iTrial,:) = reshape(stimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:), [1 responseSize]);
        else
            classificationData(nTrials+iTrial,:) = reshape(stimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:), [1 responseSize]);
        end
    end
elseif (thresholdParams.nIntervals == 2)
    for iTrial = 1:nTrials/2
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(iTrial,responseSize+1:end) = reshape(stimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:), [1 responseSize]);
            classificationData(nTrials/2+iTrial,1:responseSize) = reshape(stimData.responseInstanceArray(nTrials/2+iTrial).theMosaicPhotoCurrents(:), [1 responseSize]);
        else
            classificationData(iTrial,responseSize+1:end) = reshape(stimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:), [1 responseSize]);
            classificationData(nTrials/2+iTrial,1:responseSize) = reshape(stimData.responseInstanceArray(nTrials/2+iTrial).theMosaicIsomerizations(:), [1 responseSize]);
        end
    end
end