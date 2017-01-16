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
nTrials = size(stimData.responseInstanceArray.theMosaicIsomerizations,1);
responseSize = numel(stimData.responseInstanceArray.theMosaicIsomerizations(1,:,:,:));
if (thresholdParams.nIntervals == 1)
    for iTrial = 1:nTrials
        % Put data into the right form for SVM.
        % This loop overwrites the stimlus data each time through, a
        % little risky, coding wise, but will work unless someone
        % modifies the data generation tutorial to produce a different
        % number of noisy instances for different test directions or
        % contrasts.
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(nTrials+iTrial,:) = reshape(squeeze(stimData.responseInstanceArray.theMosaicPhotoCurrents(iTrial,:,:,:)), [1 responseSize]);
        else
            classificationData(nTrials+iTrial,:) = reshape(squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(Trial,:,:,:)), [1 responseSize]);
        end
    end
elseif (thresholdParams.nIntervals == 2)
    halfTrials = floor(nTrials/2);
    for iTrial = 1:halfTrials
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(iTrial,responseSize+1:end) = reshape(squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(iTrial,:,:,:)), [1 responseSize]);
            classificationData(halfTrials+iTrial,1:responseSize) = reshape(squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(halfTrials+iTrial,:,:,:)), [1 responseSize]);
        else
            classificationData(iTrial,responseSize+1:end) = reshape(squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(iTrial,:,:,:)), [1 responseSize]);
            classificationData(halfTrials+iTrial,1:responseSize) = reshape(squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(halfTrials+iTrial,:,:,:)), [1 responseSize]);
        end
    end
end