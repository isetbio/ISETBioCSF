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
% responseSize = numel(stimData.responseInstanceArray(1).theMosaicPhotoCurrents(:));

if (strcmp(thresholdParams.signalSource,'irspikes'))
    responseSize = numel(sum(stimData.responseInstanceArray(1).theIRSpikes,3)); 
else
    responseSize = numel(stimData.responseInstanceArray(1).theMosaicPhotoCurrents(:));
end

if (thresholdParams.nIntervals == 1)
    for iTrial = 1:nTrials
        % Put data into the right form for SVM.
        % This loop overwrites the stimlus data each time through, a
        % little risky, coding wise, but will work unless someone
        % modifies the data generation tutorial to produce a different
        % number of noisy instances for different test directions or
        % contrasts.
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(nTrials+iTrial,:) = stimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
        elseif (strcmp(thresholdParams.signalSource,'irspikes'))
            classificationData(iTrial,:) = stimData.responseInstanceArray(iTrial).theIRSpikes(:);
        else
            classificationData(nTrials+iTrial,:) = stimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
        end
    end
elseif (thresholdParams.nIntervals == 2)
    for iTrial = 1:nTrials/2
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(iTrial,responseSize+1:end) = stimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
            classificationData(nTrials/2+iTrial,1:responseSize) = stimData.responseInstanceArray(nTrials/2+iTrial).theMosaicPhotoCurrents(:);
        elseif (strcmp(thresholdParams.signalSource,'irspikes'))
            classificationData(iTrial,responseSize+1:end) = reshape(sum(stimData.responseInstanceArray(iTrial).theIRSpikes,3),responseSize,1);
            classificationData(nTrials/2+iTrial,1:responseSize) = reshape(sum(stimData.responseInstanceArray(iTrial).theIRSpikes,3),responseSize,1);
        else
            classificationData(iTrial,responseSize+1:end) = stimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
            classificationData(nTrials/2+iTrial,1:responseSize) = stimData.responseInstanceArray(nTrials/2+iTrial).theMosaicIsomerizations(:);
        end
    end
end