function [classificationData,classes] = classificationDataNoStimDataInitialize(noStimData,thresholdParams)
% [classificationData,classes] = classificationDataNoStimDataInitialize(noStimData,thresholdParams)
%
% Initialize classification data with zero contrast response instances.
%
% We can do this to simulate a one interval or a two interval task.  In the
% one interval task, the blanks and modulation instances are labelled as the
% two classes.  In the two inteval task, we concatenate [blank modulation]
% as one class and [modulation blank] as the other.  The same amount of
% data is used in each case, but the number of training instances is half
% for the two interval case, but with effective response vectors that are
% twice as long.
%
% See also
%   classificationDataStimDataInsert

nTrials = numel(noStimData.responseInstanceArray);
if (strcmp(thresholdParams.signalSource,'irspikes'))
    responseSize = numel(sum(noStimData.responseInstanceArray(1).theIRSpikes,3));   
else
    responseSize = numel(noStimData.responseInstanceArray(1).theMosaicPhotoCurrents(:));
end
if (thresholdParams.nIntervals == 1)
    classificationData = zeros(2*nTrials, responseSize);
    classes = zeros(2*nTrials, 1);
    for iTrial = 1:nTrials
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(iTrial,:) = noStimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
        elseif (strcmp(thresholdParams.signalSource,'irspikes'))
            classificationData(iTrial,:) = noStimData.responseInstanceArray(iTrial).theIRSpikes(:);
        else
            classificationData(iTrial,:) = noStimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
        end
        
        % Set up classes variable
        classes(iTrial,1) = 0;
        classes(nTrials+iTrial,1) = 1;
    end
elseif (thresholdParams.nIntervals == 2)
    classificationData = zeros(nTrials, 2*responseSize);
    classes = zeros(nTrials, 1);
    for iTrial = 1:nTrials/2
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            classificationData(iTrial,1:responseSize) = noStimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
            classificationData(nTrials/2+iTrial,responseSize+1:end) = noStimData.responseInstanceArray(nTrials/2+iTrial).theMosaicPhotoCurrents(:);
        elseif (strcmp(thresholdParams.signalSource,'irspikes'))
            classificationData(iTrial,1:responseSize) = reshape(sum(noStimData.responseInstanceArray(iTrial).theIRSpikes,3),responseSize,1);
            classificationData(nTrials/2+iTrial,responseSize+1:end) = reshape(sum(noStimData.responseInstanceArray(nTrials/2+iTrial).theIRSpikes,3),responseSize,1);
        else
            classificationData(iTrial,1:responseSize) = noStimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
            classificationData(nTrials/2+iTrial,responseSize+1:end) = noStimData.responseInstanceArray(nTrials/2+iTrial).theMosaicIsomerizations(:);
        end
        
        % Set up classes variable
        classes(iTrial,1) = 0;
        classes(nTrials/2+iTrial,1) = 1;
    end
end
