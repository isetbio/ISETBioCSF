function blockedTrialIndices = computeBlockedTrialIndices(trialBlockSize, nTrials)

    if (trialBlockSize >= 1) 
        blockedTrialIndices = {};
        trialBlocks = ceil(nTrials/trialBlockSize);
        for iTrialBlock = 1:trialBlocks
            firstTrial = trialBlockSize*(iTrialBlock-1) + 1;
            lastTrial = trialBlockSize*(iTrialBlock-1) + trialBlockSize;
            if (lastTrial > nTrials)
                 lastTrial = nTrials;
            end
            blockedTrialIndices{iTrialBlock} = firstTrial:lastTrial;
        end
    else
        blockedTrialIndices{1} = 1:nTrials;
    end
end