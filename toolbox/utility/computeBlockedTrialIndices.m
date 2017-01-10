function [trialBlockSize, blockedTrialIndices] = computeBlockedTrialIndices(trialBlocks, nTrials)

    trialBlockSize = floor(nTrials/trialBlocks);
    blockedTrialIndices = {};
    
    if (trialBlockSize >= 1) 
        for iTrialBlock = 1:trialBlocks
            firstTrial = trialBlockSize*(iTrialBlock-1) + 1;
            lastTrial = trialBlockSize*(iTrialBlock-1) + trialBlockSize;
            if (iTrialBlock == trialBlocks)
                lastTrial = nTrials;
            end
            blockedTrialIndices{iTrialBlock} = firstTrial:lastTrial;
        end
        % flip order so that the last (possibly larger block) is first
        blockedTrialIndices = fliplr(blockedTrialIndices);
    else
        blockedTrialIndices{1} = 1:nTrials;
    end
end

