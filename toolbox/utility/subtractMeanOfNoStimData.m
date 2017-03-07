function [noStimData, stimData] = subtractMeanOfNoStimData(noStimData, stimData, signalSource, repsDimension, temporalDimension)

    % Compute mean (over repetitions and time) response from noStimData
    if (strcmp(signalSource,'photocurrents'))
        meanNoStimPhotocurrentsLevels = mean(mean(noStimData.responseInstanceArray.theMosaicPhotocurrents,temporalDimension),repsDimension);
        noStimData.responseInstanceArray.theMosaicPhotocurrents = ...
            bsxfun(@minus,noStimData.responseInstanceArray.theMosaicPhotocurrents, meanNoStimPhotocurrentsLevels);
        stimData.responseInstanceArray.theMosaicPhotocurrents = ...
            bsxfun(@minus,stimData.responseInstanceArray.theMosaicPhotocurrents, meanNoStimPhotocurrentsLevels);
    else
        meanNoStimIsomerizationsLevels = mean(mean(noStimData.responseInstanceArray.theMosaicIsomerizations,temporalDimension),repsDimension);
        noStimData.responseInstanceArray.theMosaicIsomerizations = ...
            bsxfun(@minus,noStimData.responseInstanceArray.theMosaicIsomerizations, meanNoStimIsomerizationsLevels);
        stimData.responseInstanceArray.theMosaicIsomerizations = ...
            bsxfun(@minus,stimData.responseInstanceArray.theMosaicIsomerizations, meanNoStimIsomerizationsLevels);
    end
end
