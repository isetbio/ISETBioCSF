function [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams)
%
% 1/17/17   NPC  wrote it


V1filterBank = thresholdParams.V1filterBank;
standardize = thresholdParams.STANDARDIZE;

if (ndims(noStimData.responseInstanceArray.theMosaicIsomerizations) ~= 3)
    error('transformDataWithV1FilterBank not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicIsomerizations));
end

visualizeTransform = true;
if (visualizeTransform)    
    repIndex = 1;
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        photocurrentsResponseNoStim = squeeze(noStimData.responseInstanceArray.theMosaicPhotocurrents(repIndex,:,:));
        photocurrentsResponseStim = squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(repIndex,:,:));
        photocurrentsRange = [min([min(photocurrentsResponseNoStim) min(photocurrentsResponseStim)]) max([max(photocurrentsResponseNoStim) max(photocurrentsResponseStim)])];
    else
        isomerizationsResponseNoStim = squeeze(noStimData.responseInstanceArray.theMosaicIsomerizations(repIndex,:,:));
        isomerizationsResponseStim = squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(repIndex,:,:));
        isomerizationsRange = [min([min(isomerizationsResponseNoStim) min(isomerizationsResponseStim)]) max([max(isomerizationsResponseNoStim) max(isomerizationsResponseStim)])];
    end
end

% Compute mean (over repetitions and time) response from noStimData
repsDimension = 1;
spatialDimension = 2;
temporalDimension = 3;
if (strcmp(thresholdParams.signalSource,'photocurrents'))
    meanNoStimPhotocurrentsLevels = mean(mean(noStimData.responseInstanceArray.theMosaicPhotocurrents,temporalDimension),repsDimension);
else
    meanNoStimIsomerizationsLevels = mean(mean(noStimData.responseInstanceArray.theMosaicIsomerizations,temporalDimension),repsDimension);
end

if (strcmp(thresholdParams.signalSource,'photocurrents'))
    noStimData.responseInstanceArray.theMosaicPhotocurrents = ...
        bsxfun(@minus,noStimData.responseInstanceArray.theMosaicPhotocurrents, meanNoStimPhotocurrentsLevels);
    stimData.responseInstanceArray.theMosaicPhotocurrents = ...
        bsxfun(@minus,stimData.responseInstanceArray.theMosaicPhotocurrents, meanNoStimPhotocurrentsLevels);
    nTrials = size(noStimData.responseInstanceArray.theMosaicPhotocurrents,repsDimension);
    nTimeBins = size(noStimData.responseInstanceArray.theMosaicPhotocurrents,temporalDimension);
else
    noStimData.responseInstanceArray.theMosaicIsomerizations = ...
        bsxfun(@minus,noStimData.responseInstanceArray.theMosaicIsomerizations, meanNoStimIsomerizationsLevels);
    stimData.responseInstanceArray.theMosaicIsomerizations = ...
        bsxfun(@minus,stimData.responseInstanceArray.theMosaicIsomerizations, meanNoStimIsomerizationsLevels);
    nTrials = size(noStimData.responseInstanceArray.theMosaicIsomerizations,repsDimension);
    nTimeBins = size(noStimData.responseInstanceArray.theMosaicIsomerizations,temporalDimension);
end

% Compute the energy response of the V1 filter bank
netWeight = sqrt(...
    (sum(V1filterBank.cosPhasePoolingWeights(:) .* V1filterBank.cosPhasePoolingWeights(:)))^2 + ...
    (sum(V1filterBank.sinPhasePoolingWeights(:) .* V1filterBank.sinPhasePoolingWeights(:)))^2 ...
    );

V1filterBank.cosPhasePoolingWeights = repmat(V1filterBank.cosPhasePoolingWeights, [nTrials 1 nTimeBins]);
V1filterBank.sinPhasePoolingWeights = repmat(V1filterBank.sinPhasePoolingWeights, [nTrials 1 nTimeBins]);

if (strcmp(thresholdParams.signalSource,'photocurrents'))
    cosFilterLinearActivation = squeeze(sum(noStimData.responseInstanceArray.theMosaicPhotocurrents .* V1filterBank.cosPhasePoolingWeights, spatialDimension));
    sinFilterLinearActivation = squeeze(sum(noStimData.responseInstanceArray.theMosaicPhotocurrents .* V1filterBank.sinPhasePoolingWeights, spatialDimension));
    noStimData.responseInstanceArray.theMosaicPhotocurrents = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2)/netWeight;

    cosFilterLinearActivation = squeeze(sum(stimData.responseInstanceArray.theMosaicPhotocurrents .* V1filterBank.cosPhasePoolingWeights, spatialDimension));
    sinFilterLinearActivation = squeeze(sum(stimData.responseInstanceArray.theMosaicPhotocurrents .* V1filterBank.sinPhasePoolingWeights, spatialDimension));
    stimData.responseInstanceArray.theMosaicPhotocurrents = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2)/netWeight;

    if (standardize)
        % zero mean, unit std
        noStimData.responseInstanceArray.theMosaicPhotocurrents = standardizeResponses(noStimData.responseInstanceArray.theMosaicPhotocurrents);
        stimData.responseInstanceArray.theMosaicPhotocurrents = standardizeResponses(stimData.responseInstanceArray.theMosaicPhotocurrents);
    end
    
    photocurrentsBasedV1Range = [...
        min([min(stimData.responseInstanceArray.theMosaicPhotocurrents(:)) min(noStimData.responseInstanceArray.theMosaicPhotocurrents(:))]) ...
        max([max(stimData.responseInstanceArray.theMosaicPhotocurrents(:)) max(noStimData.responseInstanceArray.theMosaicPhotocurrents(:))])];
else
    cosFilterLinearActivation = squeeze(sum(noStimData.responseInstanceArray.theMosaicIsomerizations .* V1filterBank.cosPhasePoolingWeights, spatialDimension));
    sinFilterLinearActivation = squeeze(sum(noStimData.responseInstanceArray.theMosaicIsomerizations .* V1filterBank.sinPhasePoolingWeights, spatialDimension));
    noStimData.responseInstanceArray.theMosaicIsomerizations = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2)/netWeight;

    cosFilterLinearActivation = squeeze(sum(stimData.responseInstanceArray.theMosaicIsomerizations .* V1filterBank.cosPhasePoolingWeights, spatialDimension));
    sinFilterLinearActivation = squeeze(sum(stimData.responseInstanceArray.theMosaicIsomerizations .* V1filterBank.sinPhasePoolingWeights, spatialDimension));
    stimData.responseInstanceArray.theMosaicIsomerizations  = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2)/netWeight;
    
    if (standardize)
        % zero mean, unit std
        noStimData.responseInstanceArray.theMosaicIsomerizations = standardizeResponses(noStimData.responseInstanceArray.theMosaicIsomerizations);
        stimData.responseInstanceArray.theMosaicIsomerizations = standardizeResponses(stimData.responseInstanceArray.theMosaicIsomerizations);
    end  
    
    isomerizationsBasedV1Range = [ ...
        min([min(stimData.responseInstanceArray.theMosaicIsomerizations(:)) min(noStimData.responseInstanceArray.theMosaicIsomerizations(:))]) 
        max([max(stimData.responseInstanceArray.theMosaicIsomerizations(:)) max(noStimData.responseInstanceArray.theMosaicIsomerizations(:))])];
end


if (visualizeTransform)
    
    hFig = figure(1234); clf;
    set(hFig, 'Position', [10 10 1600 600]);

    subplot(1,4,1);
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        plot(1:nTimeBins, photocurrentsResponseNoStim, 'k-')
        set(gca, 'YLim', photocurrentsRange, 'XLim', [1 nTimeBins]);
        title('photocurrents (noStim)');
    else
        stairs(1:nTimeBins, isomerizationsResponseNoStim', 'k-')
        set(gca, 'YLim', isomerizationsRange,  'XLim', [1 nTimeBins]);
        title('isomerizations (noStim)');
    end

    subplot(1,4,2);
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        plot(1:nTimeBins, photocurrentsResponseStim, 'k-')
        set(gca, 'YLim', photocurrentsRange,  'XLim', [1 nTimeBins]);
        title('photocurrents (Stim)');
    else
        stairs(1:nTimeBins, isomerizationsResponseStim', 'k-')
        set(gca, 'YLim', isomerizationsRange,  'XLim', [1 nTimeBins]);
        title('isomerizations (Stim)');
    end

    subplot(1,4,3);
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        plot(1:nTimeBins, noStimData.responseInstanceArray.theMosaicPhotocurrents, 'k-')
        set(gca, 'YLim', photocurrentsBasedV1Range,  'XLim', [1 nTimeBins]);
        title('photocurrents-based V1 response (noStim)');
    else
        stairs(1:nTimeBins, noStimData.responseInstanceArray.theMosaicIsomerizations', 'k-')
        set(gca, 'YLim', isomerizationsBasedV1Range,  'XLim', [1 nTimeBins]);
        title('isomerizations-based V1 response (noStim)');
    end

    subplot(1,4,4);
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        plot(1:nTimeBins, stimData.responseInstanceArray.theMosaicPhotocurrents, 'k-')
        set(gca, 'YLim', photocurrentsBasedV1Range,  'XLim', [1 nTimeBins]);
        title(sprintf('iphotocurrents-based V1 response (c = %2.4f)', stimData.testContrast));
    else
        stairs(1:nTimeBins, stimData.responseInstanceArray.theMosaicIsomerizations', 'k-')
        set(gca, 'YLim', isomerizationsBasedV1Range,  'XLim', [1 nTimeBins]);
        title(sprintf('isomerizations-based V1 response (c = %2.4f)', stimData.testContrast));
    end

    drawnow;
end

end

function responses = standardizeResponses(responses)

    stdDimensionIndex = 2;
    s = std(responses, 0, stdDimensionIndex);
    index = find(s ~= 0);
    %responses = responses(index,:);
    s = s(index);
    m = mean(responses,stdDimensionIndex);
    responses = (responses - repmat(m,1,size(responses,stdDimensionIndex))) ./ repmat(s,1, size(responses,stdDimensionIndex));
end


