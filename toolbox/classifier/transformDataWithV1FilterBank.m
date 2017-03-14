function [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams)
% [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams)
% Compute from the raw signal responses (isomerizations/photocurrents) the
% energy response of a V1 quadrature pair filter bank
%

fprintf('Transforming data via projection to the spatial components of a V1-based filter bank (energy)\n');

if (strcmp(thresholdParams.signalSource,'photocurrents'))
    if (ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents) ~= 3)
        error('transformDataWithV1FilterBank not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents));
    end
else
    if (ndims(noStimData.responseInstanceArray.theMosaicIsomerizations) ~= 3)
        error('transformDataWithV1FilterBank not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicIsomerizations));
    end
end

V1filterBank = thresholdParams.V1filterBank;
if (~ismember(V1filterBank.activationFunction, {'energy', 'fullWaveRectifier'}))
    error('V1filterBank.activationFunction must be set to either ''energy'' or ''fullWaveRectifier''.\n')
end

% Subtract the noStimData so that zero modulation gives zero response for both isomerizations and photocurrrents
repsDimension = 1;
spatialDimension = 2;
temporalDimension = 3;
[noStimData, stimData] = subtractMeanOfNoStimData(noStimData, stimData, thresholdParams.signalSource, repsDimension, temporalDimension);

% Compute the energy response of the V1 filter bank
if (strcmp(thresholdParams.signalSource,'photocurrents'))
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicPhotocurrents, V1filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicPhotocurrents, V1filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(V1filterBank.activationFunction, 'energy')
        noStimData.responseInstanceArray.theMosaicPhotocurrents = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
    elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
        noStimData.responseInstanceArray.theMosaicPhotocurrents = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
    end
    
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicPhotocurrents, V1filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicPhotocurrents, V1filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(V1filterBank.activationFunction, 'energy')
        stimData.responseInstanceArray.theMosaicPhotocurrents = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
    elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
        stimData.responseInstanceArray.theMosaicPhotocurrents = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
    end
    
    % Clear some RAM space
    clear 'cosFilterLinearActivation'
    clear 'sinFilterLinearActivation'
    
    if (thresholdParams.STANDARDIZE)
        % zero mean, unit std
        noStimData.responseInstanceArray.theMosaicPhotocurrents = standardizeResponses(noStimData.responseInstanceArray.theMosaicPhotocurrents);
        stimData.responseInstanceArray.theMosaicPhotocurrents = standardizeResponses(stimData.responseInstanceArray.theMosaicPhotocurrents);
    end    
else
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicIsomerizations, V1filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicIsomerizations, V1filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(V1filterBank.activationFunction, 'energy')
        noStimData.responseInstanceArray.theMosaicIsomerizations = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
    elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
        noStimData.responseInstanceArray.theMosaicIsomerizations = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
    end
    
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicIsomerizations, V1filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicIsomerizations, V1filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(V1filterBank.activationFunction, 'energy')
        stimData.responseInstanceArray.theMosaicIsomerizations  = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
    elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
        stimData.responseInstanceArray.theMosaicIsomerizations  = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
    end
    
    % Clear some RAM space
    clear 'cosFilterLinearActivation'
    clear 'sinFilterLinearActivation'
    
    if (thresholdParams.STANDARDIZE)
        % zero mean, unit std
        noStimData.responseInstanceArray.theMosaicIsomerizations = standardizeResponses(noStimData.responseInstanceArray.theMosaicIsomerizations);
        stimData.responseInstanceArray.theMosaicIsomerizations = standardizeResponses(stimData.responseInstanceArray.theMosaicIsomerizations);
    end      
end

% Visualize transformed signals
visualizeTransformedSignals = true;
if (visualizeTransformedSignals) 
    hFig = figure(1234); clf;
    set(hFig, 'Position', [10 10 400 800]);

    subplot(2,1,1);
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        photocurrentsBasedV1Range = [...
            min([min(stimData.responseInstanceArray.theMosaicPhotocurrents(:)) min(noStimData.responseInstanceArray.theMosaicPhotocurrents(:))]) ...
            max([max(stimData.responseInstanceArray.theMosaicPhotocurrents(:)) max(noStimData.responseInstanceArray.theMosaicPhotocurrents(:))])];
        plot(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicPhotocurrents, 'k-')
        set(gca, 'YLim', photocurrentsBasedV1Range,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)]);
        title(sprintf('NO-STIM\nphotocurrents-based V1 filter response'));
    else
        isomerizationsBasedV1Range = [ ...
            min([min(stimData.responseInstanceArray.theMosaicIsomerizations(:)) min(noStimData.responseInstanceArray.theMosaicIsomerizations(:))]) 
            max([max(stimData.responseInstanceArray.theMosaicIsomerizations(:)) max(noStimData.responseInstanceArray.theMosaicIsomerizations(:))])];
        plot(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicIsomerizations, 'k-')
        set(gca, 'YLim', isomerizationsBasedV1Range,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)]);
        title(sprintf('NO-STIM\nisomerizations-based V1 filter response'));
    end
    ylabel('V1 filter bank energy');
    set(gca, 'FontSize', 14);
    
    subplot(2,1,2);
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        plot(stimData.responseInstanceArray.timeAxis, stimData.responseInstanceArray.theMosaicPhotocurrents, 'k-')
        set(gca, 'YLim', photocurrentsBasedV1Range,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)]);
        title(sprintf('C = %2.5f%%\nphotocurrents-based V1 filter response', stimData.testContrast*100));
    else
        plot(stimData.responseInstanceArray.timeAxis, stimData.responseInstanceArray.theMosaicIsomerizations, 'k-')
        set(gca, 'YLim', isomerizationsBasedV1Range,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)]);
        title(sprintf('C = %2.5f%%\nisomerizations-based V1 filterresponse', stimData.testContrast*100));
    end
    ylabel(sprintf('V1 filter bank output (%s)', V1filterBank.activationFunction));
    xlabel('time (ms)'); 
    set(gca, 'FontSize', 14);
    
    drawnow;
    NicePlot.exportFigToPDF('test.pdf', hFig, 300);    
end % visualize transformed signals

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


