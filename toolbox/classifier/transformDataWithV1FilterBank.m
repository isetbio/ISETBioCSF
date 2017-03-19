function [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams, paramsList, visualizeSignals)
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
if (visualizeSignals) 
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        hFig = visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicPhotocurrents, stimData.responseInstanceArray.theMosaicPhotocurrents, thresholdParams.signalSource, stimData.testContrast*100, 'V1 filter bank');
    else
        hFig = visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicIsomerizations, stimData.responseInstanceArray.theMosaicIsomerizations, thresholdParams.signalSource, stimData.testContrast*100, 'V1 filter bank');
    end    

    % Save figure
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('%s-based_%s_outputs', thresholdParams.signalSource, thresholdParams.method);
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
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


