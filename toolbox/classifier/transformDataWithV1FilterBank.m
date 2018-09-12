function [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams, paramsList, visualizeTheTransformedSignals)
% [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams)
% Compute from the raw signal responses (isomerizations/photocurrents) the
% energy response of a V1 quadrature pair filter bank
%

V1filterBank = thresholdParams.spatialPoolingKernel;
if (~ismember(V1filterBank.activationFunction, {'linear', 'energy', 'fullWaveRectifier'}))
    error('V1filterBank.activationFunction must be set to either ''linear'', ''energy'', ''halfWaveRectifier''or ''fullWaveRectifier''.\n')
end

fprintf('Transforming data via projection to the spatial components of a V1-based filter (type: %s, activationFunction: %s)\n', V1filterBank.type, V1filterBank.activationFunction);

if (strcmp(thresholdParams.signalSource,'photocurrents'))
    if (isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
        error('photocurrents has not been computed');
    end
    if (ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents) ~= 3)
        error('transformDataWithV1FilterBank not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents));
    end
else
    if (ndims(noStimData.responseInstanceArray.theMosaicIsomerizations) ~= 3)
        error('transformDataWithV1FilterBank not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicIsomerizations));
    end
end

% Subtract the noStimData so that zero modulation gives zero response for both isomerizations and photocurrrents
repsDimension = 1;
spatialDimension = 2;
temporalDimension = 3;
[noStimData, stimData] = subtractMeanOfNoStimData(noStimData, stimData, thresholdParams.signalSource, repsDimension, temporalDimension);

% Compute the transformed signals
if (strcmp(thresholdParams.signalSource,'photocurrents'))
    [noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
     stimData.responseInstanceArray.theMosaicPhotocurrents, ...
     noStimDataPCAapproximatedPhotocurrents, stimDataPCAapproximatedPhotocurrents] = computeV1FilterTransformation(V1filterBank, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE);
else
    [noStimData.responseInstanceArray.theMosaicIsomerizations, ...
     stimData.responseInstanceArray.theMosaicIsomerizations, ...
     noStimDataPCAapproximatedIsomerizations, stimDataPCAapproximatedIsomerizations] = computeV1FilterTransformation(V1filterBank, ...
            noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE);
end


% Visualize the transformed signals only if we are working on the full temporal response (i.e. no temporal PCA)
if (visualizeTheTransformedSignals) && (isinf(V1filterBank.temporalPCAcoeffs))
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
%         hFig = visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, ...
%             noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
%             stimData.responseInstanceArray.theMosaicPhotocurrents, ...
%             thresholdParams.signalSource, stimData.testContrast*100, ...
%             sprintf('V1 filter (%s,%s)', V1filterBank.type, V1filterBank.activationFunction));
        hFig = visualizeTransformedEnsembleSignals(V1filterBank, ...
            noStimData.responseInstanceArray.timeAxis, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            thresholdParams.signalSource, stimData.testContrast*100, ...
            'photocurrents');
    else
%         hFig = visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, ...
%             noStimData.responseInstanceArray.theMosaicIsomerizations, ...
%             stimData.responseInstanceArray.theMosaicIsomerizations, ...
%             thresholdParams.signalSource, stimData.testContrast*100, ...
%             sprintf('V1 filter (%s,%s)', V1filterBank.type, V1filterBank.activationFunction));
        hFig = visualizeTransformedEnsembleSignals(V1filterBank, ...
            noStimData.responseInstanceArray.timeAxis, ...
            noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            thresholdParams.signalSource, stimData.testContrast*100, ...
            'isomerizations');
    end    

    % Save figure
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('%s-based_%s_%s_%s_%.2fshrinkageFactor_outputs', thresholdParams.signalSource, thresholdParams.method, thresholdParams.spatialPoolingKernelParams.type, thresholdParams.spatialPoolingKernelParams.activationFunction, thresholdParams.spatialPoolingKernelParams.shrinkageFactor);
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end % visualize transformed signals

end


function [noStimData, stimData, noStimDataPCAapproximation, stimDataPCAapproximation] = ...
    computeV1FilterTransformation(V1filterBank, noStimData, stimData, repsDimension, spatialDimension, temporalDimension, standardizeData)
    
    if (strcmp(V1filterBank.type, 'V1envelope'))
        % Linear activation function
        noStimData = squeeze(sum(bsxfun(@times, noStimData, V1filterBank.envelopePoolingWeights), spatialDimension));
        stimData = squeeze(sum(bsxfun(@times, stimData, V1filterBank.envelopePoolingWeights), spatialDimension));
    else
        cosFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterBank.cosPhasePoolingWeights), spatialDimension));
        sinFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterBank.sinPhasePoolingWeights), spatialDimension));
        if strcmp(V1filterBank.activationFunction, 'energy')
            noStimData = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
            noStimData = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        elseif (strcmp(V1filterBank.activationFunction,'halfWaveRectifier'))
            noStimData = max(0,cosFilterLinearActivation) + max(0,sinFilterLinearActivation);
        elseif (strcmp(V1filterBank.activationFunction,'linear'))
            noStimData = cosFilterLinearActivation + sinFilterLinearActivation;
        else
            error('Activation function (''%s''), must be one of the following:{''energy'', ''linear'', ''halfWaveRectifier'', or ''fullWaveRectifier''}\n', V1filterBank.activationFunction);
        end

        cosFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterBank.cosPhasePoolingWeights), spatialDimension));
        sinFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterBank.sinPhasePoolingWeights), spatialDimension));
        if strcmp(V1filterBank.activationFunction, 'energy')
            stimData = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
            stimData = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        elseif (strcmp(V1filterBank.activationFunction,'halfWaveRectifier'))
            stimData = max(0,cosFilterLinearActivation) + max(0,sinFilterLinearActivation);
        elseif (strcmp(V1filterBank.activationFunction,'linear'))
            stimData = cosFilterLinearActivation + sinFilterLinearActivation;
        else
            error('Activation function (''%s''), must be one of the following:{''energy'', ''linear'', ''halfWaveRectifier'', or ''fullWaveRectifier''}\n', V1filterBank.activationFunction);
        end
    end
    
    noStimDataPCAapproximation = [];
    stimDataPCAapproximation =  [];
        
  %  figure(1); clf;
  %  subplot(1,2,1);
  %  plot(stimData', 'k-');
    
    if (isfield(V1filterBank, 'temporalPCAcoeffs')) && (~isinf(V1filterBank.temporalPCAcoeffs))
        nTrials = size(noStimData,repsDimension);
        tBins = size(noStimData,2);
        theData = noStimData;
        theData = cat(1, theData, stimData);
        theData = reshape(theData, [nTrials*2 tBins]);
        [theData, temporalPCAs, varianceExplained] = transformDataWithPCA(theData, V1filterBank.temporalPCAcoeffs, standardizeData);
        fprintf('Variance explained: %2.4f (sum of %d elements)\n', sum(varianceExplained), numel(varianceExplained));
        for compIndex = 1:V1filterBank.temporalPCAcoeffs
            fprintf('Variance explained by component #%d: %2.4f (accum: %2.4f)\n', compIndex, varianceExplained(compIndex), sum(varianceExplained(1:compIndex)));
        end
        
        noStimData = theData(1:nTrials,:);
        stimData = theData(nTrials + (1:nTrials),:);
        
        noStimDataPCAapproximation = (temporalPCAs * noStimData')';
        stimDataPCAapproximation = (temporalPCAs * stimData')';

   %     subplot(1,2,2)
   %     plot(stimDataPCAapproximation', 'r-');
   %     drawnow
    end
    
end



