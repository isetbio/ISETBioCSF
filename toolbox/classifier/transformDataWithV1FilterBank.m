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

V1filterBank = thresholdParams.spatialPoolingKernel;
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
    
    [noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
     stimData.responseInstanceArray.theMosaicPhotocurrents, noStimDataPCAapproximatedPhotocurrents, stimDataPCAapproximatedPhotocurrents] = computeV1FilterTransformation(V1filterBank, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE);
else
    [noStimData.responseInstanceArray.theMosaicIsomerizations, ...
     stimData.responseInstanceArray.theMosaicIsomerizations, noStimDataPCAapproximatedIsomerizations, stimDataPCAapproximatedIsomerizations] = computeV1FilterTransformation(V1filterBank, ...
            noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE);
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


function [noStimData, stimData, noStimDataPCAapproximation, stimDataPCAapproximation] = ...
    computeV1FilterTransformation(V1filterBank, noStimData, stimData, repsDimension, spatialDimension, temporalDimension, standardizeData)
    
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(V1filterBank.activationFunction, 'energy')
        noStimData = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
    elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
        noStimData = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
    end
    
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(V1filterBank.activationFunction, 'energy')
        stimData = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
    elseif (strcmp(V1filterBank.activationFunction,'fullWaveRectifier'))
        stimData = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
    end
    
    % Clear some RAM space
    clear 'cosFilterLinearActivation'
    clear 'sinFilterLinearActivation'
    
    noStimDataPCAapproximation = [];
    stimDataPCAapproximation =  [];
        
    figure(1); clf;
    subplot(1,2,1);
    plot(stimData', 'k-')
    if (isfield(V1filterBank, 'temporalPCAcoeffs')) && (V1filterBank.temporalPCAcoeffs > 0)
        nTrials = size(noStimData,repsDimension);
        tBins = size(noStimData,2);
        theData = noStimData;
        theData = cat(1, theData, stimData);
        theData = reshape(theData, [nTrials*2 tBins]);
        [theData, temporalPCAs] = transformDataWithPCA(theData, V1filterBank.temporalPCAcoeffs, standardizeData);
        
        noStimData = theData(1:nTrials,:);
        stimData = theData(nTrials + (1:nTrials),:);
        
        noStimDataPCAapproximation = (temporalPCAs * noStimData')';
        stimDataPCAapproximation = (temporalPCAs * stimData')';

        subplot(1,2,2)
        plot(stimDataPCAapproximation', 'r-');
        drawnow
    end
    
end



