function [noStimData, stimData] = transformDataWithV1FilterEnsemble(noStimData, stimData, thresholdParams, paramsList, visualizeTheTransformedSignals)
% [noStimData, stimData] = transformDataWithV1FilterEnsemble(noStimData, stimData, thresholdParams)
% Compute from the raw signal responses (isomerizations/photocurrents) the
% energy response of a V1 quadrature pair filter bank
%

    V1filterEnsemble = thresholdParams.spatialPoolingKernel;

    fprintf('Transforming data via projection to the spatial components of a population of V1-based filters\n');

    fprintf('Input response size: %d %d %d\n', size(stimData,1), size(stimData,2), size(stimData,3));
    
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
            noStimDataPCAapproximatedPhotocurrents, stimDataPCAapproximatedPhotocurrents] = computeV1FilterEnsembleTransformation(V1filterEnsemble, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE);
    else
        [noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            noStimDataPCAapproximatedIsomerizations, stimDataPCAapproximatedIsomerizations] = computeV1FilterEnsembleTransformation(V1filterEnsemble, ...
            noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE);
    end
    
    % Visualize the transformed signals only if we are working on the full temporal response (i.e. no temporal PCA)
    if (visualizeTheTransformedSignals) && (isinf(V1filterBank.temporalPCAcoeffs))
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            hFig = visualizeTransformedEnsembleSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicPhotocurrents, stimData.responseInstanceArray.theMosaicPhotocurrents, thresholdParams.signalSource, stimData.testContrast*100, sprintf('V1 filter (%s,%s)', V1filterBank.type, V1filterBank.activationFunction));
        else
            hFig = visualizeTransformedEnsembleSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicIsomerizations, stimData.responseInstanceArray.theMosaicIsomerizations, thresholdParams.signalSource, stimData.testContrast*100, sprintf('V1 filter (%s,%s)', V1filterBank.type, V1filterBank.activationFunction));
        end
        
        % Save figure
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('%s-based_%s_%s_%s_%.2fshrinkageFactor_outputs', thresholdParams.signalSource, thresholdParams.method, thresholdParams.spatialPoolingKernelParams.type, thresholdParams.spatialPoolingKernelParams.activationFunction, thresholdParams.spatialPoolingKernelParams.shrinkageFactor);
        rwObject.write(fileName, data, paramsList, theProgram, ...
            'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
    end % visualize transformed signals

    fprintf('Output response size: %d %d %d\n', size(stimData,1), size(stimData,2), size(stimData,3));
    
end


function [noStimData, stimData, noStimDataPCAapproximation, stimDataPCAapproximation] = ...
    computeV1FilterEnsembleTransformation(V1filterEnsemble, noStimData, stimData, repsDimension, spatialDimension, temporalDimension, standardizeData)
    
    unitsNum = numel(V1filterEnsemble);
    
    trialsNum = size(noStimData,repsDimension);
    binsNum = size(noStimData,temporalDimension);
    
    noStimDataUnits = zeros(unitsNum, trialsNum, binsNum);
    stimDataUnits = zeros(unitsNum, trialsNum, binsNum);
    
    for unitIndex = 1:unitsNum
        cosFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterEnsemble{unitIndex}.cosPhasePoolingWeights), spatialDimension));
        sinFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterEnsemble{unitIndex}.sinPhasePoolingWeights), spatialDimension));
        if strcmp(V1filterEnsemble{unitIndex}.activationFunction, 'energy')
            noStimDataUnits(unitIndex,:,:) = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(V1filterEnsemble{unitIndex}.activationFunction,'fullWaveRectifier'))
            noStimDataUnits(unitIndex,:,:)  = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        else
            error('Activation function (''%s''), must be either energy ot fullWaveRectifier\n', V1filterEnsemble{unitIndex}.activationFunction);
        end
        
        cosFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterEnsemble{unitIndex}.cosPhasePoolingWeights), spatialDimension));
        sinFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterEnsemble{unitIndex}.sinPhasePoolingWeights), spatialDimension));
        if strcmp(V1filterEnsemble{unitIndex}.activationFunction, 'energy')
            stimDataUnits(unitIndex,:,:) = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(V1filterEnsemble{unitIndex}.activationFunction,'fullWaveRectifier'))
            stimDataUnits(unitIndex,:,:)  = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        else
            error('Activation function (''%s''), must be either energy ot fullWaveRectifier\n', V1filterEnsemble{unitIndex}.activationFunction);
        end
        
        % Put back in expected shape [trials spatial temporal]
        noStimData = permute(noStimDataUnits, [2 1 3]); 
        stimData = permute(stimDataUnits, [2 1 3]); 
        
        noStimDataPCAapproximation = [];
        stimDataPCAapproximation =  [];
        
        % No PCA for now
    end
    
end
