function [noStimData, stimData] = transformDataWithV1FilterEnsemble(noStimData, stimData, thresholdParams, paramsList, visualizeTheTransformedSignals, parforWorkersNum)
% [noStimData, stimData] = transformDataWithV1FilterEnsemble(noStimData, stimData, thresholdParams)
% Compute from the raw signal responses (isomerizations/photocurrents) the
% energy response of a V1 quadrature pair filter bank
%

    V1filterEnsemble = thresholdParams.spatialPoolingKernel;

    fprintf('Transforming data via projection to the spatial components of a population of V1-based filters\n');
    
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
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE, parforWorkersNum);
    else
        [noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            noStimDataPCAapproximatedIsomerizations, stimDataPCAapproximatedIsomerizations] = computeV1FilterEnsembleTransformation(V1filterEnsemble, ...
            noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            repsDimension, spatialDimension, temporalDimension, thresholdParams.STANDARDIZE, parforWorkersNum);
    end
    
    % Visualize the transformed signals only if we are working on the full temporal response (i.e. no temporal PCA)
    if (visualizeTheTransformedSignals) && (isinf(V1filterEnsemble{1}.temporalPCAcoeffs))
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            hFigs = visualizeTransformedEnsembleSignals(V1filterEnsemble, ...
                noStimData.responseInstanceArray.timeAxis, ...
                noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
                stimData.responseInstanceArray.theMosaicPhotocurrents, ...
                thresholdParams.signalSource, stimData.testContrast*100, ...
                sprintf('%s(%s)', V1filterEnsemble{1}.type, V1filterEnsemble{1}.activationFunction));
        else
            hFigs = visualizeTransformedEnsembleSignals(V1filterEnsemble, ...
                noStimData.responseInstanceArray.timeAxis, ...
                noStimData.responseInstanceArray.theMosaicIsomerizations, ...
                stimData.responseInstanceArray.theMosaicIsomerizations, ...
                thresholdParams.signalSource, stimData.testContrast*100, ...
                sprintf('%s(%s)', V1filterEnsemble{1}.type, V1filterEnsemble{1}.activationFunction));
        end
        
        for bandwidthIndex = 1:numel(hFigs)
            % Save figure
            theProgram = mfilename;
            rwObject = IBIOColorDetectReadWriteBasic;
            data = 0;
            fileName = sprintf('%s-based_%s_%s_%s_bandwidthIndex_%.0f', thresholdParams.signalSource, thresholdParams.method, thresholdParams.spatialPoolingKernelParams.type, thresholdParams.spatialPoolingKernelParams.activationFunction, bandwidthIndex);
            rwObject.write(fileName, data, paramsList, theProgram, ...
                'type', 'NicePlotExportPDF', 'FigureHandle', hFigs(bandwidthIndex), 'FigureType', 'pdf');
        end
        
    end % visualize transformed signals
end


function [noStimData, stimData, noStimDataPCAapproximation, stimDataPCAapproximation] = ...
    computeV1FilterEnsembleTransformation(V1filterEnsemble, noStimData, stimData, repsDimension, spatialDimension, temporalDimension, standardizeData, parforWorkersNum)
    
    unitsNum = numel(V1filterEnsemble);
    trialsNum = size(noStimData,repsDimension);
    binsNum = size(noStimData,temporalDimension);
    
    % Open a parpool with parforWorkersNum
    poolobj = gcp('nocreate');
    previousPoolSize = 1;
    if (~isempty(poolobj)) 
        % there is a parpool open already.
        % get that pool's # of workers and delete
        previousPoolSize = poolobj.NumWorkers;
        delete(poolobj);
        fprintf('Previous parpool had %d workers. Deleting and staring one with %d workers\n', previousPoolSize, parforWorkersNum);
    end
    % Start a new parpool with desired workers num
    poolobj = parpool(parforWorkersNum);
    
    % Transform the noStimData
    noStimDataUnits = zeros(unitsNum, trialsNum, binsNum);
    parfor unitIndex = 1:unitsNum
        cosFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterEnsemble{unitIndex}.cosPhasePoolingWeights), spatialDimension));
        sinFilterLinearActivation = squeeze(sum(bsxfun(@times, noStimData, V1filterEnsemble{unitIndex}.sinPhasePoolingWeights), spatialDimension));
        if strcmp(V1filterEnsemble{unitIndex}.activationFunction, 'energy')
            noStimDataUnits(unitIndex,:,:) = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(V1filterEnsemble{unitIndex}.activationFunction,'fullWaveRectifier'))
            noStimDataUnits(unitIndex,:,:)  = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        else
            error('Activation function (''%s''), must be either energy ot fullWaveRectifier\n', V1filterEnsemble{unitIndex}.activationFunction);
        end
    end
    
    % Put back in expected shape [trials spatial temporal]
    noStimData = permute(noStimDataUnits, [2 1 3]);  
    
    % Transform the stimData
    stimDataUnits = zeros(unitsNum, trialsNum, binsNum);
    parfor unitIndex = 1:unitsNum
        cosFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterEnsemble{unitIndex}.cosPhasePoolingWeights), spatialDimension));
        sinFilterLinearActivation = squeeze(sum(bsxfun(@times, stimData, V1filterEnsemble{unitIndex}.sinPhasePoolingWeights), spatialDimension));
        if strcmp(V1filterEnsemble{unitIndex}.activationFunction, 'energy')
            stimDataUnits(unitIndex,:,:) = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(V1filterEnsemble{unitIndex}.activationFunction,'fullWaveRectifier'))
            stimDataUnits(unitIndex,:,:)  = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        else
            error('Activation function (''%s''), must be either energy ot fullWaveRectifier\n', V1filterEnsemble{unitIndex}.activationFunction);
        end
    end
    
    % Put back in expected shape [trials spatial temporal]
    stimData = permute(stimDataUnits, [2 1 3]); 
    
    % no PCA
    noStimDataPCAapproximation = [];
    stimDataPCAapproximation =  [];
    
    % Restore previous parpool size
    if (previousPoolSize > 1)
        fprintf('Restoring previous parpool size\n');
        delete(poolobj);
        parpool(previousPoolSize);
    end
    
end
