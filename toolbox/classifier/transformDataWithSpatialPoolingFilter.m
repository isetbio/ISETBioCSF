function [noStimData, stimData] = transformDataWithSpatialPoolingFilter(noStimData, stimData, thresholdParams)

    fprintf('Transforming data via projection to the spatial pooling kernel\n');
    
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        if (ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents) ~= 3)
            error('transformDataWithSpatialPoolingFilternot yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents));
        end
    else
        if (ndims(noStimData.responseInstanceArray.theMosaicIsomerizations) ~= 3)
            error('transformDataWithSpatialPoolingFilter not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicIsomerizations));
        end
    end
    
    % Subtract the noStimData so that zero modulation gives zero response for both isomerizations and photocurrrents
    repsDimension = 1;
    spatialDimension = 2;
    temporalDimension = 3;
    [noStimData, stimData] = subtractMeanOfNoStimData(noStimData, stimData, thresholdParams.signalSource, repsDimension, temporalDimension);

        
    spatialPoolingKernel = thresholdParams.spatialPoolingKernel;
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        spatialKernelLinearActivation = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicPhotocurrents, spatialPoolingKernel.poolingWeights), spatialDimension));
        noStimData.responseInstanceArray.theMosaicPhotocurrents = spatialKernelLinearActivation;
        
        spatialKernelLinearActivation = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicPhotocurrents, spatialPoolingKernel.poolingWeights), spatialDimension));
        stimData.responseInstanceArray.theMosaicPhotocurrents = spatialKernelLinearActivation;
    
        if (thresholdParams.STANDARDIZE)
            % zero mean, unit std
            noStimData.responseInstanceArray.theMosaicPhotocurrents = standardizeResponses(noStimData.responseInstanceArray.theMosaicPhotocurrents);
            stimData.responseInstanceArray.theMosaicPhotocurrents = standardizeResponses(stimData.responseInstanceArray.theMosaicPhotocurrents);
        end 
    else
        spatialKernelLinearActivation = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicIsomerizations, spatialPoolingKernel.poolingWeights), spatialDimension));
        noStimData.responseInstanceArray.theMosaicIsomerizations = spatialKernelLinearActivation;
        
        spatialKernelLinearActivation = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicIsomerizations, spatialPoolingKernel.poolingWeights), spatialDimension));
        stimData.responseInstanceArray.theMosaicIsomerizations = spatialKernelLinearActivation;
    
        if (thresholdParams.STANDARDIZE)
            % zero mean, unit std
            noStimData.responseInstanceArray.theMosaicIsomerizations = standardizeResponses(noStimData.responseInstanceArray.theMosaicIsomerizations);
            stimData.responseInstanceArray.theMosaicIsomerizations = standardizeResponses(stimData.responseInstanceArray.theMosaicIsomerizations);
        end 
    end
    
    % Visualize transformed signals
    visualizeSignals = true;
    if (visualizeSignals) 
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicPhotocurrents, stimData.responseInstanceArray.theMosaicPhotocurrents, thresholdParams.signalSource, stimData.testContrast*100, 'Spatial pooling filter');
        else
            visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicIsomerizations, stimData.responseInstanceArray.theMosaicIsomerizations, thresholdParams.signalSource, stimData.testContrast*100, 'Spatial pooling filter');
        end
    end
end

