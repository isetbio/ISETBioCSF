function [noStimData, stimData] = transformDataWithSpatialPoolingFilter(noStimData, stimData, thresholdParams, paramsList, visualizeSignals)

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
    if (visualizeSignals) 
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            hFig = visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicPhotocurrents, stimData.responseInstanceArray.theMosaicPhotocurrents, thresholdParams.signalSource, stimData.testContrast*100, 'Spatial pooling filter');
        else
            hFig = visualizeTransformedSignals(noStimData.responseInstanceArray.timeAxis, noStimData.responseInstanceArray.theMosaicIsomerizations, stimData.responseInstanceArray.theMosaicIsomerizations, thresholdParams.signalSource, stimData.testContrast*100, 'Spatial pooling filter');
        end
        
        % Save figure
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('%s-based_%s_%.2fshrinkageFactor_outputs', thresholdParams.signalSource, thresholdParams.method, thresholdParams.spatialPoolingKernelParams.shrinkageFactor);
        rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
    end
end

