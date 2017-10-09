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
        noStimData.responseInstanceArray.theMosaicPhotocurrents = squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicPhotocurrents, spatialPoolingKernel.poolingWeights), spatialDimension));
        stimData.responseInstanceArray.theMosaicPhotocurrents = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicPhotocurrents, spatialPoolingKernel.poolingWeights), spatialDimension));
        noStimData.noiseFreePhotocurrents = sum(bsxfun(@times, noStimData.noiseFreePhotocurrents, spatialPoolingKernel.poolingWeights(:)),1);
        stimData.noiseFreePhotocurrents= sum(bsxfun(@times, stimData.noiseFreePhotocurrents, spatialPoolingKernel.poolingWeights(:)),1);
    else
        noStimData.responseInstanceArray.theMosaicIsomerizations= squeeze(sum(bsxfun(@times, noStimData.responseInstanceArray.theMosaicIsomerizations, spatialPoolingKernel.poolingWeights), spatialDimension));
        stimData.responseInstanceArray.theMosaicIsomerizations = squeeze(sum(bsxfun(@times, stimData.responseInstanceArray.theMosaicIsomerizations, spatialPoolingKernel.poolingWeights), spatialDimension));
        noStimData.noiseFreeIsomerizations = sum(bsxfun(@times, noStimData.noiseFreeIsomerizations, spatialPoolingKernel.poolingWeights(:)),1);
        stimData.noiseFreeIsomerizations = sum(bsxfun(@times, stimData.noiseFreeIsomerizations, spatialPoolingKernel.poolingWeights(:)),1);
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
        fileName = sprintf('%s-based_%s_%s_%s_%.2fshrinkageFactor_outputs', thresholdParams.signalSource, thresholdParams.method, thresholdParams.spatialPoolingKernelParams.type, thresholdParams.spatialPoolingKernelParams.activationFunction, thresholdParams.spatialPoolingKernelParams.shrinkageFactor);
        rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
    end
end

