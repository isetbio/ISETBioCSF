function [noStimData, stimData] = transformDataWithGaussianPooling(noStimData, stimData, thresholdParams, paramsList, visualizeTheTransformedSignals, parforWorkersNum)
% [noStimData, stimData] =  transformDataWithGaussPooling(noStimData, stimData, thresholdParams)
% Replace each cone's output with the spatially-pooled output using a
% Gaussian kernel specified in thresholdParams
%


fprintf('Transforming data via Gaussian pooling\n');

if (strcmp(thresholdParams.signalSource,'photocurrents'))
    if (isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
        error('photocurrents has not been computed');
    end
    if (ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents) ~= 3)
        error('transformDataWithGaussPooling not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents));
    end
else
    if (ndims(noStimData.responseInstanceArray.theMosaicIsomerizations) ~= 3)
        error('transformDataWithGaussPooling not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicIsomerizations));
    end
end

repsDimension = 1;
spatialDimension = 2;
temporalDimension = 3;

if (visualizeTheTransformedSignals)
    
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        hFig = visualizeGaussianPooledSignals(thresholdParams.spatialPoolingKernel, ...
            noStimData.responseInstanceArray.timeAxis, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            thresholdParams.signalSource, stimData.testContrast*100, ...
            'photocurrents');
    else
        hFig = visualizeGaussianPooledSignals(thresholdParams.spatialPoolingKernel, ...
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
    fileName = sprintf('%s-based_%s_inputs', thresholdParams.signalSource, thresholdParams.method);
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end


% Compute the transformed signals
if (strcmp(thresholdParams.signalSource,'photocurrents'))
    [noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
     stimData.responseInstanceArray.theMosaicPhotocurrents] = computeGaussianPooling(thresholdParams.spatialPoolingKernel, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            repsDimension, spatialDimension, temporalDimension);
else
    [noStimData.responseInstanceArray.theMosaicIsomerizations, ...
     stimData.responseInstanceArray.theMosaicIsomerizations] = computeGaussianPooling(thresholdParams.spatialPoolingKernel, ...
            noStimData.responseInstanceArray.theMosaicIsomerizations, ...
            stimData.responseInstanceArray.theMosaicIsomerizations, ...
            repsDimension, spatialDimension, temporalDimension);
end


% Visualize the transformed signals only if we are working on the full temporal response (i.e. no temporal PCA)
if (visualizeTheTransformedSignals)
    if (strcmp(thresholdParams.signalSource,'photocurrents'))
        hFig = visualizeGaussianPooledSignals(thresholdParams.spatialPoolingKernel, ...
            noStimData.responseInstanceArray.timeAxis, ...
            noStimData.responseInstanceArray.theMosaicPhotocurrents, ...
            stimData.responseInstanceArray.theMosaicPhotocurrents, ...
            thresholdParams.signalSource, stimData.testContrast*100, ...
            'photocurrents');
    else
        hFig = visualizeGaussianPooledSignals(thresholdParams.spatialPoolingKernel, ...
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
    fileName = sprintf('%s-based_%s_outputs', thresholdParams.signalSource, thresholdParams.method);
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end % visualize transformed signals

end


function [noStimData, stimData] = ...
    computeGaussianPooling(spatialPoolingKernel, noStimData, stimData, repsDimension, spatialDimension, temporalDimension)

    conesNum = size(noStimData,spatialDimension);
    poolingWeights = spatialPoolingKernel.poolingWeights;
    
    tmp = noStimData;
    parfor coneIndex = 1:conesNum
        noStimData(:,coneIndex,:) = sum(bsxfun(@times, tmp, poolingWeights(coneIndex,:)), spatialDimension);
    end
    
    tmp = stimData;
    parfor coneIndex = 1:conesNum
        stimData(:,coneIndex,:) = sum(bsxfun(@times, tmp, poolingWeights(coneIndex,:)), spatialDimension);
    end
end



