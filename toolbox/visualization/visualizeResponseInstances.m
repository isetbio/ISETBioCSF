function visualizeResponseInstances(theMosaic, stimData, noStimData, responseNormalization, condIndex, condsNum, instancesToVisualize, format)
         
    instancesNum = size(stimData.responseInstanceArray.theMosaicIsomerizations,1);
    if (instancesNum < 1)
        return;
    end
    
    if (isa(theMosaic, 'coneMosaicHex'))
        nonNullCones = theMosaic.pattern(theMosaic.pattern>1);
        for coneIndex = 2:4
            submosaicConeIndices{coneIndex} = find(nonNullCones==coneIndex);
        end
    else
        for coneIndex = 2:4
            submosaicConeIndices{coneIndex} = find(theMosaic.pattern==coneIndex);
        end
    end

    timeBins = numel(stimData.responseInstanceArray.timeAxis);
    if (timeBins == 1)
        if (isa(theMosaic, 'coneMosaicHex'))
            coneDims = 2;
        else
            coneDims = [2 3];
        end
    elseif (timeBins > 1)
        if (isa(theMosaic, 'coneMosaicHex'))
            coneDims = 2;
        else
            coneDims = [2 3];
        end
    else
        error('timeBins = %d', timeBins)
    end
    
    photocurrents = [];
    noiseFreePhotocurrents = [];
    if (strcmp(responseNormalization, 'LMSabsoluteResponseBased')) || (strcmp(responseNormalization, 'LMabsoluteResponseBased')) || (strcmp(responseNormalization, 'MabsoluteResponseBased')) 
        % Max from L- and M-cone mosaics
        if (strcmp(responseNormalization, 'LMSabsoluteResponseBased')) 
            normalizationConeIndices = [submosaicConeIndices{2}; submosaicConeIndices{3}; submosaicConeIndices{4}];
        elseif (strcmp(responseNormalization, 'LMabsoluteResponseBased')) 
            normalizationConeIndices = [submosaicConeIndices{2}; submosaicConeIndices{3};];
        elseif (strcmp(responseNormalization, 'MabsoluteResponseBased')) 
            % Max from M-cone mosaic only
            normalizationConeIndices = [submosaicConeIndices{3}];
        else
            error('unknown normalization: ''%s''.', responseNormalization);
        end
        
        if (timeBins == 1) && isa(theMosaic, 'coneMosaicHex')
            stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations';
            noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations';
        end
        
        [absorptions, minAbsorptions, maxAbsorptions] = coneIndicesBasedScaling(stimData.responseInstanceArray.theMosaicIsomerizations, coneDims, normalizationConeIndices, true);
        [noiseFreeIsomerizations, minNoiseFreeIsomerizations, maxNoiseFreeIsomerizations] = coneIndicesBasedScaling(stimData.noiseFreeIsomerizations, coneDims, normalizationConeIndices, false);
        if (~isempty(stimData.noiseFreePhotocurrents))
            [photocurrents, minPhotocurrents, maxPhotocurrents] = coneIndicesBasedScaling(stimData.responseInstanceArray.theMosaicPhotocurrents, coneDims, normalizationConeIndices, true);
            [noiseFreePhotocurrents, minNoiseFreePhotocurrents, maxNoiseFreePhotocurrents] = coneIndicesBasedScaling(stimData.noiseFreePhotocurrents, coneDims, normalizationConeIndices, false);
        end
    elseif strcmp(responseNormalization, 'submosaicBasedZscore')
        if (timeBins == 1) && isa(theMosaic, 'coneMosaicHex')
            stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations';
            noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations';
        end
        
         [absorptions, minAbsorptions, maxAbsorptions] = submosaicBasedZscore(stimData.responseInstanceArray.theMosaicIsomerizations, noStimData.responseInstanceArray.theMosaicIsomerizations, coneDims, submosaicConeIndices, true);
         [noiseFreeIsomerizations, minNoiseFreeIsomerizations, maxNoiseFreeIsomerizations] = submosaicBasedZscore(stimData.noiseFreeIsomerizations, noStimData.noiseFreeIsomerizations,coneDims, submosaicConeIndices, false);
         if (~isempty(stimData.noiseFreePhotocurrents))
             [photocurrents, minPhotocurrents, maxPhotocurrents] = submosaicBasedZscore(stimData.responseInstanceArray.theMosaicPhotocurrents, noStimData.responseInstanceArray.theMosaicPhotocurrents, coneDims,  submosaicConeIndices, true);
             [noiseFreePhotocurrents, minNoiseFreePhotocurrents, maxNoiseFreePhotocurrents] = submosaicBasedZscore(stimData.noiseFreePhotocurrents, noStimData.noiseFreePhotocurrents, coneDims,  submosaicConeIndices, false);
         end
    else
        error('Unknown responseNormalization method: ''%s''.', responseNormalization);
    end

    if (instancesToVisualize > instancesNum)
        instancesToVisualize = instancesNum;
    end
     
    if (isa(theMosaic, 'coneMosaicHex'))
        for instanceIndex = 1:instancesToVisualize
            if (instanceIndex==1)
                tmp = squeeze(absorptions(instanceIndex,:,:));
                if (timeBins == 1)
                    tmp = tmp';
                end
                tmp = theMosaic.reshapeHex2DmapToHex3Dmap(tmp);
                absorptionsHex = zeros(instancesToVisualize, size(tmp,1), size(tmp,2), size(tmp,3), 'single');
                absorptionsHex(instanceIndex,:,:,:) = tmp;
                
                if (~isempty(photocurrents))
                    tmp = theMosaic.reshapeHex2DmapToHex3Dmap(squeeze(photocurrents(instanceIndex,:,:)));
                    photocurrentsHex = zeros(instancesToVisualize, size(tmp,1), size(tmp,2), size(tmp,3), 'single');
                    photocurrentsHex(instanceIndex,:,:,:) = tmp;
                else
                  photocurrentsHex = [];
                end
            else
                tmp = squeeze(absorptions(instanceIndex,:,:));
                if (timeBins == 1)
                    tmp = tmp';
                end
                absorptionsHex(instanceIndex,:,:,:) = theMosaic.reshapeHex2DmapToHex3Dmap(tmp);
                if (~isempty(photocurrents)) 
                    photocurrentsHex(instanceIndex,:,:,:) = theMosaic.reshapeHex2DmapToHex3Dmap(squeeze(photocurrents(instanceIndex,:,:)));
                end
            end
        end
        
        absorptions = absorptionsHex;
        photocurrents = photocurrentsHex;
        noiseFreeIsomerizations = theMosaic.reshapeHex2DmapToHex3Dmap(noiseFreeIsomerizations);
        if (~isempty(noiseFreePhotocurrents))
            noiseFreePhotocurrents = theMosaic.reshapeHex2DmapToHex3Dmap(noiseFreePhotocurrents);
        end
        activeConesActivations = find(theMosaic.pattern > 1);
        [iRows,iCols] = ind2sub(size(theMosaic.pattern), activeConesActivations);
    end
        
    absorptionsTimeAxis = stimData.responseInstanceArray.timeAxis;
    photocurrentsTimeAxis = stimData.responseInstanceArray.photocurrentTimeAxis;
         
    mosaicXaxis = (squeeze(theMosaic.patternSupport(1,:,1)) + theMosaic.center(1))*1e6;
    mosaicYaxis = (squeeze(theMosaic.patternSupport(:,1,2)) + theMosaic.center(2))*1e6;
    cMapLevels = 1024;
    
    if (isa(theMosaic, 'coneMosaicHex'))
        mosaicXaxis = mosaicXaxis(iCols);
        mosaicYaxis = mosaicYaxis(iRows);
        isHexActivation = true;
        iTheta = (0:10:360)/180*pi;
    else
        iTheta = (0:90:360)/180*pi;
        isHexActivation = false;
    end
    
    % Note that pigment.pdWidth defines the size of a square collective
    % aperture. Here we compute the equivalent circular aperture
    dx = sqrt((theMosaic.pigment.pdWidth^2)/pi)*2;
    apertureOutline.x = dx/2 * cos(iTheta)*1e6;
    apertureOutline.y = dx/2 * sin(iTheta)*1e6;
        
    g = max([1 round(mosaicXaxis/100)]);
    
    xTicks = theMosaic.center(1)*1e6 + g*(-100:50:100);
    yTicks = theMosaic.center(2)*1e6 + g*(-100:50:100);
    xTickLabels = sprintf('%2.0f um\n', xTicks);
    yTickLabels = sprintf('%2.0f um\n', yTicks);
           
    colorbarTicks = 0:0.25:1.0;
    
    if (strcmp(format, 'montage')) && (numel(absorptionsTimeAxis) > 1)
        % Plot of absorptions as a montage for the first instance only
            instanceIndex = 1;
            randomBaseFigNum = round(rand*10000);
            hFig(1) = figure(randomBaseFigNum+condIndex); clf;
            set(hFig(1), 'Position', [10 10 1600 1000], 'Color', [0 0 0], 'Name', sprintf('Absorptions (first instance, integrationTime: %2.2fms); %s', theMosaic.integrationTime*1000, stimData.stimulusLabel));

            if (~isempty(photocurrents))
                randomBaseFigNum = randomBaseFigNum+1000;
                hFig(2) = figure(randomBaseFigNum+condIndex); clf;
                set(hFig(2), 'Position', [10+200 10+20 1600 1000], 'Color', [0 0 0], 'Name', sprintf('Photocurrents (first instance); %s', stimData.stimulusLabel));
            end

            subplotCols = max([1 floor(1.25*sqrt(numel(absorptionsTimeAxis)))]);
            subplotRows = ceil(numel(absorptionsTimeAxis) / subplotCols);
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', subplotRows, ...
                   'colsNum', subplotCols, ...
                   'heightMargin',   0.03, ...
                   'widthMargin',    0.03, ...
                   'leftMargin',     0.03, ...
                   'rightMargin',    0.04, ...
                   'bottomMargin',   0.02, ...
                   'topMargin',      0.02);

    for k = 1:numel(hFig)
        figure(hFig(k));

        for tBin = 1:numel(absorptionsTimeAxis)
            row = floor((tBin-1)/subplotCols) + 1;
            col = mod(tBin-1, subplotCols)+1;
            subplot('Position', subplotPosVectors(row,col).v);
            if (k == 1)
                activation = squeeze(absorptions(instanceIndex, :,:,tBin));
                instanceLabel = sprintf('%2.2fms', absorptionsTimeAxis(tBin)*1000);
                minActivation = minAbsorptions;
                maxActivation = maxAbsorptions;
            else
                activation = squeeze(photocurrents(instanceIndex, :,:,tBin));
                instanceLabel = sprintf('%2.2fms', photocurrentsTimeAxis(tBin)*1000);
                minActivation = minPhotocurrents;
                maxActivation = maxPhotocurrents;
            end
            if (isHexActivation)
                activation = activation(activeConesActivations);
            end
        
            
            signalName = '';
            xTicks = [];
            displayColorBar = false;
            if (tBin == numel(absorptionsTimeAxis))
                displayColorBar = true;
            end
            renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                        responseNormalization, signalName, ...
                        instanceLabel, ...
                        colorbarTicks, minActivation+colorbarTicks*(maxActivation-minActivation), displayColorBar,  ...
                        xTicks, yTicks, xTickLabels, yTickLabels);
        end % tBin
    end % k
    else
        % Plots of absorptions, photocurrents, and their averages
        randomBaseFigNum = round(rand*10000);
        hFig = figure(randomBaseFigNum+condIndex); clf;
        set(hFig, 'Position', [10 10 1100 1050], 'Color', [0 0 0], 'Name', stimData.stimulusLabel);

        if (isempty(photocurrents))
            subplotRows = 1;
        else
            subplotRows = 2;
        end
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', subplotRows, ...
           'colsNum', 2, ...
           'heightMargin',   0.06, ...
           'widthMargin',    0.13, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.10, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
    
        for instanceIndex = 1:instancesToVisualize
         for tBin = 1: numel(absorptionsTimeAxis)  
             
            % Instance absorptions on the left
            subplot('Position', subplotPosVectors(1,1).v);
            activation = squeeze(absorptions(instanceIndex, :,:,tBin));
            if (isHexActivation)
                activation = activation(activeConesActivations);
            end
            renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                responseNormalization, sprintf('absorptions (intTime: %2.2fms)',  theMosaic.integrationTime*1000), ...
                sprintf('instance %d/%d (t: %2.1fms)', instanceIndex, instancesNum, absorptionsTimeAxis(tBin)*1000),  ...
                colorbarTicks, minAbsorptions + colorbarTicks*(maxAbsorptions-minAbsorptions), true, ...
                xTicks, yTicks, xTickLabels, yTickLabels);
           
            % Noise-free isomerizations on the right
            subplot('Position', subplotPosVectors(1,2).v);
            activation = squeeze(noiseFreeIsomerizations(:,:,tBin));
            if (isHexActivation)
                activation = activation(activeConesActivations);
            end
            renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                responseNormalization, sprintf('noise-free absorptions (intTime: %2.2fms)',  theMosaic.integrationTime*1000), ...
                sprintf('cond: %d/%d (t: %2.1fms)', condIndex, condsNum, absorptionsTimeAxis(tBin)*1000), ...
                colorbarTicks, minNoiseFreeIsomerizations + colorbarTicks*(maxNoiseFreeIsomerizations-minNoiseFreeIsomerizations), true, ...
                xTicks, yTicks, xTickLabels, yTickLabels);
 

            if (~isempty(photocurrents))
                % Instance photocurrents on the left
                subplot('Position', subplotPosVectors(2,1).v);
                activation = squeeze(photocurrents(instanceIndex, :,:,tBin));
                if (isHexActivation)
                    activation = activation(activeConesActivations);
                end
                renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                    responseNormalization, 'photocurrents', ...
                    sprintf('instance %d/%d (t: %2.1fms)', instanceIndex, instancesNum, photocurrentsTimeAxis(tBin)*1000), ...
                    colorbarTicks, minPhotocurrents + colorbarTicks*(maxPhotocurrents-minPhotocurrents), true, ...
                    xTicks, yTicks, xTickLabels, yTickLabels);

                % Noise-free photocurrents on the right
                subplot('Position', subplotPosVectors(2,2).v);
                activation = squeeze(noiseFreePhotocurrents(:,:,tBin));
                if (isHexActivation)
                    activation = activation(activeConesActivations);
                end
                renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                    responseNormalization, 'noise-free photocurrents', ...
                    sprintf('cond: %d/%d (t: %2.1fms)', condIndex, condsNum, photocurrentsTimeAxis(tBin)*1000), ...
                    colorbarTicks, minNoiseFreePhotocurrents + colorbarTicks*(maxNoiseFreePhotocurrents-minNoiseFreePhotocurrents), true, ...
                    xTicks, yTicks, xTickLabels, yTickLabels);  
            end
            
            colormap(gray(cMapLevels));
            drawnow;
         end % tBin
        end % instanceIndex 
    end

end 
    
function renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                responseNormalization, signalName, instanceLabel, ...
                colorbarTicks, colorbarTickLabels, displayColorbar, xTicks, yTicks, xTickLabels, yTickLabels)      
    if (isHexActivation)
        edgeColor = 'none';
        lineWidth = 1.0;
        renderPatchArray(apertureOutline, mosaicXaxis, mosaicYaxis, activation, edgeColor, lineWidth);
    else
        imagesc(mosaicXaxis, mosaicYaxis, activation);
    end
    axis 'image'; axis 'xy'; box 'on';
    set(gca, 'XTick', xTicks, 'YTick', [], 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5]);
    set(gca, 'CLim', [0 1], 'FontSize', 14, 'Color', [0 0 0]);
    title(sprintf('%s\n%s', signalName, instanceLabel), 'Color', [0.8 0.8 0.5]);
    colormap(gray(1024));

    if (displayColorbar)
        % Add colorbar
        originalPosition = get(gca, 'position');

        if strcmp(responseNormalization, 'submosaicBasedZscore')
            colorbarLabel = sprintf('z-score');
        else
            colorbarLabel = sprintf('(R*/cone/integrationTIme)');
        end
        hCbar = colorbar('Ticks', colorbarTicks, 'TickLabels', sprintf('%2.3f\n',colorbarTickLabels));
        hCbar.Orientation = 'vertical'; 
        hCbar.Label.String = colorbarLabel;
        hCbar.FontSize = 14; 
        hCbar.FontName = 'Menlo'; 
        hCbar.FontWeight = 'Bold'; 
        hCbar.Color = [0.5 0.5 0.5];
        % The addition changes the figure size, so undo this change
        newPosition = get(gca, 'position');
        set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
    end

end
    
function renderPatchArray(pixelOutline, xCoords, yCoords, faceColorsNormalizedValues,  edgeColor, lineWidth)
    verticesPerCone = numel(pixelOutline.x);
    verticesList = zeros(verticesPerCone * numel(xCoords), 2);
    facesList = [];
    colors = [];
    for coneIdx = 1:numel(xCoords)
        idx = (coneIdx-1)*verticesPerCone + (1:verticesPerCone);
        verticesList(idx,1) = pixelOutline.x(:) + xCoords(coneIdx);
        verticesList(idx,2) = pixelOutline.y(:) + yCoords(coneIdx);
        facesList = cat(1, facesList, idx);
        colors = cat(1, colors, repmat(faceColorsNormalizedValues(coneIdx), [verticesPerCone 1]));
    end

    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.LineWidth = lineWidth;
    patch(S);
end
    

function [responseDataZscore, minZscore, maxZscore] = submosaicBasedZscore(responseData, noResponseData, coneDims, submosaicConeIndices, isInstanceData)
    if (numel(coneDims) == 2)
        originalResponseDataDims = size(responseData);
        if (isInstanceData)
            responseData = reshape(responseData, [size(responseData,1) size(responseData,2)*size(responseData,3) size(responseData,4)]);
            noResponseData = reshape(noResponseData, [size(noResponseData,1) size(noResponseData,2)*size(noResponseData,3) size(noResponseData,4)]);
        else
            responseData = reshape(responseData, [size(responseData,1)*size(responseData,2) size(responseData,3)]);
            noResponseData = reshape(noResponseData, [size(noResponseData,1)*size(noResponseData,2) size(noResponseData,3)]);
        end
    end
    
    responseDataZscore = responseData;
    for coneIndex = 2:4
        if (isInstanceData)
            subMosaicResponseData = responseData(:, submosaicConeIndices{coneIndex},:);
            subMosaicNoResponseData = noResponseData(:, submosaicConeIndices{coneIndex},:);
        else
            % noise-free data
            subMosaicResponseData = responseData(submosaicConeIndices{coneIndex},:);
            subMosaicNoResponseData = noResponseData(submosaicConeIndices{coneIndex},:);
        end
        
        % subtract mean over all cones of a particular type (across all instances and time) for the noResponse data
        meanSubMosaicNoResponse = mean(subMosaicNoResponseData(:));
        subMosaicResponseData = subMosaicResponseData - meanSubMosaicNoResponse;
        
        % divide by std over all cones of a particular type (across all instances and time) for the noResponse data
        if (isInstanceData)
            stdSubMosaicNoResponse = std(subMosaicNoResponseData(:)); 
            subMosaicResponseData = subMosaicResponseData/stdSubMosaicNoResponse;
        end
        
        if (isInstanceData)
            responseDataZscore(:, submosaicConeIndices{coneIndex},:) = subMosaicResponseData;
        else
            % noise-free data
            responseDataZscore(submosaicConeIndices{coneIndex},:) = subMosaicResponseData;
        end
    end % coneIndex

    % Normalize to [0 .. 1]
    if (isInstanceData)
        maxZscore = prctile(responseDataZscore(:), 98); 
        minZscore = prctile(responseDataZscore(:), 2);
    else
        maxZscore = max(responseDataZscore(:)); 
        minZscore = min(responseDataZscore(:)); 
    end
    responseDataZscore = (responseDataZscore-minZscore)/(maxZscore-minZscore);
    responseDataZscore(responseDataZscore<0) = 0;
    responseDataZscore(responseDataZscore>1) = 1;
    
    % Back to original shape
    if (numel(coneDims) == 2)
        responseDataZscore = reshape(responseDataZscore, originalResponseDataDims);
    end
end
 
 
 function [scaledResp, minResp, maxResp] = coneIndicesBasedScaling(responseData, coneDims, coneIndices, isInstanceData)
    if (numel(coneDims) == 2)
        originalResponseDataDims = size(responseData);
        if (isInstanceData)
            responseData = reshape(responseData, [size(responseData,1) size(responseData,2)*size(responseData,3) size(responseData,4)]);
        else
            responseData = reshape(responseData, [size(responseData,1)*size(responseData,2) size(responseData,3)]);
        end
    end
    
    if (isInstanceData)
        subMosaicResponseData = responseData(:, coneIndices,:);
    else
        % noise-free data
        subMosaicResponseData = responseData(coneIndices,:);
    end
        
    % Normalize to [0 .. 1]
    maxResp = max(subMosaicResponseData(:));
    minResp = min(subMosaicResponseData(:));
    scaledResp = (responseData-minResp)/(maxResp-minResp);
    
    % Back to original shape
    if (numel(coneDims) == 2)
        scaledResp = reshape(scaledResp, originalResponseDataDims);
    end
 end