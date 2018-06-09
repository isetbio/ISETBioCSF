function hFigsInfo = visualizeResponseInstances(theMosaic, stimData, noStimData, visualizeOuterSegmentFilters, responseNormalization, condIndex, condsNum, format)


    instancesNum = size(stimData.responseInstanceArray.theMosaicIsomerizations,1);
    if (instancesNum < 1)
        return;
    end
    
    hFigsInfo = {};
    
    if (visualizeOuterSegmentFilters)
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'osImpulseResponses',...
            'hFig', plotImpulseResponseFunctions(stimData));
    end
    
    if (~isempty(noStimData.responseInstanceArray.theMosaicIsomerizations))
        % transform isomerization counts to isomerization rate
        stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations / theMosaic.integrationTime;
        noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations / theMosaic.integrationTime;
        stimData.responseInstanceArray.theMosaicIsomerizations = stimData.responseInstanceArray.theMosaicIsomerizations / theMosaic.integrationTime;
        noStimData.responseInstanceArray.theMosaicIsomerizations = noStimData.responseInstanceArray.theMosaicIsomerizations / theMosaic.integrationTime;
    end

    timeAxis = 1000*noStimData.responseInstanceArray.timeAxis;
    
    if (numel(timeAxis) == 1)
        stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations(:);
        noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations(:);
        stimData.noiseFreePhotocurrents = stimData.noiseFreePhotocurrents(:);
        noStimData.noiseFreePhotocurrents = noStimData.noiseFreePhotocurrents(:);
    end
   
    [hFigs, peakConeIndex] = visualizeNoiseFreeResponses(theMosaic, timeAxis, stimData, noStimData);
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXTResponses',...
            'hFig', hFigs{1});
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXYResponses',...
            'hFig', hFigs{2});
        
    
    for submosaicIndex = 1:3
        if (~isempty(peakConeIndex{submosaicIndex}))
            theSelectedConeIndex = peakConeIndex{submosaicIndex};
            if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
                noStimPhotocurrentsResponseInstances(submosaicIndex,:,:)  = squeeze(noStimData.responseInstanceArray.theMosaicPhotocurrents(:,theSelectedConeIndex,:));
                noStimNoiseFreePhotocurrentsResponse(submosaicIndex,:)    = noStimData.noiseFreePhotocurrents(theSelectedConeIndex,:);
                stimPhotocurrentsResponseInstances(submosaicIndex,:,:)    = squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(:,theSelectedConeIndex,:));
                stimNoiseFreePhotocurrentsResponse(submosaicIndex,:)      = stimData.noiseFreePhotocurrents(theSelectedConeIndex,:);     
            end
            
            if (~isempty(noStimData.responseInstanceArray.theMosaicIsomerizations))
                noStimIsomerizationsResponseInstances(submosaicIndex,:,:) = squeeze(noStimData.responseInstanceArray.theMosaicIsomerizations(:,theSelectedConeIndex,:));
                noStimNoiseFreeIsomerizationsResponse(submosaicIndex,:)   = noStimData.noiseFreeIsomerizations(theSelectedConeIndex,:);
                stimIsomerizationsResponseInstances(submosaicIndex,:,:)   = squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(:,theSelectedConeIndex,:));
                stimNoiseFreeIsomerizationsResponse(submosaicIndex,:)     = stimData.noiseFreeIsomerizations(theSelectedConeIndex,:);
            end
        end
    end
    
    isomerizationQuantizationLevelsNum = 30;
    photocurrentsQuantizationLevelsNum = 100;
    
    if (~isempty(noStimData.responseInstanceArray.theMosaicIsomerizations))
        isomerizationLevels = [0 30]; 
        hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
            timeAxis, ...
            noStimIsomerizationsResponseInstances*theMosaic.integrationTime, ...
            stimIsomerizationsResponseInstances*theMosaic.integrationTime, ...
            noStimNoiseFreeIsomerizationsResponse*theMosaic.integrationTime, ...
            stimNoiseFreeIsomerizationsResponse*theMosaic.integrationTime, ...
            isomerizationLevels, ...
            isomerizationQuantizationLevelsNum, ...
            'density', 'isomerizations', 'R*/cone', 5001);
    else
        hFig = [];
    end 
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'peakConeIsomerizationResponseInstances',...
            'hFig', hFig);
    
    if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
        
        photocurrentsLevels = [-90 -50]; 
        hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
            timeAxis, ...
            noStimPhotocurrentsResponseInstances, ...
            stimPhotocurrentsResponseInstances, ...
            noStimNoiseFreePhotocurrentsResponse, ...
            stimNoiseFreePhotocurrentsResponse, ...
            photocurrentsLevels, ...
            photocurrentsQuantizationLevelsNum, ...
            'line', 'photocurrents', 'pA', 5002);
    else
        hFig = [];
    end
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'peakConePhotocurrentResponseInstances',...
            'hFig', hFig);
        
   end

function hFig = plotImpulseResponseFunctions(stimData)

    if (numel(stimData.osImpulseResponseTimeAxis) > 1)
        % Visualize the os impulse response functions
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 700 500], 'Color', [1 1 1]);
        hold on

        timeAxis = 1000*stimData.osImpulseResponseTimeAxis;
        plot(timeAxis, stimData.osImpulseResponses(:,1), 'r-', 'LineWidth', 1.5);
        plot(timeAxis, stimData.osImpulseResponses(:,2), 'g-', 'LineWidth', 1.5);
        plot(timeAxis, stimData.osImpulseResponses(:,3), 'b-', 'LineWidth', 1.5);
        xlabel('time (msec)');
        set(gca, 'FontSize', 14);
        grid on;
        title(sprintf('os impulse response functions'));
        drawnow;
    else
        hFig = [];
    end
end
    

function [hFig, peakConeIndex] = visualizeNoiseFreeResponses(theMosaic, timeAxis, stimData, noStimData)

    hFig{1} = [];
    hFig{2} = [];
    
    isomerizationsRange = [...
        min([min(noStimData.noiseFreeIsomerizations(:)) min(stimData.noiseFreeIsomerizations(:))]) ...
        max([max(noStimData.noiseFreeIsomerizations(:)) max(stimData.noiseFreeIsomerizations(:))]) ];
    
    if (isomerizationsRange(1) / isomerizationsRange(2) < 0.9)
        isomerizationsRange(2) = isomerizationsRange(2) * 1.1;
        isomerizationsRange(1) = isomerizationsRange(1) / 1.1;
    end
    
    if (isomerizationsRange(1) > isomerizationsRange(2)*0.50)
        isomerizationsRange(1) = isomerizationsRange(2)*0.50;
    end

    if (~isempty(noStimData.noiseFreePhotocurrents))
        photocurrentsRange = [...
            min([min(noStimData.noiseFreePhotocurrents(:)) min(stimData.noiseFreePhotocurrents(:))]) ...
            max([max(noStimData.noiseFreePhotocurrents(:)) max(stimData.noiseFreePhotocurrents(:))]) ];

        deltaRange = photocurrentsRange(2)-photocurrentsRange(1);
        if (deltaRange < 2)
            midRange = 0.5*(photocurrentsRange(2)+photocurrentsRange(1));
            photocurrentsRange(1) = midRange - 1;
            photocurrentsRange(2) = midRange + 1;
        end
    end
    
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',   0.08, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.08, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.03);
      
    
    nonNullCones = theMosaic.pattern(theMosaic.pattern > 1);
    
    LMSindices = [];
    noStimLMSnoiseFreeIsomerizations = [];
    noStimLMSnoiseFreePhotocurrents = [];
    stimLMSnoiseFreeIsomerizations = [];
    stimLMSnoiseFreePhotocurrents = [];
    
    peakConeIndex = cell(1,3);
    for submosaicIndex = 2:4
        submosaicConeIndices = find(nonNullCones==submosaicIndex);
        if (~isempty(submosaicConeIndices))
            if (isempty(LMSindices))
                LMSindices = submosaicConeIndices(:);
                noStimLMSnoiseFreeIsomerizations = noStimData.noiseFreeIsomerizations(submosaicConeIndices,:);
                stimLMSnoiseFreeIsomerizations = stimData.noiseFreeIsomerizations(submosaicConeIndices,:);
                if (~isempty(noStimData.noiseFreePhotocurrents))
                    noStimLMSnoiseFreePhotocurrents  = noStimData.noiseFreePhotocurrents(submosaicConeIndices,:);
                    stimLMSnoiseFreePhotocurrents  = stimData.noiseFreePhotocurrents(submosaicConeIndices,:);
                end
            else
                LMSindices = cat(1, LMSindices, submosaicConeIndices(:));
                noStimLMSnoiseFreeIsomerizations = cat(1, noStimLMSnoiseFreeIsomerizations, noStimData.noiseFreeIsomerizations(submosaicConeIndices,:));
                stimLMSnoiseFreeIsomerizations = cat(1, stimLMSnoiseFreeIsomerizations, stimData.noiseFreeIsomerizations(submosaicConeIndices,:));
                if (~isempty(noStimData.noiseFreePhotocurrents))
                    noStimLMSnoiseFreePhotocurrents  = cat(1, noStimLMSnoiseFreePhotocurrents,  noStimData.noiseFreePhotocurrents(submosaicConeIndices,:));
                    stimLMSnoiseFreePhotocurrents  = cat(1, stimLMSnoiseFreePhotocurrents,  stimData.noiseFreePhotocurrents(submosaicConeIndices,:));
                end
            end
            
            x = stimData.noiseFreeIsomerizations(submosaicConeIndices,:);
            [~, idx] = max(x(:));
            [idx, ~] = ind2sub(size(x), idx);
            peakConeIndex{submosaicIndex-1} = submosaicConeIndices(idx);
        end
    end


    isomerizationModulation = 100*(stimLMSnoiseFreeIsomerizations - noStimLMSnoiseFreeIsomerizations)/noStimLMSnoiseFreeIsomerizations(1);
    isomerizationsRangeModulation = max(abs(isomerizationModulation(:)))*[-1 1];
    
   
    hFig{1} = figure(6001); clf;
    set(hFig{1}, 'Position', [10 10 770 900], 'Color', [1 1 1]);
    subplot('Position', subplotPosVectors(1,1).v);
    
    imagesc(timeAxis, 1:size(noStimData.noiseFreeIsomerizations,1), noStimLMSnoiseFreeIsomerizations);
    set(gca, 'CLim', isomerizationsRange, 'FontSize', 14);
    ylabel('cone #');
    hcb = colorbar('northoutside');
    colorTitleHandle = get(hcb,'Title');
    set(colorTitleHandle ,'String', sprintf('noise-free isomerization rates (R*/cone/sec)'));
    set(hcb, 'FontSize', 12);
    
    subplot('Position', subplotPosVectors(1,2).v);
    imagesc(timeAxis, 1:size(stimData.noiseFreeIsomerizations,1), isomerizationModulation);
    set(gca, 'CLim', isomerizationsRangeModulation, 'YTickLabel', {}, 'FontSize', 14);
    
    hcb = colorbar('northoutside');
    colorTitleHandle = get(hcb,'Title');
    set(colorTitleHandle ,'String', 'noise-free isomerization rate modulation (%%)');
    set(hcb, 'FontSize', 12);
    
    if (~isempty(noStimData.noiseFreePhotocurrents))
        deltaPhotocurrents = stimLMSnoiseFreePhotocurrents-noStimLMSnoiseFreePhotocurrents(1);
        deltaPhotocurrentsRange = [min(deltaPhotocurrents(:)) max(deltaPhotocurrents(:))];
        subplot('Position', subplotPosVectors(2,1).v);
        imagesc(timeAxis, 1:size(noStimData.noiseFreePhotocurrents,1), noStimLMSnoiseFreePhotocurrents);
        set(gca, 'CLim', photocurrentsRange, 'FontSize', 14);
        ylabel('cone #');
        xlabel('time (msec)');
        hcb = colorbar('northoutside');
        colorTitleHandle = get(hcb,'Title');
        set(colorTitleHandle ,'String', sprintf('noise-free photocurrents (pAmps)'));
        set(hcb, 'FontSize', 12);

        subplot('Position', subplotPosVectors(2,2).v);
        imagesc(timeAxis, 1:size(stimData.noiseFreePhotocurrents,1), deltaPhotocurrents);
        set(gca, 'CLim', deltaPhotocurrentsRange, 'YTickLabel', {}, 'FontSize', 14);
        xlabel('time (msec)');
        hcb = colorbar('northoutside');
        colorTitleHandle = get(hcb,'Title');
        set(colorTitleHandle ,'String', sprintf('noise-free delta photocurrents (pAmps)'));
        set(hcb, 'FontSize', 12);
    end
    
    colormap(gray(1024));
    drawnow;
    
    
    if (isa(theMosaic, 'coneMosaicHex'))
       
        stimData.noiseFreeIsomerizationsModulations = 100*(stimData.noiseFreeIsomerizations-noStimData.noiseFreeIsomerizations)./noStimData.noiseFreeIsomerizations;
        isomerizationsRangeModulation = max(abs(stimData.noiseFreeIsomerizationsModulations(:)))*[-1 1];
    
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',   0.07, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.06, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.05);
       
        hFig{2} = figure(6002); clf;
        set(hFig{2}, 'Position', [10 10 900 1000], 'Color', [1 1 1]);
    
        subplot('Position', subplotPosVectors(1,2).v);
        timeBinOfPeakIsomerizationResponse = renderHexActivationMap(theMosaic, stimData.noiseFreeIsomerizationsModulations, isomerizationsRangeModulation, []);
        set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 18);
        hcb = colorbar('northoutside');
        colorTitleHandle = get(hcb,'Title');
        set(colorTitleHandle ,'String', sprintf('test stimulus\nisomerization modulation (%%)'));
        set(hcb, 'FontSize', 18);
        colormap(gca, brewermap(512, '*RdBu'));
        
        subplot('Position', subplotPosVectors(1,1).v);
        renderHexActivationMap(theMosaic, noStimData.noiseFreeIsomerizations, isomerizationsRange, timeBinOfPeakIsomerizationResponse);
        set(gca, 'XTickLabel', {}, 'FontSize', 18);
        hcb = colorbar('northoutside');
        colorTitleHandle = get(hcb,'Title');
        set(colorTitleHandle ,'String', sprintf('null stimulus\nisomerization response (R*/cone/sec)'));
        set(hcb, 'FontSize', 18);
        ylabel('microns', 'FontWeight', 'bold');
        colormap(gca, gray(1024));
        
        if (~isempty(noStimData.noiseFreePhotocurrents))
            deltaPhotocurrents = stimData.noiseFreePhotocurrents-noStimData.noiseFreePhotocurrents;
            subplot('Position', subplotPosVectors(2,2).v);
            timeBinOfPeakPhotocurrentResponse = renderHexActivationMap(theMosaic, deltaPhotocurrents, deltaPhotocurrentsRange, []);
            set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 18);
            hcb = colorbar('northoutside');
            colorTitleHandle = get(hcb,'Title');
            set(colorTitleHandle ,'String', sprintf('noise-free, test stimulus delta photocurrent (pAmps)'));
            set(hcb, 'FontSize', 18);
            colormap(gca, gray(1024));
        
            subplot('Position', subplotPosVectors(2,1).v);
            renderHexActivationMap(theMosaic, noStimData.noiseFreePhotocurrents, photocurrentsRange, timeBinOfPeakPhotocurrentResponse);
            set(gca, 'XTickLabel', {},  'FontSize', 18);
            hcb = colorbar('northoutside');
            colorTitleHandle = get(hcb,'Title');
            set(colorTitleHandle ,'String', sprintf('noise-free, null stimulus photocurrent response (pAmps)'));
            set(hcb, 'FontSize', 18);
            ylabel('microns', 'FontWeight', 'bold');
            colormap(gca, gray(1024));
        end
        
        
        drawnow;
    end
end

function timeBinOfPeakResponse = renderHexActivationMap(theMosaic, signal, signalRange, timeBinOfPeakResponse)

    [~, idx] = max(signal(:));
    if (isempty(timeBinOfPeakResponse))
        [coneIndexOfPeakResponse, timeBinOfPeakResponse] = ind2sub(size(signal), idx);
    end
    
    activeConeIndices = find(theMosaic.pattern > 1);
    hexMap = theMosaic.reshapeHex2DmapToHex3Dmap(signal);
    hexMap = squeeze(hexMap(:, :,timeBinOfPeakResponse));
    hexMap = hexMap(activeConeIndices);
    [iRows,iCols] = ind2sub(size(theMosaic.pattern), activeConeIndices);
    
    mosaicXaxis = (squeeze(theMosaic.patternSupport(1,:,1)) + theMosaic.center(1))*1e6;
    mosaicYaxis = (squeeze(theMosaic.patternSupport(:,1,2)) + theMosaic.center(2))*1e6;
    mosaicXaxis = mosaicXaxis(iCols);
    mosaicYaxis = mosaicYaxis(iRows);
    iTheta = (0:30:360)/180*pi;
    % Note that pigment.pdWidth defines the size of a square collective
    % aperture. Here we compute the equivalent circular aperture
    dx = sqrt((theMosaic.pigment.pdWidth^2)/pi)*2;
    apertureOutline.x = dx/2 * cos(iTheta)*1e6;
    apertureOutline.y = dx/2 * sin(iTheta)*1e6;
    
    edgeColor = 'none';
    lineWidth = 1.0;
    renderPatchArray(apertureOutline, mosaicXaxis, mosaicYaxis, hexMap, edgeColor, lineWidth);
    axis 'image'; axis 'xy'; %box 'on';
    set(gca, 'CLim', signalRange, 'Color', [0 0 0]);
end


function renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                responseNormalization, signalName, instanceLabel, eyeMovementPathsToThisPoint,...
                colorbarTicks, colorbarTickLabels, displayColorbar, xTicks, yTicks, xTickLabels, yTickLabels)      
    if (isHexActivation)
        edgeColor = 'none';
        lineWidth = 1.0;
        renderPatchArray(apertureOutline, mosaicXaxis, mosaicYaxis, activation, edgeColor, lineWidth);
    else
        imagesc(mosaicXaxis, mosaicYaxis, activation);
    end
    hold on;
    instancesVisualized = size(eyeMovementPathsToThisPoint,1);
    instanceColors = jet(instancesVisualized+3);
    for k = 1:instancesVisualized 
        plot(squeeze(eyeMovementPathsToThisPoint(k,:,1)), squeeze(eyeMovementPathsToThisPoint(k,:,2)), '-', 'LineWidth', 1.5, 'Color', squeeze(instanceColors(k,:)));
    end
    hold off
    axis 'image'; axis 'xy'; box 'on';
    set(gca, 'XTick', xTicks, 'YTick', [], 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5]);
    set(gca, 'CLim', [0 1], 'FontSize', 12, 'Color', [0 0 0]);
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
        hCbar = colorbar('Ticks', colorbarTicks, 'TickLabels', sprintf('%2.2f\n',colorbarTickLabels));
        hCbar.Orientation = 'vertical'; 
        hCbar.Label.String = colorbarLabel;
        hCbar.FontSize = 12; 
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