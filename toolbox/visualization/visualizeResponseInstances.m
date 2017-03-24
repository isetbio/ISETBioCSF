function hFigs = visualizeResponseInstances(theMosaic, stimData, noStimData, responseNormalization, condIndex, condsNum, format)

    instancesNum = size(stimData.responseInstanceArray.theMosaicIsomerizations,1);
    if (instancesNum < 1)
        return;
    end
    
    % NEW
    hFigs1 = visualizeNoiseFreeResponses(theMosaic, stimData, noStimData);
    hFigs2 = visualizeSingleConeResponseInstances(theMosaic, stimData, noStimData);

    hFigs = [hFigs1 hFigs2];
    
    % OLD
    %visualizeXYTResponseInstances(theMosaic, stimData, noStimData, responseNormalization, condIndex, condsNum, format);
end

function hFig = visualizeSingleConeResponseInstances(theMosaic, stimData, noStimData)
    
    % transform isomerization counts to isomerization rate
    stimData.responseInstanceArray.theMosaicIsomerizations = stimData.responseInstanceArray.theMosaicIsomerizations / theMosaic.integrationTime;
    noStimData.responseInstanceArray.theMosaicIsomerizations = noStimData.responseInstanceArray.theMosaicIsomerizations / theMosaic.integrationTime;
    stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations / theMosaic.integrationTime;
    noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations / theMosaic.integrationTime;
    
    isomerizationsRange = [...
        min([min(noStimData.responseInstanceArray.theMosaicIsomerizations(:)) min(stimData.responseInstanceArray.theMosaicIsomerizations(:))]) ...
        max([max(noStimData.responseInstanceArray.theMosaicIsomerizations(:)) max(stimData.responseInstanceArray.theMosaicIsomerizations(:))]) ];
    
    photocurrentsRange = [...
        min([min(noStimData.responseInstanceArray.theMosaicPhotocurrents(:)) min(stimData.responseInstanceArray.theMosaicPhotocurrents(:))]) ...
        max([max(noStimData.responseInstanceArray.theMosaicPhotocurrents(:)) max(stimData.responseInstanceArray.theMosaicPhotocurrents(:))]) ];
    
    nonNullCones = theMosaic.pattern(theMosaic.pattern>1);
    for submosaicIndex = 1:3
        % find the cone indices for this submosaic
        submosaicConeIndices = find(nonNullCones==submosaicIndex+1);
        
        % find the current submosaic's cone that has the max noise-free delta (stim - noStim) photocurrent signal
        submosaicNoiseFreeDeltaPhotocurrents = stimData.noiseFreePhotocurrents(submosaicConeIndices,:) - noStimData.noiseFreePhotocurrents(submosaicConeIndices,:);
        
        [maxDeltaPhotocurrent, idx] = max(submosaicNoiseFreeDeltaPhotocurrents(:));
        [coneIndexOfPeakResponse, timeBinOfPeakResponse(submosaicIndex)] = ind2sub(size(submosaicNoiseFreeDeltaPhotocurrents), idx);
        
        submosaicInstances        = noStimData.responseInstanceArray.theMosaicIsomerizations(:,submosaicConeIndices,:);
        submosaicNoiseFreeSignals = noStimData.noiseFreeIsomerizations(submosaicConeIndices,:);
        isomerizationNoStimResponseInstances(submosaicIndex,:,:) = submosaicInstances(:, coneIndexOfPeakResponse, :);
        isomerizationNoStimNoiseFreeResponse(submosaicIndex,:) = submosaicNoiseFreeSignals(coneIndexOfPeakResponse,:);
        
        submosaicInstances        = noStimData.responseInstanceArray.theMosaicPhotocurrents(:,submosaicConeIndices,:);
        submosaicNoiseFreeSignals = noStimData.noiseFreePhotocurrents(submosaicConeIndices,:);
        photocurrentNoStimResponseInstances(submosaicIndex,:,:) = submosaicInstances(:, coneIndexOfPeakResponse, :);
        photocurrentNoStimNoiseFreeResponse(submosaicIndex,:) = submosaicNoiseFreeSignals(coneIndexOfPeakResponse,:);
        
        submosaicInstances        = stimData.responseInstanceArray.theMosaicIsomerizations(:,submosaicConeIndices,:);
        submosaicNoiseFreeSignals = stimData.noiseFreeIsomerizations(submosaicConeIndices,:);
        isomerizationStimResponseInstances(submosaicIndex,:,:) = submosaicInstances(:, coneIndexOfPeakResponse, :);
        isomerizationStimNoiseFreeResponse(submosaicIndex,:) = submosaicNoiseFreeSignals(coneIndexOfPeakResponse,:);
        
        submosaicInstances        = stimData.responseInstanceArray.theMosaicPhotocurrents(:,submosaicConeIndices,:);
        submosaicNoiseFreeSignals = stimData.noiseFreePhotocurrents(submosaicConeIndices,:);
        photocurrentStimResponseInstances(submosaicIndex,:,:) = submosaicInstances(:, coneIndexOfPeakResponse, :);
        photocurrentStimNoiseFreeResponse(submosaicIndex,:) = submosaicNoiseFreeSignals(coneIndexOfPeakResponse,:);
    end
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 2, ...
           'heightMargin',   0.09, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.07, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.04);
       
    nTrials = size(isomerizationNoStimResponseInstances,2);
    opacity = 0.2;
    
    for signalIndex = 1:2
        hFig{signalIndex} = figure(4+signalIndex); clf;
        set(hFig{signalIndex}, 'Position', [10 10 900 950], 'Color', [1 1 1]);
    
        for submosaicIndex = 1:3
            
            % NULL STIMULUS
            subplot('Position', subplotPosVectors(submosaicIndex,1).v);
            hold on
            if (signalIndex == 1)
                for instanceIndex = 1:nTrials
                    scatter(1000*noStimData.responseInstanceArray.timeAxis, squeeze(isomerizationNoStimResponseInstances(submosaicIndex,instanceIndex,:)), 'ks', 'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2], 'LineWidth', 1.0);
                end
                plot(1000*noStimData.responseInstanceArray.timeAxis, squeeze(isomerizationNoStimNoiseFreeResponse(submosaicIndex,:)), 'r-', 'LineWidth', 1.5);
                alpha(opacity);
            else
                for instanceIndex = 1:nTrials
                    h = plot(1000*noStimData.responseInstanceArray.timeAxis, squeeze(photocurrentNoStimResponseInstances(submosaicIndex,instanceIndex,:)), '-', 'Color', [0.2 0.2 0.2],'LineWidth', 1.0);
                    h.Color(4) = opacity;
                end
                plot(1000*noStimData.responseInstanceArray.timeAxis, squeeze(photocurrentNoStimNoiseFreeResponse(submosaicIndex,:)), 'r-', 'LineWidth', 1.5);
            end
            hold off
            
            if (signalIndex == 1)
                set(gca, 'YLim', isomerizationsRange, 'XLim', 1000*[noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)], 'FontSize', 12);
                if (submosaicIndex == 1)
                    title(sprintf('NULL stim (nTrials = %d)\npeak L-cone isomerization responses', nTrials));
                elseif (submosaicIndex == 2)
                    title(sprintf('NULL stim (nTrials = %d)\npeak M-cone isomerization responses', nTrials));
                elseif (submosaicIndex == 3)
                    title(sprintf('NULL stim (nTrials = %d)\npeak S-cone isomerization responses', nTrials));
                end
                ylabel('R*/cone/sec');
            else
                set(gca, 'YLim',photocurrentsRange, 'XLim', 1000*[noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)], 'FontSize', 12);
                if (submosaicIndex == 1)
                    title(sprintf('NULL stim (nTrials = %d)\npeak L-cone photocurrent responses', nTrials));
                elseif (submosaicIndex == 2)
                    title(sprintf('NULL stim (nTrials = %d)\npeak M-cone photocurrent responses', nTrials));
                elseif (submosaicIndex == 3)
                    title(sprintf('NULL stim (nTrials = %d)\npeak S-cone photocurrent responses', nTrials));
                end
                ylabel('pAmps');
            end
            
            if (submosaicIndex == 3)
                xlabel('time (msec)');
            end
        
            % STIMULUS
            subplot('Position', subplotPosVectors(submosaicIndex,2).v);
            hold on
            if (signalIndex == 1)
                for instanceIndex = 1:nTrials
                    scatter(1000*noStimData.responseInstanceArray.timeAxis, squeeze(isomerizationStimResponseInstances(submosaicIndex,instanceIndex,:)), 'ks', 'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2], 'LineWidth', 1.0);
                end
                plot(1000*noStimData.responseInstanceArray.timeAxis, squeeze(isomerizationStimNoiseFreeResponse(submosaicIndex,:)), 'r-', 'LineWidth', 1.5);
                alpha(opacity);
            else
                for instanceIndex = 1:nTrials
                    h = plot(1000*noStimData.responseInstanceArray.timeAxis, squeeze(photocurrentStimResponseInstances(submosaicIndex,instanceIndex,:)), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);
                    h.Color(4) = opacity;
                end
                plot(1000*noStimData.responseInstanceArray.timeAxis, squeeze(photocurrentStimNoiseFreeResponse(submosaicIndex,:)), 'r-', 'LineWidth', 1.5);
            end
            hold off
            
            if (signalIndex == 1)
                set(gca, 'YLim', isomerizationsRange, 'XLim', 1000*[noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)], 'FontSize', 12);
                set(gca, 'YTickLabel', {});
                if (submosaicIndex == 1)
                    title(sprintf('Stim (nTrials = %d)\npeak L-cone isomerization responses', nTrials));
                elseif (submosaicIndex == 2)
                    title(sprintf('Stim (nTrials = %d)\npeak M-cone isomerization responses (R*/cone/sec)', nTrials));
                elseif (submosaicIndex == 3)
                    title(sprintf('Stim (nTrials = %d)\npeak S-cone isomerization responses (R*/cone/sec)', nTrials));
                end
            else
                set(gca, 'YLim', photocurrentsRange, 'XLim', 1000*[noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(end)], 'FontSize', 12);
                set(gca, 'YTickLabel', {});
                if (submosaicIndex == 1)
                    title(sprintf('Stim (nTrials = %d)\npeak L-cone photocurrent responses', nTrials));
                elseif (submosaicIndex == 2)
                    title(sprintf('Stim (nTrials = %d)\npeak M-cone photocurrent responses', nTrials));
                elseif (submosaicIndex == 3)
                    title(sprintf('Stim (nTrials = %d)\npeak S-cone photocurrent responses', nTrials));
                end
            end
            
            if (submosaicIndex == 3)
                xlabel('time (msec)');
            end
        end % submosaicIndex
        drawnow;
    end % signalIndex
end


function hFig = visualizeNoiseFreeResponses(theMosaic, stimData, noStimData)

    % Visualize the os impulse response functions
    hFig{1} = figure(1); clf;
    set(hFig{1}, 'Position', [10 10 700 500], 'Color', [1 1 1]);
    hold on
    plot(1000*stimData.osImpulseResponseTimeAxis, stimData.osImpulseResponses(:,1), 'r-', 'LineWidth', 1.5);
    plot(1000*stimData.osImpulseResponseTimeAxis, stimData.osImpulseResponses(:,2), 'g-', 'LineWidth', 1.5);
    plot(1000*stimData.osImpulseResponseTimeAxis, stimData.osImpulseResponses(:,3), 'b-', 'LineWidth', 1.5);
    xlabel('time (msec)');
    set(gca, 'FontSize', 14);
    grid on;
    title(sprintf('os impulse response functions'));
    drawnow;
    
    % transform isomerization counts to isomerization rate
    stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations / theMosaic.integrationTime;
    noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations / theMosaic.integrationTime;
    
    isomerizationsRange = [...
        min([min(noStimData.noiseFreeIsomerizations(:)) min(stimData.noiseFreeIsomerizations(:))]) ...
        max([max(noStimData.noiseFreeIsomerizations(:)) max(stimData.noiseFreeIsomerizations(:))]) ];
    photocurrentsRange = [...
        min([min(noStimData.noiseFreePhotocurrents(:)) min(stimData.noiseFreePhotocurrents(:))]) ...
        max([max(noStimData.noiseFreePhotocurrents(:)) max(stimData.noiseFreePhotocurrents(:))]) ];
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',   0.1, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.04);
       
    hFig{2} = figure(2); clf;
    set(hFig{2}, 'Position', [10 10 1100 820], 'Color', [1 1 1]);
    subplot('Position', subplotPosVectors(1,1).v);
    imagesc(1000*noStimData.responseInstanceArray.timeAxis, 1:size(noStimData.noiseFreeIsomerizations,1), noStimData.noiseFreeIsomerizations);
    set(gca, 'CLim', isomerizationsRange, 'FontSize', 12);
    title(sprintf('NULL stim\nnoise-free isomerization rate (R*/cone/sec)'));
    ylabel('cone #');
    
    subplot('Position', subplotPosVectors(1,2).v);
    imagesc(1000*stimData.responseInstanceArray.timeAxis, 1:size(stimData.noiseFreeIsomerizations,1), stimData.noiseFreeIsomerizations);
    set(gca, 'CLim', isomerizationsRange, 'YTickLabel', {}, 'FontSize', 12);
    title(sprintf('Stim\nnoise-free isomerization rate (R*/cone/sec)'));
    colorbar
    
    subplot('Position', subplotPosVectors(2,1).v);
    imagesc(1000*noStimData.responseInstanceArray.timeAxis, 1:size(noStimData.noiseFreePhotocurrents,1), noStimData.noiseFreePhotocurrents);
    set(gca, 'CLim', photocurrentsRange, 'FontSize', 12);
    ylabel('cone #');
    xlabel('time (seconds)');
    title(sprintf('NULL stim\nnoise-free photocurrents (pAmps)'));
    
    subplot('Position', subplotPosVectors(2,2).v);
    imagesc(1000*stimData.responseInstanceArray.timeAxis, 1:size(stimData.noiseFreePhotocurrents,1), stimData.noiseFreePhotocurrents);
    set(gca, 'CLim', photocurrentsRange, 'YTickLabel', {}, 'FontSize', 12);
    xlabel('time (seconds)');
    title(sprintf('Stim\nnoise-free photocurrents (pAmps)'));
    colorbar
    colormap(gray(1024));
    drawnow;

    if (isa(theMosaic, 'coneMosaicHex'))
        
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.04);
       
        hFig{3} = figure(3); clf;
        set(hFig{3}, 'Position', [10 10 1340 1290], 'Color', [1 1 1]);
    
        subplot('Position', subplotPosVectors(1,2).v);
        timeBinOfPeakIsomerizationResponse = renderHexActivationMap(theMosaic, stimData.noiseFreeIsomerizations, isomerizationsRange, []);
        title(sprintf('Stim\nnoise-free isomerization rate (R*/cone/sec) [t = %2.2f msec]', 1000*stimData.responseInstanceArray.timeAxis(timeBinOfPeakIsomerizationResponse)), 'FontSize', 14);
        set(gca, 'XTickLabel', {});
        set(gca, 'YTickLabel', {});
        
        subplot('Position', subplotPosVectors(1,1).v);
        timeBinOfPeakIsomerizationResponse = renderHexActivationMap(theMosaic, noStimData.noiseFreeIsomerizations, isomerizationsRange, timeBinOfPeakIsomerizationResponse);
        title(sprintf('NULL stim\nnoise-free isomerization rate (R*/cone/sec) [t = %2.2f msec]', 1000*stimData.responseInstanceArray.timeAxis(timeBinOfPeakIsomerizationResponse)), 'FontSize', 14);
        set(gca, 'XTickLabel', {});
        
        subplot('Position', subplotPosVectors(2,2).v);
        timeBinOfPeakPhotocurrentResponse = renderHexActivationMap(theMosaic, stimData.noiseFreePhotocurrents, photocurrentsRange, []);
        title(sprintf('Stim\nnoise-free photocurrents (pAmps) [t = %2.2f msec]', 1000*stimData.responseInstanceArray.timeAxis(timeBinOfPeakPhotocurrentResponse)), 'FontSize', 14);
        set(gca, 'YTickLabel', {});
        
        subplot('Position', subplotPosVectors(2,1).v);
        timeBinOfPeakPhotocurrentResponse = renderHexActivationMap(theMosaic, noStimData.noiseFreePhotocurrents, photocurrentsRange, timeBinOfPeakPhotocurrentResponse);
        title(sprintf('NULL stim\nnoise-free photocurrents (pAmps) [t = %2.2f msec]', 1000*stimData.responseInstanceArray.timeAxis(timeBinOfPeakPhotocurrentResponse)), 'FontSize', 14);
        
        colormap(gray(1024));
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
    axis 'image'; axis 'xy'; box 'on';
    set(gca, 'CLim', signalRange, 'Color', [0 0 0], 'FontSize', 12);
    colorbar
end

function visualizeXYTResponseInstances(theMosaic, stimData, noStimData, responseNormalization, condIndex, condsNum, format)
    instancesNum = size(stimData.responseInstanceArray.theMosaicIsomerizations,1);
    
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

     
    if (isa(theMosaic, 'coneMosaicHex'))
        for instanceIndex = 1:instancesNum
            if (instanceIndex==1)
                tmp = squeeze(absorptions(instanceIndex,:,:));
                if (timeBins == 1)
                    tmp = tmp';
                end
                tmp = theMosaic.reshapeHex2DmapToHex3Dmap(tmp);
                absorptionsHex = zeros(instancesNum, size(tmp,1), size(tmp,2), size(tmp,3), 'single');
                absorptionsHex(instanceIndex,:,:,:) = tmp;
                
                if (~isempty(photocurrents))
                    tmp = theMosaic.reshapeHex2DmapToHex3Dmap(squeeze(photocurrents(instanceIndex,:,:)));
                    photocurrentsHex = zeros(instancesNum, size(tmp,1), size(tmp,2), size(tmp,3), 'single');
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
        end % instanceIndex
        
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
            set(hFig(1), 'Position', [10 10 2290 1650], 'Color', [0 0 0], 'Name', sprintf('Absorptions (first instance, integrationTime: %2.2fms); %s', theMosaic.integrationTime*1000, stimData.stimulusLabel));

            if (~isempty(photocurrents))
                randomBaseFigNum = randomBaseFigNum+1000;
                hFig(2) = figure(randomBaseFigNum+condIndex); clf;
                set(hFig(2), 'Position', [10+200 10+20 2290 1650], 'Color', [0 0 0], 'Name', sprintf('Photocurrents (first instance); %s', stimData.stimulusLabel));
            end

            subplotCols = max([1 floor(1.25*sqrt(numel(absorptionsTimeAxis)))]);
            subplotRows = ceil(numel(absorptionsTimeAxis) / subplotCols);
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', subplotRows, ...
                   'colsNum', subplotCols, ...
                   'heightMargin',   0.02, ...
                   'widthMargin',    0.01, ...
                   'leftMargin',     0.01, ...
                   'rightMargin',    0.01, ...
                   'bottomMargin',   0.01, ...
                   'topMargin',      0.01);

    for k = 1:numel(hFig)
        figure(hFig(k));

        for tBin = 1:numel(absorptionsTimeAxis)
            row = floor((tBin-1)/subplotCols) + 1;
            col = mod(tBin-1, subplotCols)+1;
            subplot('Position', subplotPosVectors(row,col).v);
            if (k == 1)
                activation = squeeze(absorptions(instanceIndex, :,:,tBin));
                instanceLabel = sprintf('%2.1fms (i:1, em:1-%d)', absorptionsTimeAxis(tBin)*1000, instancesNum);
                minActivation = minAbsorptions;
                maxActivation = maxAbsorptions;
            else
                activation = squeeze(photocurrents(instanceIndex, :,:,tBin));
                instanceLabel = sprintf('%2.1fms (i:1, em:1-%d)', photocurrentsTimeAxis(tBin)*1000, instancesNum);
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
            eyeMovementPathsToThisPoint = stimData.responseInstanceArray.theMosaicEyeMovements(:,1:tBin,:);
            renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                        responseNormalization, signalName, instanceLabel, ...
                        eyeMovementPathsToThisPoint, ...
                        colorbarTicks, minActivation+colorbarTicks*(maxActivation-minActivation), displayColorBar,  ...
                        xTicks, yTicks, xTickLabels, yTickLabels);
            drawnow;
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
           'rightMargin',    0.04, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
    
        for instanceIndex = 1:instancesNum
         for tBin = 1: numel(absorptionsTimeAxis)  
             
            % Instance absorptions on the left
            subplot('Position', subplotPosVectors(1,1).v);
            activation = squeeze(absorptions(instanceIndex, :,:,tBin));
            if (isHexActivation)
                activation = activation(activeConesActivations);
            end
            eyeMovementPathsToThisPointForThisInstance = stimData.responseInstanceArray.theMosaicEyeMovements(instanceIndex,1:tBin,:);
            renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                responseNormalization, sprintf('absorptions (intTime: %2.2fms)',  theMosaic.integrationTime*1000), ...
                sprintf('instance %d/%d (t: %2.1fms)', instanceIndex, instancesNum, absorptionsTimeAxis(tBin)*1000),  ...
                eyeMovementPathsToThisPointForThisInstance, ...
                colorbarTicks, minAbsorptions + colorbarTicks*(maxAbsorptions-minAbsorptions), true, ...
                xTicks, yTicks, xTickLabels, yTickLabels);
           
            % Noise-free isomerizations on the right
            subplot('Position', subplotPosVectors(1,2).v);
            activation = squeeze(noiseFreeIsomerizations(:,:,tBin));
            if (isHexActivation)
                activation = activation(activeConesActivations);
            end
            eyeMovementPathsToThisPointForAllInstances = stimData.responseInstanceArray.theMosaicEyeMovements(:,1:tBin,:);
            renderPlot(isHexActivation, mosaicXaxis, mosaicYaxis, activation, apertureOutline, ...
                responseNormalization, sprintf('noise-free absorptions (intTime: %2.2fms)',  theMosaic.integrationTime*1000), ...
                sprintf('cond: %d/%d (t: %2.1fms)', condIndex, condsNum, absorptionsTimeAxis(tBin)*1000), ...
                eyeMovementPathsToThisPointForAllInstances, ...
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
                    eyeMovementPathsToThisPointForThisInstance, ...
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
                    eyeMovementPathsToThisPointForAllInstances, ...
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