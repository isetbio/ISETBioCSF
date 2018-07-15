function hFigsInfo = visualizeResponseInstances(theMosaic, ...
    stimData, noStimData, visualizeOuterSegmentFilters, ...
    responseNormalization, condIndex, condsNum, format)

    % Which best cone types to visualize instances from
    coneTypeIndicesToVisualize = [1 2 3];  % the best L, M- and S-cone
    coneTypeIndicesToVisualize = [1];      % only the L-cone
    
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
   
    % Noise free mosaic responses XT and XYat the time bin of peak response
    signalsToVisualize = 'isomerizationsAndPhotocurrents';
    signalsToVisualize = 'isomerizationsOnly';
    
    [hFigs, peakConeIndex] = visualizeNoiseFreeXTandXYResponses(theMosaic, timeAxis, stimData, noStimData, 6000, signalsToVisualize);
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXTResponses',...
            'hFig', hFigs{1});
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXYResponses',...
            'hFig', hFigs{2});
        
        
    hFigs = visualizeNoiseFreeAndFirstInstanceWithActivationProfile(...
             stimData, noStimData, theMosaic);   
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXYResponsesFirstInstanceAndActivationProfilesSTIM',...
            'hFig', hFigs{1});   
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXYResponsesFirstInstanceAndActivationProfilesNULL',...
            'hFig', hFigs{2});
        
    % Single L, M, and S-cone (best responses) showing distribution of all response instances
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
    
    if (~isempty(noStimData.responseInstanceArray.theMosaicIsomerizations))
        isomerizationLevels = [0 30];
        isomerizationQuantizationLevelsNum = round(isomerizationLevels(2)-isomerizationLevels(1));
        hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
            timeAxis, ...
            noStimIsomerizationsResponseInstances*theMosaic.integrationTime, ...
            stimIsomerizationsResponseInstances*theMosaic.integrationTime, ...
            noStimNoiseFreeIsomerizationsResponse*theMosaic.integrationTime, ...
            stimNoiseFreeIsomerizationsResponse*theMosaic.integrationTime, ...
            isomerizationLevels, ...
            isomerizationQuantizationLevelsNum, ...
            coneTypeIndicesToVisualize, ...
            'density', 'isomerizations', 'R*/cone', 5001);
    else
        hFig = [];
    end 
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'peakConeIsomerizationResponseInstances',...
            'hFig', hFig);
    
        

    if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
        photocurrentsQuantizationLevelsNum = 100;
        photocurrentsLevels = [-90 -50]; 
        hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
            timeAxis, ...
            noStimPhotocurrentsResponseInstances, ...
            stimPhotocurrentsResponseInstances, ...
            noStimNoiseFreePhotocurrentsResponse, ...
            stimNoiseFreePhotocurrentsResponse, ...
            photocurrentsLevels, ...
            photocurrentsQuantizationLevelsNum, ...
            coneTypeIndicesToVisualize, ...
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
    

function [hFig, peakConeIndex] = visualizeNoiseFreeXTandXYResponses(theMosaic, timeAxis, stimData, noStimData, figNo, signals)

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

    % Round to nearest 200 R*/cone/sec
    isomerizationsRange = round(isomerizationsRange*200)/200;
    
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
    
    
    if (strcmp(signals, 'isomerizationsAndPhotocurrents'))
        rowsNum = 2;
    else
        rowsNum = 1;
    end

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum, ...
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
            
            differentialStimNoStimIsomerizationResponses = (stimData.noiseFreeIsomerizations(submosaicConeIndices,:) - noStimData.noiseFreeIsomerizations(submosaicConeIndices,:));
            differentialStimNoStimIsomerizationResponses = 100 * differentialStimNoStimIsomerizationResponses ./ noStimData.noiseFreeIsomerizations(submosaicConeIndices,:);
            [bestModulation, idx] = max(abs(differentialStimNoStimIsomerizationResponses(:)));
            [bestCone, bestTimeBin] = ind2sub(size(differentialStimNoStimIsomerizationResponses), idx);
            peakConeIndex{submosaicIndex-1} = submosaicConeIndices(bestCone);
        end
    end


    isomerizationModulation = 100*(stimLMSnoiseFreeIsomerizations - noStimLMSnoiseFreeIsomerizations)./noStimLMSnoiseFreeIsomerizations;
    isomerizationsRangeModulation = max(abs(isomerizationModulation(:)))*[-1 1];
    
    hFig{1} = figure(figNo+1); clf;
    if (rowsNum == 2)
        set(hFig{1}, 'Position', [10 10 770 900], 'Color', [1 1 1]);
    else
        set(hFig{1}, 'Position', [10 10 770 450], 'Color', [1 1 1]);
    end
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
    
    if (~isempty(noStimData.noiseFreePhotocurrents)) && (rowsNum == 2)
        deltaPhotocurrents = stimLMSnoiseFreePhotocurrents - noStimLMSnoiseFreePhotocurrents;
        [~,idx] = max(abs(deltaPhotocurrents(:)));
        [~,timeBinOfPeakPhotocurrentResponse] = ind2sub(size(deltaPhotocurrents), idx);
        deltaPhotocurrentsRange = [min(deltaPhotocurrents(:,timeBinOfPeakPhotocurrentResponse)) max(deltaPhotocurrents(:,timeBinOfPeakPhotocurrentResponse))];
        
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
    
    
    % Plot the 2D XY noise-free responses and modulations
    if (isa(theMosaic, 'coneMosaicHex'))
       
        stimData.noiseFreeIsomerizationsModulations = 100*(stimData.noiseFreeIsomerizations-noStimData.noiseFreeIsomerizations)./noStimData.noiseFreeIsomerizations;
        isomerizationsRangeModulation = max(abs(stimData.noiseFreeIsomerizationsModulations(:)))*[-1 1];
    
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum,  ...
           'colsNum', 2, ...
           'heightMargin',   0.05, ...
           'widthMargin',    0.07, ...
           'leftMargin',     0.015, ...
           'rightMargin',    0.05, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.04);
       
        hFig{2} = figure(figNo+2); clf;
        if (~isempty(noStimData.noiseFreePhotocurrents)) && (rowsNum == 2)
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_TWO_CONDIITIONS');
        else
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_SINGLE_CONDIITION');
        end
        
        [~,idx] = max(stimData.noiseFreeIsomerizationsModulations(:));
        [bestRespondingConeIndex,timeBinOfPeakIsomerizationResponse] = ind2sub(size(stimData.noiseFreeIsomerizationsModulations), idx);
        visualizedActivationPattern = theMosaic.reshapeHex1DmapToHex2Dmap(squeeze(stimData.noiseFreeIsomerizationsModulations(:,timeBinOfPeakIsomerizationResponse)));
        
        % The response modulations for the test stimulus on the right
        activationLUT = brewermap(512, '*RdBu');
        ax1 = plotMosaicActivation(theMosaic, visualizedActivationPattern, ...
             isomerizationsRangeModulation, activationLUT, [1 1 1], subplotPosVectors, [1 2], ...
             sprintf('test stimulus\nmean isomerization modulation'), ...
             sprintf('modulation (%%)'), ...
             'displaylabelX', (rowsNum == 1),  'displayTicksX', (rowsNum == 1), ...
             'displaylabelY', false, 'displayTicksY', false);
        
        % The response for the null stimulus on the left
        visualizedActivationPattern = theMosaic.reshapeHex1DmapToHex2Dmap(squeeze(noStimData.noiseFreeIsomerizations(:,timeBinOfPeakIsomerizationResponse)));
        activationLUT = gray(1024);
        ax2 = plotMosaicActivation(theMosaic, visualizedActivationPattern*theMosaic.integrationTime, ...
             isomerizationsRange*theMosaic.integrationTime, activationLUT, [0 0 0], subplotPosVectors, [1 1], ...
             sprintf('null stimulus\n mean isomerizations / %2.0f msec', 1000*theMosaic.integrationTime), ...
             sprintf('R* / cone /%2.0fmsec', 1000*theMosaic.integrationTime), ...
             'displaylabelX', (rowsNum == 1), 'displayTicksX', (rowsNum == 1), ...
             'displaylabelY', false, 'displayTicksY', false);
        
        if (~isempty(noStimData.noiseFreePhotocurrents)) && (rowsNum == 2)
            deltaPhotocurrents = (stimData.noiseFreePhotocurrents - noStimData.noiseFreePhotocurrents);
            [~,idx] = max(abs(deltaPhotocurrents(:)));
            [~,timeBinOfPeakPhotocurrentResponse] = ind2sub(size(deltaPhotocurrents), idx);
            deltaPhotocurrentsRange = max(abs(deltaPhotocurrents(:)))*[-1 1];
            visualizedActivationPattern = theMosaic.reshapeHex1DmapToHex2Dmap(squeeze(deltaPhotocurrents(:,timeBinOfPeakPhotocurrentResponse)));

            activationLUT = brewermap(512, '*RdBu');
            ax3 = plotMosaicActivation(theMosaic, visualizedActivationPattern, ...
                deltaPhotocurrentsRange, activationLUT, [1 1 1], subplotPosVectors, [2 2], ...
                sprintf('test stimulus\n mean photocurrent differential'), ...
                sprintf('pAmps'), ...
                'displaylabelY', false, 'displayTicksY', false);

        
            visualizedActivationPattern = theMosaic.reshapeHex1DmapToHex2Dmap(squeeze(noStimData.noiseFreePhotocurrents(:,timeBinOfPeakPhotocurrentResponse)));
            activationLUT = gray(1024);
            ax4 = plotMosaicActivation(theMosaic, visualizedActivationPattern, ...
                photocurrentsRange, activationLUT, [0 0 0], subplotPosVectors, [2 1], ...
                sprintf('null stimulus\n mean photocurrent'), ...
                sprintf('pAmps'), ...
                'displaylabelY', false, 'displayTicksY', false);
            
        end
        
        if (~isempty(noStimData.noiseFreePhotocurrents)) && (rowsNum == 2)
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_TWO_CONDIITIONS', 'theAxes', ax1);
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_TWO_CONDIITIONS', 'theAxes',ax2);
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_TWO_CONDIITIONS', 'theAxes',ax3);
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_TWO_CONDIITIONS', 'theAxes',ax4);
        else
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_SINGLE_CONDIITION', 'theAxes',ax1);
            formatFigureForPaper(hFig{2}, 'figureType', 'RESPONSE_INSTANCE_SINGLE_CONDIITION', 'theAxes',ax2);
        end
            
        drawnow;
    end
end


function ax = plotMosaicActivation(cMosaic, visualizedActivationPattern, ...
            signalRange, activationLUT, backgroundColor, subplotPosVectors, subplotPos, ...
            titleText, titleForColorBar, varargin)

    p = inputParser;
    p.addParameter('displaylabelX',true, @islogical);   
    p.addParameter('displaylabelY',true, @islogical); 
    p.addParameter('displayTicksX',true, @islogical); 
    p.addParameter('displayTicksY',true, @islogical); 
    p.parse(varargin{:});

    ax = subplot('Position', subplotPosVectors(subplotPos(1),subplotPos(2)).v); 
    cMosaic.renderActivationMap(ax, visualizedActivationPattern, ...
             'signalRange', signalRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', true, ...
             'labelColorBarTicks', true, ...
             'titleForColorBar', titleForColorBar, ...
             'colorMap', activationLUT, ...
             'backgroundColor', backgroundColor);
    
    if (~p.Results.displaylabelX); xlabel(ax, ''); end
    if (~p.Results.displaylabelY); ylabel(ax, ''); end
    if (~p.Results.displayTicksX); set(ax, 'XTickLabels', {}); end
    if (~p.Results.displayTicksY); set(ax, 'YTickLabels', {}); end
    title(ax, titleText);
end

function hFigs = visualizeNoiseFreeAndFirstInstanceWithActivationProfile(...
             stimData, noStimData, theMosaic)
    
    % Transform rates to counts
    stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations * theMosaic.integrationTime;
    noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations * theMosaic.integrationTime;
    stimData.responseInstanceArray.theMosaicIsomerizations = stimData.responseInstanceArray.theMosaicIsomerizations * theMosaic.integrationTime;
    noStimData.responseInstanceArray.theMosaicIsomerizations = noStimData.responseInstanceArray.theMosaicIsomerizations * theMosaic.integrationTime;
    
    % Determine visualized time bin of max response
    [~,kidx] = max(stimData.noiseFreeIsomerizations(:));
    [~,timeBinOfPeakIsomerizationResponse] = ind2sub(size(stimData.noiseFreeIsomerizations), kidx);
    
    % Visualized response instances
    visualizedResponseInstances = (squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(:,:,timeBinOfPeakIsomerizationResponse)))';
    % Visualized noise-free isomerizations
    visualizedNoiseFreeNoStimIsomerizations = squeeze(noStimData.noiseFreeIsomerizations(:,timeBinOfPeakIsomerizationResponse));
    
    % Determine visualized response range
    isomerizationsRange = determineVisualizedIsomerizationsRange(theMosaic, stimData, timeBinOfPeakIsomerizationResponse);
    
    % Instance to visualize
    visualizedResponseInstance = 1;
    
    signalRange = [-12 12]; coneLinePlotType = 'differential_activations';
    hFigs{1} = renderFigure(1222, theMosaic, visualizedResponseInstance, ...
        stimData, noStimData, timeBinOfPeakIsomerizationResponse, ...
        isomerizationsRange, signalRange, coneLinePlotType, 'diff cone excitation');    
    
    signalRange = isomerizationsRange; coneLinePlotType = 'activations';
    hFigs{2} = renderFigure(1333, theMosaic, visualizedResponseInstance, ...
        noStimData, noStimData, timeBinOfPeakIsomerizationResponse, ...
        isomerizationsRange, signalRange, coneLinePlotType, 'cone excitation');    
    
end

function hFig = renderFigure(figNo, theMosaic, visualizedResponseInstance, ...
    stimData, noStimData, timeBinOfPeakIsomerizationResponse, ...
    isomerizationsRange, signalRange, coneLinePlotType, yLabelTitle)
    
    % Generate colormaps for modulations and for excitations
    %modulationsColorMap = brewermap(512, '*RdBu');
    excitationsColorMap = gray(1024);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 730 500], 'Color', [1 1 1]);
    ax = subplot('Position', [0.08 0.335 0.44 0.68]);
    theMosaic.renderActivationMap(ax, squeeze(stimData.noiseFreeIsomerizations(:,timeBinOfPeakIsomerizationResponse)) ,...
        'visualizedConeAperture', 'geometricArea', ...
        'signalRange', isomerizationsRange, ...
        'colorMap', excitationsColorMap, ...
        'showColorBar', ~true, ...
        'labelColorBarTicks', ~true, ...
        'outlineConesAlongHorizontalMeridian', true, ...
        'showXLabel', false, ...
        'showYLabel', false, ...
        'backgroundColor', 0*[0.5 0.5 0.5]);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    
    
    ax = subplot('Position', [0.55 0.335 0.44 0.68]);
    theMosaic.renderActivationMap(ax, squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(visualizedResponseInstance,:,timeBinOfPeakIsomerizationResponse))', ...
            'visualizedConeAperture', 'geometricArea', ...
            'signalRange', isomerizationsRange, ...
            'colorMap', excitationsColorMap, ...
            'showColorBar', ~true, ...
            'labelColorBarTicks', ~true, ...
            'outlineConesAlongHorizontalMeridian', true, ...
            'showXLabel', false, ...
            'showYLabel', false, ...
            'backgroundColor', 0*[0.5 0.5 0.5]);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
     
    ax = subplot('Position', [0.08 0.095 0.44 0.24]);
    generateConeLinePlot(ax, theMosaic, ...
        squeeze(stimData.noiseFreeIsomerizations(:,timeBinOfPeakIsomerizationResponse)), ...
        squeeze(noStimData.noiseFreeIsomerizations(:,timeBinOfPeakIsomerizationResponse)), ...
        signalRange, coneLinePlotType, yLabelTitle);
    
    ax = subplot('Position', [0.55 0.095 0.44 0.24]);
    generateConeLinePlot(ax, theMosaic, ...
        squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(visualizedResponseInstance,:,timeBinOfPeakIsomerizationResponse))', ...
        squeeze(noStimData.noiseFreeIsomerizations(:,timeBinOfPeakIsomerizationResponse)), ...
        signalRange, coneLinePlotType, '');
end


function generateConeLinePlot(ax, theMosaic, activation, nullStimActivation, signalRange, coneLinePlotType, yLabelTitle)
    if (any(size(activation) ~= size(theMosaic.pattern)))    
       activation = theMosaic.reshapeHex2DmapToHex3Dmap(activation);
       nullStimActivation = theMosaic.reshapeHex2DmapToHex3Dmap(nullStimActivation);
    end
    sampledHexMosaicXaxis = squeeze(theMosaic.patternSupport(1, :, 1)) + ...
        theMosaic.center(1);
    sampledHexMosaicYaxis = squeeze(theMosaic.patternSupport(:, 1, 2)) + ...
        theMosaic.center(2);
    
    dx = diameterForCircularApertureFromWidthForSquareAperture(...
            theMosaic.pigment.width) * 1e6 / theMosaic.micronsPerDegree;
        
    idx = find(theMosaic.pattern > 1);
    [iRows, iCols] = ind2sub(size(theMosaic.pattern), idx);  
    coneXcoords = (sampledHexMosaicXaxis(iCols))';
    coneYcoords = sampledHexMosaicYaxis(iRows);
    coneActivations = activation(idx);
    nullStimActivation = nullStimActivation(idx);
    coneXcoordsDegs = coneXcoords * 1e6 / theMosaic.micronsPerDegree;
    coneYcoordsDegs = coneYcoords * 1e6 / theMosaic.micronsPerDegree;
    conePosRange = max([max(abs(coneXcoordsDegs)) max(abs(coneYcoordsDegs))]);
    
    % Find cones lying near the y=0 axis
    indicesOfConesAlongXaxis = find(abs(coneYcoordsDegs) < dx);
    coneXcoordsDegs = coneXcoordsDegs(indicesOfConesAlongXaxis);
    coneYcoordsDegs = coneYcoordsDegs(indicesOfConesAlongXaxis);
    coneActivations = coneActivations(indicesOfConesAlongXaxis);
    nullStimActivation  = nullStimActivation(indicesOfConesAlongXaxis);
    identitiesOfConesAlongXaxis = theMosaic.pattern(idx(indicesOfConesAlongXaxis));

    hold(ax, 'on');
    for coneIndex = 2:4
        switch(coneIndex)
            case 2
                edgeColor = [1 0 0];
                faceColor = [1 0.5 0.5];
            case 3
                edgeColor = [0 0.8 0];
                faceColor = [0.5 1 0.5];
            case 4
                edgeColor = [0 0 1];
                faceColor = [0.5 0.8 1];
        end
        iidx = find(identitiesOfConesAlongXaxis == coneIndex);
        if (strcmp(coneLinePlotType, 'differential_activations'))
            signal = coneActivations(iidx)-nullStimActivation(iidx);
            yticks = -10:5:10;
        else
            signal = coneActivations(iidx);
            yticks = 0:10:30;
        end
        plot(ax, coneXcoordsDegs(iidx), signal, ...
            'o', 'Color', edgeColor, 'MarkerSize', 8, 'LineWidth', 1.5, ...
            'MarkerFaceColor',faceColor,'MarkerEdgeColor',edgeColor);
    end
    xtickLabels = {'-0.2', '', '-0.1', '', '0', '', '+0.1', '', '+0.2'};
    set(ax, 'XLim', conePosRange*[-1 1], 'YLim', signalRange, ...
        'FontSize', 16, 'XTick', -0.2:0.05:0.2, 'XTickLabel', xtickLabels, ...
        'YTick', yticks);
    if (~isempty(yLabelTitle))
        ylabel(yLabelTitle, 'FontWeight', 'bold')
    else
        set(gca, 'YTickLabel', {});
    end
    grid 'on'; box 'on';
    xlabel('space (deg)', 'FontWeight', 'bold');
end

function isomerizationsRange = determineVisualizedIsomerizationsRange(theMosaic, stimData, timeBinOfPeakIsomerizationResponse)
    % Determine response range based on the responses of L/M cones
    nonNullCones = theMosaic.pattern(theMosaic.pattern > 1);
    lmConeIndices = find(nonNullCones==2 | nonNullCones==3);
    
    noiseFreeIsomerizationsRange = [...
        min(stimData.noiseFreeIsomerizations(lmConeIndices,timeBinOfPeakIsomerizationResponse)) ...
        max(stimData.noiseFreeIsomerizations(lmConeIndices,timeBinOfPeakIsomerizationResponse)) ...
        ];
    instancesIsomerizationsRange = [ ...
        min(min(squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(:,lmConeIndices,timeBinOfPeakIsomerizationResponse)))) ...
        max(max(squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(:,lmConeIndices,timeBinOfPeakIsomerizationResponse)))) ...
    ];

    isomerizationsRange = noiseFreeIsomerizationsRange;
    
    % add to the range so the visualized range matches the low-end of the 
    % instancesIsomerizationsRange
    delta = isomerizationsRange(1)-instancesIsomerizationsRange(1);
    if (delta > 0)
        isomerizationsRange = isomerizationsRange + delta * [-1 1];
    end
    isomerizationsRange = round(isomerizationsRange);
end
