function hFigsInfo = visualizeResponseInstances(theMosaic, ...
    stimData, noStimData, visualizeOuterSegmentFilters, ...
    condIndex, condsNum)
    
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
    
    % Visualize the mosaic and some EM paths
    hFig = visualizeMosaicAndSomeEMpaths(theMosaic, stimData);
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'mosaicAndSomeEMpaths',...
            'hFig', hFig);
        
        
    % Visualize mosaic activations over time
    hFig = visualizeMosaicActivationsOverTime(theMosaic, stimData);
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'mosaicActivationsMovieFrames',...
            'hFig', hFig);
        
        
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
    %signalsToVisualize = 'isomerizationsOnly';
    
    [hFigs, peakConeIndex] = visualizeNoiseFreeXTandXYResponses(...
        theMosaic, timeAxis, stimData, noStimData, 6000, signalsToVisualize);
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXTResponses',...
            'hFig', hFigs{1});
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'noiseFreeXYResponses',...
            'hFig', hFigs{2});
        
           
    % Response instace visualization
    
    if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
        signalsToVisualizeInstances = {'Isomerization', 'Photocurrent'};
    else
        signalsToVisualizeInstances = {'Isomerization'};
    end
    
    % Instance to visualize
    visualizedResponseInstance = 1;
        
    for signalIndex = 1:numel(signalsToVisualizeInstances)
        signalName = signalsToVisualizeInstances{signalIndex};
        
        if (strcmp(signalName, 'Isomerization'))
            visualizedSignalRange = [-12 12];
        elseif (strcmp(signalName, 'Photocurrent'))
            visualizedSignalRange = [-9 9];
        end

        if (condIndex == condsNum)
            makeVideos = true;
        else
            makeVideos = false;
        end
        
        makeVideos = false;
        
        hFigs = visualizeNoiseFreeAndFirstInstanceWithActivationProfile(...
                 stimData, noStimData, theMosaic, signalName, ...
                 visualizedSignalRange, visualizedResponseInstance, ...
                 makeVideos);   

        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('noiseFreeXY%sResponseSTIM', signalName),...
                'hFig', hFigs{1});  
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('firstInstanceXY%sResponseSTIM',signalName),...
                'hFig', hFigs{2}); 
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('noiseFree%sActivationProfileSTIM',signalName),...
                'hFig', hFigs{3});  
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('firstInstance%sActivationProfileSTIM',signalName),...
                'hFig', hFigs{4});

        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('firstInstanceXTActivation%sProfileSTIM',signalName),...
                'hFig', hFigs{5});


        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('noiseFreeXY%sResponseNULL',signalName),...
                'hFig', hFigs{6});  
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('firstInstanceXY%sResponseNULL',signalName),...
                'hFig', hFigs{7}); 
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('noiseFree%sActivationProfileNULL',signalName),...
                'hFig', hFigs{8});  
        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('firstInstance%sActivationProfileNULL',signalName),...
                'hFig', hFigs{9});

        hFigsInfo{numel(hFigsInfo)+1} = struct(...
                'filename', sprintf('firstInstanceXT%sActivationProfileNULL',signalName),...
                'hFig', hFigs{10});
            
    end % signalIndex

    
    % Single L, M, and S-cone (best responses) showing distribution of all response instances
    %
    % Range of instances to be outlined, in percent
    instancesPercentRange = [10 90];
    
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
        isomerizationResponsRange = [0 30];
        hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
            timeAxis, ...
            noStimIsomerizationsResponseInstances*theMosaic.integrationTime, ...
            stimIsomerizationsResponseInstances*theMosaic.integrationTime, ...
            noStimNoiseFreeIsomerizationsResponse*theMosaic.integrationTime, ...
            stimNoiseFreeIsomerizationsResponse*theMosaic.integrationTime, ...
            isomerizationResponsRange, ...
            instancesPercentRange, ...
            'isomerizations', 'cone excitation', 5001);
    else
        hFig = [];
    end 
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'peakConeIsomerizationResponseInstances',...
            'hFig', hFig);

    if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))  
        photocurrentResponsRange = [-90 -65]; 
        hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
            timeAxis, ...
            noStimPhotocurrentsResponseInstances, ...
            stimPhotocurrentsResponseInstances, ...
            noStimNoiseFreePhotocurrentsResponse, ...
            stimNoiseFreePhotocurrentsResponse, ...
            photocurrentResponsRange, ...
            instancesPercentRange, ...
            'photocurrents', 'p-current (pA)', 5002);
    else
        hFig = [];
    end
    
    hFigsInfo{numel(hFigsInfo)+1} = struct(...
            'filename', 'peakConePhotocurrentResponseInstances',...
            'hFig', hFig);
end


function hFig = visualizeMosaicActivationsOverTime(theMosaic, stimData)

    showDots = false;
    spaceLimsDegs = 5/60*[-1 1];
    spaceLimsMeters = spaceLimsDegs * theMosaic.micronsPerDegree * 1e-6;
    timeLims = [0 200];
    isomerizationsLimits = [0 30];
    photocurrentLimits = [-83 -68];
    
    targetPosDegs = [0 0];
    [targetConeIndexInFullArray, conePositions, coneTypes, targetConeIndex] = indicesForConesAtPositions(theMosaic, targetPosDegs);
    
    
    timeAxis = 1000*stimData.responseInstanceArray.timeAxis;
    isomerizationInstances = stimData.responseInstanceArray.theMosaicIsomerizations;
    meanIsomerizations = stimData.noiseFreeIsomerizations;
    photocurrentInstances = stimData.responseInstanceArray.theMosaicPhotocurrents;
    meanPhotocurrents = stimData.noiseFreePhotocurrents;
    
    
    visualizedResponseInstance = 2; % 13; 24
    meanIsomerizations = squeeze(meanIsomerizations(targetConeIndex,:));
    
    meanPhotocurrent = squeeze(meanPhotocurrents(targetConeIndex,:));
    isomerizationInstanceMean = squeeze(mean(isomerizationInstances(:, targetConeIndex,:),1));
    photocurrentInstanceMean = squeeze(mean(photocurrentInstances(:, targetConeIndex,:),1));
    isomerizationInstance = squeeze(isomerizationInstances(visualizedResponseInstance, targetConeIndex,:));
    photocurrentInstance = squeeze(photocurrentInstances(visualizedResponseInstance, targetConeIndex,:));
    
    theEMpathMicrons = squeeze(stimData.responseInstanceArray.theMosaicEyeMovementsMicrons(visualizedResponseInstance,:,:));
    theEMpathMeters = theEMpathMicrons / 1e6;
    theEMpathDegs = theEMpathMicrons/theMosaic.micronsPerDegree;
    
    
    retinalImage = stimData.thePeakOI.RGBimage.^0.5;
    retinalImage = bsxfun(@minus, retinalImage, retinalImage(1,1,:));
    retinalImage = squeeze(retinalImage(:,:,2));
    m1 = min(retinalImage(:));
    m2 = max(retinalImage(:));
    retinalImage = (retinalImage - m1)/(m2-m1);
    retinalImage = stimData.thePeakOI.RGBimage;
    
    emPathColor = [1 0 0];
    cMap = brewermap(8, 'Reds');
    emPathColor0 = cMap(1,:);
    emPathColor1 = cMap(2,:);
    emPathColor2 = cMap(3,:);
    emPathColor3 = cMap(4,:);
    
    visualizationTimes = [35 40 45 50];
    hFig = figure(345); clf;
    set(hFig, 'Position', [10 10 950  700], 'Color', [1 1 1]);
    
    % Plot the optical image
    width = 0.2*0.78;
    height = 0.26*0.78;
    dX = width*1.005;
    dY = 0
    for k = 1:numel(visualizationTimes)
        ax = axes('Position',[0.03+(k-1)*dX 0.9-height width height]);
        idx(k) = opticalImageEyePathComboPlot(ax, showDots, timeAxis, visualizationTimes(k), theEMpathDegs, emPathColor0,...
            stimData.thePeakOI.xAxisDegs, stimData.thePeakOI.yAxisDegs, retinalImage, spaceLimsDegs);
    end
    

    signalRange = [min(isomerizationInstances(:)) max(isomerizationInstances(:))];
    for k = 1:numel(visualizationTimes)
        ax = axes('Position',[0.03+(k-1)*dX 0.62-height width height]);
        showActivationMap(ax, theMosaic, timeAxis, visualizationTimes(k),  ...
            squeeze(isomerizationInstances(visualizedResponseInstance, :,:)), spaceLimsDegs, spaceLimsMeters,signalRange);
    end   

    colormap(gray(256));
    
%     subplot(2,4,5);
%     plot(timeAxis, theEMpathDegs(:,1), 'r-'); hold on;
%     plot(timeAxis, theEMpathDegs(:,2), 'b-');
%     plot(timeAxis(idx0), theEMpathDegs(idx0,1), '*', 'MarkerSize', 12);
%     plot(timeAxis(idx1), theEMpathDegs(idx1,1), '*', 'MarkerSize', 12);
%     plot(timeAxis(idx2), theEMpathDegs(idx2,1), '*', 'MarkerSize', 12);
%     plot(timeAxis(idx3), theEMpathDegs(idx3,1), '*', 'MarkerSize', 12);
%     plot(timeAxis, timeAxis*0, 'k-');
%     set(gca, 'XLim', [0 150], 'YLim', spaceLimsDegs);
    
    
    ax = axes('Position',[0.8 0.78 width height]);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'labelConeTypes', true,...
        'ticksInVisualDegs', true, ...
        'noXaxisLabel',  true, ...
        'noYaxisLabel', true, ...
        'backgroundColor', [0 0 0] ...
    );
    set(ax, 'XTickLabel', {}, 'YTickLabel',  {}, 'XLim', spaceLimsMeters, 'YLim', spaceLimsMeters);
    
    
    subplot(3,4,10)
    nTimeSamples = 1000;
    nConesNum = 8192;
    deltaT = (timeAxis(2)-timeAxis(1))/1000;
    samplingFrequency = 1/deltaT;
    noiseOnlyResponses = osAddNoise(zeros(1, nConesNum, nTimeSamples), 'sampTime', deltaT);
    noiseOnlyResponses = (squeeze(noiseOnlyResponses(1,:,:)))';
    
    showNoiseSpectrum = ~true;
    if (showNoiseSpectrum)
        [noisePS, freq] = pspectrum(noiseOnlyResponses, samplingFrequency);
    
        noisePS = mean(noisePS,2);
        plot(freq, noisePS, 'k-', 'LineWidth', 1.5);
    
        set(gca, 'FontSize', 14, 'XScale', 'log', 'YScale', 'log', 'XLim', [1 100],  'YTick', [0.01 0.03 0.10], 'YTickLabel', {'.01', '.03', '.10'}, 'XTick', [1 3 10 30 100 300 1000]);
        xlabel('\it frequency (Hz)')
        ylabel('\it pAmps^2 / Hz');
    else
        timeAxis2 = 1000*(1: nTimeSamples)*deltaT-deltaT/2;
        hPlot = plot(timeAxis2, noiseOnlyResponses(:,1:512), 'k-', 'LineWidth', 1.5);
        for k = 1:numel(hPlot)
            hPlot(k).Color(4) = 0.2;  % 5% transparent
        end
        set(gca, 'XLim', timeLims, 'YLim', 15*[-1 1], 'XTIck', 0:50:500, 'YTick', [-20:5:20], 'FontSize', 14);
        xlabel('\it time (msec)')
        ylabel('\it pAmps');
    end
    
    axis 'square';
    box on; grid on;
    
    subplot(3,4,9)
    dt = timeAxis(2)-timeAxis(1);
    yTicks =[0 :0.05:0.2];
    tt  = (1:size(stimData.osImpulseResponses,1))*dt - dt;
    plot(tt, stimData.osImpulseResponses(:,1), 'r-', 'LineWidth', 1.5); hold on
    plot(tt, stimData.osImpulseResponses(:,2), 'g-', 'LineWidth', 1.5);
    plot(tt, stimData.osImpulseResponses(:,3), 'b-', 'LineWidth', 1.5);
    set(gca, 'XLim', timeLims, 'XTIck', 0:50:500, 'YTick', yTicks, 'YTIckLabel', sprintf('%2.2f\n', yTicks), 'FontSize', 14); 
    axis 'square';
    xlabel('\it time (msec)');
    ylabel('\it pAmps');
    box on; grid on;
    
    % Plot the mosaic with one emPath
    ax = subplot(3,4,8);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'labelConeTypes', false,...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [0 0 0] ...
    );
    hold(ax, 'on');
    plot([-1 1]/1000, [0 0], 'k-', 'LineWidth', 1.5);
    plot([0 0], [-1 1]/1000, 'k-', 'LineWidth', 1.5);
    
    idx = find(timeAxis <= 100);
    plot(theEMpathMeters(idx,1), theEMpathMeters(idx,2), '-', 'LineWidth', 2.0, 'Color', emPathColor);
    hold on;
    if (showDots)
    plot(theEMpathMeters(idx0,1), theEMpathMeters(idx0,2), 'ro', 'LineWidth', 2.0, 'MarkerFaceColor', emPathColor0, 'MarkerSize', 10);
    plot(theEMpathMeters(idx1,1), theEMpathMeters(idx1,2), 'ro', 'LineWidth', 2.0, 'MarkerFaceColor', emPathColor1, 'MarkerSize', 10);
    plot(theEMpathMeters(idx2,1), theEMpathMeters(idx2,2), 'ro', 'LineWidth', 2.0, 'MarkerFaceColor', emPathColor2, 'MarkerSize', 10);
    plot(theEMpathMeters(idx3,1), theEMpathMeters(idx3,2), 'ro', 'LineWidth', 2.0, 'MarkerFaceColor', emPathColor3, 'MarkerSize', 10);
    end
    
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    xlabel('');
    ylabel('')
    set(gca, 'XLim', spaceLimsMeters, 'YLim', spaceLimsMeters);
    
    ax  = subplot(3,4,11);
    edgeColor = [0 0 0];
    faceColor = [0.8 0.8 0.8];
    plottingStyle = 'steps';
    
    hold on
    renderResponseRangeAreaPlot(ax,timeAxis, isomerizationInstance*0-10, isomerizationInstance, zeros(size(isomerizationInstance))-10, ...
        edgeColor, faceColor, 'steps')
    stairs(timeAxis, meanIsomerizations, '-', 'LineWidth', 3, 'Color', [0 0.5 1]);
    stairs(timeAxis, meanIsomerizations, '-', 'LineWidth', 1.5, 'Color', [0 0 1]);
    
    if (showDots)
    plot((timeAxis(idx0)+2.5)*[1 1], [28 isomerizationInstance(idx0)+5], 'r-', 'LineWidth', 1.5);
    plot((timeAxis(idx1)+2.5)*[1 1], [28 isomerizationInstance(idx1)+5], 'r-', 'LineWidth', 1.5);
    plot((timeAxis(idx2)+2.5)*[1 1], [28 isomerizationInstance(idx2)+5], 'r-', 'LineWidth', 1.5);
    plot((timeAxis(idx3)+2.5)*[1 1], [28 isomerizationInstance(idx3)+5], 'r-', 'LineWidth', 1.5);
    
    plot(timeAxis(idx0)+2.5, isomerizationInstance(idx0)+5, 'rv', 'MarkerFaceColor', emPathColor0,'MarkerSize', 10);
    plot(timeAxis(idx1)+2.5, isomerizationInstance(idx1)+5, 'rv', 'MarkerFaceColor', emPathColor1, 'MarkerSize', 10);
    plot(timeAxis(idx2)+2.5, isomerizationInstance(idx2)+5, 'rv', 'MarkerFaceColor', emPathColor2, 'MarkerSize', 10);
    plot(timeAxis(idx3)+2.5, isomerizationInstance(idx3)+5, 'rv', 'MarkerFaceColor', emPathColor3, 'MarkerSize', 10);
    end
    
    set(gca, 'YLim', isomerizationsLimits, 'XLim', timeLims, 'XTick', [0:50:300], 'YTick', [0:5:100], 'FontSize', 14);
    xlabel('\it time (msec)');
    ylabel('\it R*/c/sec');
    box on; grid on;
    axis 'square'
    ax  = subplot(3,4,12);
    plot(timeAxis, photocurrentInstance, 'k-', 'LineWidth', 1.5);  hold  on
    %renderResponseRangeAreaPlot(ax,timeAxis, photocurrentInstance*0+photocurrentLimits(1)-10, photocurrentInstance, ...
    %    zeros(size(photocurrentInstance))+photocurrentLimits(1)-10, edgeColor, faceColor, 'lines')
    
    
    plot(timeAxis, meanPhotocurrent, '-', 'LineWidth', 3, 'Color', [0 0.5 1]);
    plot(timeAxis, meanPhotocurrent, '-', 'LineWidth', 1.5, 'Color', [0 0 1]);
    
    set(gca, 'YLim', photocurrentLimits, 'XLim', timeLims, 'XTick', [0:50:300], 'YTick', [-90:2:-60], 'FontSize', 14);
    axis 'square';
    xlabel('\it time (msec)');
    ylabel('\it pAmps');
    box on; grid on;
    pause
end

function renderResponseRangeAreaPlot(ax,x, yLow, yHigh, yMean, edgeColor, faceColor, plottingStyle)

    if (numel(yLow) ~= numel(x))
        return;
    end
    
    v = [x(1) yLow(1)];
    yMeanTrace = yMean(1);
    xMeanTrace = x(1);
    dt = x(2)-x(1);
    
    timeSamples = numel(x);
    for k = 1:(timeSamples-1)
        if (strcmp(plottingStyle, 'lines'))
            newV = [x(k) yHigh(k)];
            yMeanTraceIncrement = yMean(k);
            xMeanTraceIncrement = x(k);
        else
            newV = [x(k) yHigh(k); x(k+1) yHigh(k)];
            yMeanTraceIncrement = [yMean(k) yMean(k)];
            xMeanTraceIncrement = [x(k) x(k+1)];
        end
        v = cat(1, v, newV);
        yMeanTrace = cat(2, yMeanTrace, yMeanTraceIncrement);
        xMeanTrace = cat(2, xMeanTrace, xMeanTraceIncrement);
    end
    
    if (strcmp(plottingStyle, 'lines'))
        v = cat(1 ,v, [x(timeSamples) yHigh(timeSamples)]);
        yMeanTrace = cat(2, yMeanTrace, yMean(timeSamples));
        xMeanTrace = cat(2, xMeanTrace, x(timeSamples));
    else
         v = cat(1 ,v, [x(timeSamples) yHigh(timeSamples); x(timeSamples)+dt yHigh(timeSamples); x(timeSamples)+dt yLow(timeSamples); x(timeSamples) yLow(timeSamples)]);
         yMeanTrace = cat(2, yMeanTrace, [yMean(k+1) yMean(k+1)]);
         xMeanTrace = cat(2, xMeanTrace, [x(k+1) x(k+1)+dt]);    
    end
    
    for k = (numel(yLow)):-1:2
        if (strcmp(plottingStyle, 'lines'))
            newV = [x(k) yLow(k)];
        else
            newV = [x(k) yLow(k-1); x(k-1) yLow(k-1)];
        end
        v = cat(1, v, newV);
    end
    
    v = cat(1,v, [x(1) yLow(1)]);
    if (strcmp(plottingStyle, 'lines'))
        timeAlignmentFactor = timeSamples/(timeSamples-2);
        v(:,1) = v(:,1)*timeAlignmentFactor;
        xMeanTrace = xMeanTrace*timeAlignmentFactor ;
    end
    
    f = 1:size(v,1);
    desaturation = 0.0;
    alpha = 0.5;
    patch(ax,'Faces',f,'Vertices',v,...
        'FaceAlpha', alpha, ...
        'FaceColor', faceColor*(1-desaturation)+desaturation*[1 1 1], ...
        'EdgeColor', edgeColor*(1-desaturation)+desaturation*[1 1 1], ...
        'LineWidth',1.5)
    hold(ax, 'on');
    plot(ax, xMeanTrace, yMeanTrace, 'k-', 'Color', edgeColor, 'LineWidth', 3);
    
end

function showActivationMap(ax, theMosaic, timeAxis, visualizationTime,  isomerizationInstance, spaceLimsDegs, spaceLimsMeters,signalRange)
    indices = find(timeAxis <= visualizationTime);
    tBin = indices(end);
    theMosaic.renderActivationMap(ax, squeeze(isomerizationInstance(:,tBin)), ...
        'signalRange',  signalRange, 'mapType', 'modulated disks', ...
        'showXLabel', false, 'showYLabel', false, 'showXTicks', false, 'showYTicks', false, ...
        'visualizedFOV', max(spaceLimsDegs)*2);
    hold on
    spaceLimsMeters = spaceLimsMeters + 1.0 * 1e-6*[-1 1];
    
    x = [spaceLimsMeters(1) spaceLimsMeters(2)  spaceLimsMeters(2) spaceLimsMeters(1)  spaceLimsMeters(1) ];
    y = [spaceLimsMeters(1) spaceLimsMeters(1)  spaceLimsMeters(2) spaceLimsMeters(2)  spaceLimsMeters(1) ];
    plot(ax,x,y,'w-', 'LineWidth', 1.5);
    text(ax, spaceLimsMeters(1), spaceLimsMeters(1)*1.1, sprintf('%2.0f msec', visualizationTime), 'FontSize', 14);
end

function idx0 = opticalImageEyePathComboPlot(ax, showDots, timeAxis, visualizationTime, theEMpathDegs, emPathColor, xAxisDegs, yAxisDegs, retinalImage, spaceLimsDegs)
    indices = find(timeAxis <= visualizationTime);
    idx0 = indices(end);
    offset = theEMpathDegs(idx0,:);
    imagesc(ax,xAxisDegs + offset(1), yAxisDegs + offset(2), retinalImage.^0.7); hold(ax, 'on');
    %plot([xAxisDegs(1) xAxisDegs(end)], [0 0], 'k-', 'LineWidth', 1.5);
    %plot([0 0], [yAxisDegs(1) yAxisDegs(end)], 'k-', 'LineWidth', 1.5);
    if (~showDots)
        plot(ax,theEMpathDegs(indices,1), theEMpathDegs(indices,2), '-', 'LineWidth', 1.5, 'Color', 'k');
        plot(ax, offset(1) + [-1 1], offset(2)*[1 1], 'k-', 'LineWidth', 1.0);
        plot(ax, offset(1) * [1 1], offset(2)+[-1 1], 'k-', 'LineWidth', 1.0);
    end
    axis(ax, 'image'); axis(ax, 'xy')
    set(ax, 'CLim', [0 1], 'XLim', spaceLimsDegs, 'YLim', spaceLimsDegs, 'FontSize', 14, ...
        'XTick', [-0.4:0.05:0.4], 'YTick', [-0.4:0.05:0.4], 'XTickLabel', {}, 'YTickLabel', {});
    box(ax, 'on'); grid(ax, 'on');
    text(ax, spaceLimsDegs(1), spaceLimsDegs(1)*1.1, sprintf('%2.0f msec', visualizationTime), 'FontSize', 14);
end

function hFig = visualizeMosaicAndSomeEMpaths(theMosaic, stimData)

    visualizedResponseInstances = 1:4;
    theEMpathMicrons = squeeze(stimData.responseInstanceArray.theMosaicEyeMovementsMicrons(visualizedResponseInstances,:,:));
    theEMpathMeters = theEMpathMicrons / 1e6;
    
    hFig = figure(987); clf;
    set(hFig, 'Position', [10 10 513 600], 'Color', [1 1 1]);
    ax = subplot('Position', [0.01 0.12 0.98 0.895]);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'labelConeTypes', false,...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [0 0 0] ...
    );
    hold(ax, 'on');
    emColors = [0 1 0; 0 1 1; 1 1 0; 1 0.5 0.3];
    for emIndex = 1:numel(visualizedResponseInstances)
        emPathX = squeeze(theEMpathMeters(emIndex,:,1));
        emPathY = squeeze(theEMpathMeters(emIndex,:,2));
        plot(emPathX, emPathY, 'g-', 'LineWidth', 4, 'Color', squeeze(emColors(emIndex,:))*0.7);
        plot(emPathX, emPathY, 'g-', 'LineWidth', 2, 'Color', squeeze(emColors(emIndex,:)));
    end
    
    hold(ax, 'off');
    xtickLabels = {'-0.2', '',  '-0.1', '', '0', '', '+0.1', '', '+0.2'};
    xTicks = (-0.2:0.05:0.2)*theMosaic.micronsPerDegree*1e-6;
    set(gca, 'XTick', xTicks, 'XTickLabel', xtickLabels, 'YTickLabel', {});
    set(gca, 'FontSize', 28);
    xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
    drawnow;
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
        xlabel('\it time (msec)');
        
        set(gca, 'FontSize', 14, 'XTick', 0:50:300, 'XLim', [0 300]);
        box on;
        grid on;
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

    isomerizationsRange(1) = 0;
    
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
             stimData, noStimData, theMosaic, signalName, ...
             visualizedSignalRange, visualizedResponseInstance, makeVideos)
    
    if (strcmp(signalName, 'Isomerization'))
        % Transform rates to counts
        stimData.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations * theMosaic.integrationTime;
        noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations * theMosaic.integrationTime;
        stimData.responseInstanceArray.theMosaicIsomerizations = stimData.responseInstanceArray.theMosaicIsomerizations * theMosaic.integrationTime;
        noStimData.responseInstanceArray.theMosaicIsomerizations = noStimData.responseInstanceArray.theMosaicIsomerizations * theMosaic.integrationTime;

        % Determine visualized time bin of max response
        diffSignal = bsxfun(@minus,stimData.noiseFreeIsomerizations, squeeze(stimData.noiseFreeIsomerizations(:,1)));
        [~,kidx] = max(abs(diffSignal(:)));
        [~,timeBinOfPeakResponse] = ind2sub(size(diffSignal), kidx);
        timeBinOfPeakResponse = timeBinOfPeakResponse - 1;
        fprintf('Isomerization peak response at %2.0f msec\n', stimData.responseInstanceArray.timeAxis(timeBinOfPeakResponse)*1000);
        
        % Determine visualized response range
        signalRange = determineVisualizedIsomerizationsRange(theMosaic, stimData, timeBinOfPeakResponse);
        
        signalLabel = 'cone excitation';
        signalNameFigNo = 0;
    elseif (strcmp(signalName, 'Photocurrent'))

        % Determine visualized time bin of max response
        diffSignal = bsxfun(@minus,stimData.noiseFreePhotocurrents, squeeze(stimData.noiseFreePhotocurrents(:,1)));
        
        [~,kidx] = max(abs(diffSignal(:)));
        [~,timeBinOfPeakResponse] = ind2sub(size(diffSignal), kidx);

        fprintf('Photocurrent peak response at %2.0f msec\n', stimData.responseInstanceArray.timeAxis(timeBinOfPeakResponse)*1000);
        
        % Determine visualized response range
        signalRange = determineVisualizedPhotocurrentsRange(theMosaic, stimData, timeBinOfPeakResponse);
        signalLabel = 'p-current (pA)';
    
        signalNameFigNo = 10000;
    end
    
    
    timeAxis = noStimData.responseInstanceArray.timeAxis;
    
    % The responses to the STIM 
    hFigs = {};
    coneLinePlotType = 'stim';
    hFigsTmp = renderFigure(1222+signalNameFigNo, theMosaic, visualizedResponseInstance, ...
        stimData, noStimData, timeBinOfPeakResponse, ...
        signalRange, signalName, timeAxis, ...
        coneLinePlotType, signalLabel, makeVideos);   
    for hFigIndex = 1:numel(hFigsTmp)
        hFigs{numel(hFigs) + 1} = hFigsTmp{hFigIndex};
    end
    
    % The actual responses to the NULL
    coneLinePlotType = 'null';
    hFigsTmp = renderFigure(1333+signalNameFigNo, theMosaic, visualizedResponseInstance, ...
        noStimData, noStimData, timeBinOfPeakResponse, ...
        signalRange, signalName, timeAxis, ...
        coneLinePlotType, signalLabel, makeVideos);    
    for hFigIndex = 1:numel(hFigsTmp)
        hFigs{numel(hFigs) + 1} = hFigsTmp{hFigIndex};
    end
    
end

function hFig = renderFigure(figNo, theMosaic, visualizedResponseInstance, ...
    stimData, noStimData, timeBinOfPeakResponse, ...
    signalRange, signalName, timeAxis, ...
    coneLinePlotType, yLabelTitle, makeVideos)
    
    % Generate colormaps for modulations and for excitations
    %modulationsColorMap = brewermap(512, '*RdBu');
    excitationsColorMap = gray(1024);

    xtickLabels = {'-0.2', '',  '-0.1', '', '0', '', '+0.1', '', '+0.2'};

    if (strcmp(signalName, 'Isomerization'))
        stimDataNoiseFreeSignalAtPeakTime = squeeze(stimData.noiseFreeIsomerizations(:,timeBinOfPeakResponse));
        noStimDataNoiseFreeSignalAtPeakTime = squeeze(noStimData.noiseFreeIsomerizations(:,timeBinOfPeakResponse));
        stimDataSignalInstanceAtPeakTime = squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(visualizedResponseInstance,:,timeBinOfPeakResponse))';
        stimDataSignalInstanceFull = squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(visualizedResponseInstance,:,:));
        noStimDataNoiseFreeSignalFull = noStimData.noiseFreeIsomerizations;
    elseif (strcmp(signalName, 'Photocurrent'))
        stimDataNoiseFreeSignalAtPeakTime = squeeze(stimData.noiseFreePhotocurrents(:,timeBinOfPeakResponse));
        noStimDataNoiseFreeSignalAtPeakTime = squeeze(noStimData.noiseFreePhotocurrents(:,timeBinOfPeakResponse));
        stimDataSignalInstanceAtPeakTime = squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(visualizedResponseInstance,:,timeBinOfPeakResponse))';
        stimDataSignalInstanceFull = squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(visualizedResponseInstance,:,:));
        noStimDataNoiseFreeSignalFull = noStimData.noiseFreePhotocurrents;
    end

    
    stimDataEMpathMicrons = squeeze(stimData.responseInstanceArray.theMosaicEyeMovementsMicrons(visualizedResponseInstance,:,:));
    noStimDataEMpathMicrons = squeeze(noStimData.responseInstanceArray.theMosaicEyeMovementsMicrons(visualizedResponseInstance,:,:));
    
    if (makeVideos) 
        % Make videos
        if (strcmp(signalName, 'Isomerization'))
            videoFileName = sprintf('%sIsomerizationsMovie', coneLinePlotType);
        else
            videoFileName = sprintf('%sPhotocurrensMovie', coneLinePlotType);
        end

        videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 30;
        videoOBJ.Quality = 100;
        videoOBJ.open();

        hFigVideo = figure(9988); clf;
        set(hFigVideo, 'Position', [10 10 600 900], 'Color', [1 1 1]);
        
        ax = subplot('Position', [0.155 0.72 0.84 0.23]);
        [timelinePlot, timelineArrowPlot] = generateXTConeLinePlot(ax, theMosaic, ...
            stimDataEMpathMicrons, ...
            noStimDataEMpathMicrons, ...
            stimDataSignalInstanceFull, ...
            noStimDataNoiseFreeSignalFull, ...
            signalRange, timeAxis, ...
            coneLinePlotType, signalName, '', 'time (sec)', ...
            false, yLabelTitle);
    
        ax = subplot('Position', [0.155 0.02 0.84 0.70]);
        for tBinIndex = 1:numel(timeAxis)
            if (strcmp(signalName, 'Isomerization'))
                stimDataSignalInstanceAtCurrentTime = ...
                    squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(visualizedResponseInstance,:,tBinIndex))';
            elseif (strcmp(signalName, 'Photocurrent'))
                stimDataSignalInstanceAtCurrentTime = ...
                    squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(visualizedResponseInstance,:,tBinIndex))';
            end
            set(timelinePlot, 'YData', timeAxis(tBinIndex)*[1 1]);
            set(timelineArrowPlot, 'YData', timeAxis(tBinIndex));
            theMosaic.renderActivationMap(ax, stimDataSignalInstanceAtCurrentTime, ...
                'visualizedConeAperture', 'geometricArea', ...
                'mapType', 'modulated disks', ...
                'signalRange', signalRange, ...
                'colorMap', excitationsColorMap, ...
                'showColorBar', ~true, ...
                'labelColorBarTicks', ~true, ...
                'outlineConesAlongHorizontalMeridian', ~true, ...
                'showXLabel', false, ...
                'showYLabel', false, ...
                'backgroundColor', 0*[0.5 0.5 0.5]);
                set(gca, 'XTick', (-0.2:0.05:0.2)*theMosaic.micronsPerDegree*1e-6, ...
                    'XTickLabel', xtickLabels, 'YTickLabel', {});
                set(gca, 'FontSize', 28);
                xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
            drawnow

            % Add video frame
            videoOBJ.writeVideo(getframe(hFigVideo));
        end

        % Close video stream
        videoOBJ.close();
        fprintf('File saved in %s\n', videoFileName);
    end
    
    
    
    
    
    
    hFig = {};
    hFigTmp = figure(1+figNo); clf;
    hFig{numel(hFig)+1} = hFigTmp;
    set(hFigTmp, 'Position', [10 10 513 600], 'Color', [1 1 1]);
    ax = subplot('Position', [0.01 0.12 0.98 0.895]);
    theMosaic.renderActivationMap(ax, stimDataNoiseFreeSignalAtPeakTime ,...
        'visualizedConeAperture', 'geometricArea', ...
        'mapType', 'modulated disks', ...
        'signalRange', signalRange, ...
        'colorMap', excitationsColorMap, ...
        'showColorBar', ~true, ...
        'labelColorBarTicks', ~true, ...
        'outlineConesAlongHorizontalMeridian', ~true, ...
        'showXLabel', false, ...
        'showYLabel', false, ...
        'backgroundColor', 0*[0.5 0.5 0.5]);
    set(gca, 'XTick', (-0.2:0.05:0.2)*theMosaic.micronsPerDegree*1e-6, 'XTickLabel', xtickLabels, 'YTickLabel', {});
    set(gca, 'FontSize', 28);
    xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
   
    drawnow

    hFigTmp = figure(2+figNo); clf;
    hFig{numel(hFig)+1} = hFigTmp;
    set(hFigTmp, 'Position', [10 10 513 600], 'Color', [1 1 1]);
    ax = subplot('Position', [0.01 0.12 0.98 0.895]);
    theMosaic.renderActivationMap(ax, stimDataSignalInstanceAtPeakTime, ...
        'visualizedConeAperture', 'geometricArea', ...
        'mapType', 'modulated disks', ...
        'signalRange', signalRange, ...
        'colorMap', excitationsColorMap, ...
        'showColorBar', ~true, ...
        'labelColorBarTicks', ~true, ...
        'outlineConesAlongHorizontalMeridian', ~true, ...
        'showXLabel', false, ...
        'showYLabel', false, ...
        'backgroundColor', 0*[0.5 0.5 0.5]);
    set(gca, 'XTick', (-0.2:0.05:0.2)*theMosaic.micronsPerDegree*1e-6, ...
        'XTickLabel', xtickLabels, 'YTickLabel', {});
    set(gca, 'FontSize', 28);
    xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
    drawnow

    hFigTmp = figure(3+figNo); clf;
    hFig{numel(hFig)+1} = hFigTmp;
    set(hFigTmp, 'Position', [10 10 600 345], 'Color', [1 1 1]);
    ax = subplot('Position', [0.15 0.25 0.83 0.71]);
    generateConeLinePlot(ax, theMosaic, ...
        stimDataNoiseFreeSignalAtPeakTime, ...
        noStimDataNoiseFreeSignalAtPeakTime, ...
        signalRange, coneLinePlotType, signalName, yLabelTitle, false);
    drawnow

    hFigTmp = figure(4+figNo); clf;
    hFig{numel(hFig)+1} = hFigTmp;
    set(hFigTmp, 'Position', [10 10 600 345], 'Color', [1 1 1]);
    ax = subplot('Position', [0.15 0.25 0.83 0.701]);
    generateConeLinePlot(ax, theMosaic, ...
        stimDataSignalInstanceAtPeakTime, ...
        noStimDataNoiseFreeSignalAtPeakTime, ...
        signalRange, coneLinePlotType, signalName, yLabelTitle, false);
    drawnow

    % Finally the XT ConeLinePlot
    hFigTmp = figure(5+figNo); clf;
    hFig{numel(hFig)+1} = hFigTmp;
    set(hFigTmp, 'Position', [10 10 600 400], 'Color', [1 1 1]);
    ax = subplot('Position', [0.15 0.22 0.83 0.75]);
    generateXTConeLinePlot(ax, theMosaic, ...
        stimDataEMpathMicrons, ...
        noStimDataEMpathMicrons, ...
        stimDataSignalInstanceFull, ...
        noStimDataNoiseFreeSignalFull, ...
        signalRange, timeAxis, 'differential_activations', signalName, ...
        'space (deg)', 'time (sec)', true, yLabelTitle);
    drawnow
        
end

function [timelinePlot, timelineArrowPlot] = generateXTConeLinePlot(ax, theMosaic, EMpath, ...
        nullStimEMpath, activation, nullStimActivation, signalRange, timeAxis, coneLinePlotType, signalName, xLabelTitle, yLabelTitle, showXTicks, colorbarTitle)

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
    coneXcoordsDegs = coneXcoords * 1e6 / theMosaic.micronsPerDegree;
    coneYcoordsDegs = coneYcoords * 1e6 / theMosaic.micronsPerDegree;
    conePosRange = max([max(abs(coneXcoordsDegs)) max(abs(coneYcoordsDegs))])-dx/2;
    % Find cones lying near the y=0 axis
    indicesOfConesAlongXaxis = find(abs(coneYcoordsDegs) < dx);
    coneXcoordsDegs = coneXcoordsDegs(indicesOfConesAlongXaxis);
    coneYcoordsDegs = coneYcoordsDegs(indicesOfConesAlongXaxis);
    identitiesOfConesAlongXaxis = theMosaic.pattern(idx(indicesOfConesAlongXaxis));
    
    
    [~,timeBinsNo] = size(activation);
    coneActivationsXT = zeros(timeBinsNo, numel(indicesOfConesAlongXaxis));
    nullStimActivationXT = zeros(timeBinsNo, numel(indicesOfConesAlongXaxis));
    
    for k = 1:timeBinsNo
        activationTmp = theMosaic.reshapeHex2DmapToHex3Dmap(squeeze(activation(:,k)));
        nullStimActivationTmp = theMosaic.reshapeHex2DmapToHex3Dmap(squeeze(nullStimActivation(:,k)));
        coneActivationsTmp = activationTmp(idx);
        nullStimActivationTmp = nullStimActivationTmp(idx);
        coneActivationsXT(k,:) = coneActivationsTmp(indicesOfConesAlongXaxis);
        nullStimActivationXT(k,:)  = nullStimActivationTmp(indicesOfConesAlongXaxis);
    end
    
    if (strcmp(coneLinePlotType, 'differential_activations'))
        signalXT = coneActivationsXT-nullStimActivationXT;
        if (strcmp(signalName, 'Isomerization'))
            cTicks = -30:5:30;
            cTickLabels = {'-30R*', '-25R*', '-20R*', '-15R*', '-10R*', '-5R*', '0R*', '5R*', '10R*', '15R*','20R*', '25R*', '30R*'};
            signalRange = max(abs(signalXT(:)))*[-0.4 0.4];
        elseif (strcmp(signalName, 'Photocurrent'))
            cTicks = -10:2:10;
            cTickLabels = {'-10pA', '-8pA', '-6pA', '-4pA', '-2pA', '0pA', '2pA', '4pA', '6pA','8pA', '10pA'};
            signalRange = max(abs(signalXT(:)))*[-0.4 0.4];
        end
        
    else
        signalXT = coneActivationsXT;
        if (strcmp(signalName, 'Isomerization'))
            cTicks = 0:10:30;
        elseif (strcmp(signalName, 'Photocurrent'))
            cTicks = -90:2:0;
        end
    end
        
    if (strcmp(coneLinePlotType,'stim')) || (strcmp(coneLinePlotType, 'differential_activations'))
        emTRajectoryXposDegs = squeeze(EMpath(:,1))/theMosaic.micronsPerDegree;
    else
        emTRajectoryXposDegs = squeeze(nullStimEMpath(:,1))/theMosaic.micronsPerDegree;
    end
        
    timeAxisRange = [timeAxis(1) timeAxis(end)];
    timeTicks = -0.2 : 0.05 : 0.2;
    timeTickLabels = {'-.20', '-.15', '-.10', '-.05', '0', '.05', '.10', '.15', '.20'};
    
    imagesc(ax, coneXcoordsDegs, timeAxis, signalXT);
    axis(ax, 'xy');
    hold(ax, 'on');
    plot(emTRajectoryXposDegs, timeAxis, '-', 'LineWidth', 5, 'Color', [0 0.7 0]);
    plot(emTRajectoryXposDegs, timeAxis, '-', 'LineWidth', 3, 'Color', [0 1 0]);
    
    timelinePlot = plot(ax, [min(coneXcoordsDegs) max(coneXcoordsDegs)], [nan nan], 'g-', 'LineWidth', 2);
    timelineArrowPlot = plot(ax, min(coneXcoordsDegs), nan, 'gs', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', [0.6 1.0 0.6]);
    hold(ax, 'off');
    if (showXTicks)
        xTickLabels = {'-0.2', '', '-0.1', '', '0', '', '+0.1', '', '+0.2'};
    else
        xTickLabels = {};
    end
    
    set(ax, 'XLim', conePosRange*[-1 1], 'YLim', [timeAxisRange(1)-2.5/1000 timeAxisRange(2)+2.5/1000], 'CLim', signalRange, ...
        'FontSize', 28, 'XTick', -0.2:0.05:0.2, 'XTickLabel', xTickLabels, ...
        'YTick', timeTicks, 'YTickLabel', timeTickLabels, 'LineWidth', 1.0);
    if (~isempty(yLabelTitle))
        ylabel(sprintf('\\it %s',yLabelTitle), 'FontWeight', 'normal', 'FontSize', 36)
    else
        set(gca, 'YTickLabel', {});
    end
    colormap(gray(1024));
    grid 'on'; box 'on';
    if (~isempty(xLabelTitle))
        xlabel(sprintf('\\it %s',xLabelTitle), 'FontWeight', 'normal', 'FontSize', 36);
    end
    
    % Colorbar
    hcb = colorbar('northoutside', 'Ticks',cTicks,'TickLabels', cTickLabels);
    
    %colorTitleHandle = get(hcb,'Title');
    %set(colorTitleHandle ,'String', colorbarTitle);
    set(hcb, 'FontSize', 22);
end

function generateConeLinePlot(ax, theMosaic, activation, nullStimActivation, signalRange, coneLinePlotType, signalName, yLabelTitle, exportNullProfile)
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
        
    idx = find(theMosaic.pattern == 2);
    a = nullStimActivation(idx);
    meanLconeActivation = mean(a(:));
    
    idx = find(theMosaic.pattern == 3);
    a = nullStimActivation(idx);
    meanMconeActivation = mean(a(:));
    
    idx = find(theMosaic.pattern == 4);
    a = nullStimActivation(idx);
    meanSconeActivation = mean(a(:));
    
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
                fname = 'LconeProfile';
                meanActivation = meanLconeActivation;
            case 3
                edgeColor = [0 0.8 0];
                faceColor = [0.5 1 0.5];
                fname = 'MconeProfile';
                meanActivation = meanMconeActivation;
            case 4
                edgeColor = [0 0 1];
                faceColor = [0.5 0.8 1];
                fname = 'SconeProfile';
                meanActivation = meanSconeActivation;
        end
        
        iidx = find(identitiesOfConesAlongXaxis == coneIndex);
        if (strcmp(coneLinePlotType, 'differential_activations'))
            signal = coneActivations(iidx)-nullStimActivation(iidx);
            if (strcmp(signalName, 'Isomerization'))
                yticks = -10:5:10;
            elseif (strcmp(signalName, 'Photocurrent'))
                yticks = -20:2:20;
            end
        else
            signal = coneActivations(iidx);
            if (strcmp(signalName, 'Isomerization'))
                yticks = 0:10:30;
            elseif (strcmp(signalName, 'Photocurrent'))
                yticks = -90:2:0;
            end
        end
        
        
        if (strcmp(signalName, 'Isomerization')) && (exportNullProfile)
            if (coneIndex < 4)
            coneXpos = coneXcoordsDegs(iidx);
            if (strcmp(fname,'LconeProfile'))
                for k = 1:numel(signal)
                    fprintf('%f %f\n', coneActivations(iidx(k)), nullStimActivation(iidx(k)));
                end
            end
            
            fullFname = sprintf('%s_%s.mat',fname, coneLinePlotType);
            save(fullFname,'coneXpos', 'signal', 'meanActivation');
            fprintf('Saved profile to %s\n', fullFname);
            end
        end
    
        plot(ax, coneXcoordsDegs(iidx), signal, ...
            'o', 'Color', edgeColor, 'MarkerSize', 14, 'LineWidth', 1.5, ...
            'MarkerFaceColor',faceColor,'MarkerEdgeColor',edgeColor);
    end
    xtickLabels = {'-0.2', '', '-0.1', '', '0', '', '+0.1', '', '+0.2'};
    set(ax, 'XLim', conePosRange*[-1 1], 'YLim', signalRange, ...
        'FontSize', 28, 'XTick', -0.2:0.05:0.2, 'XTickLabel', xtickLabels, ...
        'YTick', yticks, 'LineWidth', 1.0);
    if (~isempty(yLabelTitle))
        ylabel(sprintf('\\it %s',yLabelTitle), 'FontWeight', 'normal', 'FontSize', 36)
    else
        set(gca, 'YTickLabel', {});
    end
    grid 'on'; box 'off';
    xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
    drawnow

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
    isomerizationsRange = prctile(stimData.responseInstanceArray.theMosaicIsomerizations(:), [0 95]);
end


function photocurrentsRange = determineVisualizedPhotocurrentsRange(theMosaic, stimData, timeBinOfPeakResponse)
    % Determine response range based on the responses of L/M cones
    nonNullCones = theMosaic.pattern(theMosaic.pattern > 1);
    lmConeIndices = find(nonNullCones==2 | nonNullCones==3);
    
    noiseFreeRange = [...
        min(stimData.noiseFreePhotocurrents(lmConeIndices,timeBinOfPeakResponse)) ...
        max(stimData.noiseFreePhotocurrents(lmConeIndices,timeBinOfPeakResponse)) ...
        ];
%     instancesRange = [ ...
%         min(min(squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(:,lmConeIndices,timeBinOfPeakResponse)))) ...
%         max(max(squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(:,lmConeIndices,timeBinOfPeakResponse)))) ...
%     ];

    photocurrentsRange = noiseFreeRange +[-3 3];
    
end
