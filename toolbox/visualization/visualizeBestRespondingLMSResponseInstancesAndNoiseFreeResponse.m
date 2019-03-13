function hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
    timeAxis, noStimResponseInstances, stimResponseInstances, ...
    noStimNoiseFreeResponse, stimNoiseFreeResponse, responseRange, instancesPercentRange, ...
    signalSource, yAxisLabel, figNo)

    dt = timeAxis(2)-timeAxis(1);
    stimOnsetTime = timeAxis(end)+dt;
    compositeResponseTimeAxis = cat(2, timeAxis, timeAxis-timeAxis(1)+dt+timeAxis(end)); 
    
    % Place stimOnset at 0 msec
    compositeResponseTimeAxis = compositeResponseTimeAxis - stimOnsetTime;
    stimOnsetTime = 0;
      
    compositeResponseTimeAxis = compositeResponseTimeAxis / 1000;
    stimOnsetTime = stimOnsetTime/1000;
    
    compositeNoiseFreeResponses = cat(2, noStimNoiseFreeResponse, stimNoiseFreeResponse);
    compositeResponseInstances = cat(3, noStimResponseInstances, stimResponseInstances);
    
    timeTicks = [-250 -200 -150 -100 -50 0 50 100 150 200 250]/1000;
    if (max(compositeResponseTimeAxis) < 0.2)
        timeTickLabels = {' ', '-0.2', ' ', '-0.1', '-.05', '0', '.05', '0.1', ' ', '0.2', ' '};
    else
        timeTickLabels = {' ', '-0.2', ' ', '-0.1', ' ', '0', ' ', '0.1', ' ', '0.2', ' '};
    end
    timeLimits = max(abs(compositeResponseTimeAxis))*[-1 1];
    
    if (strcmp(signalSource, 'isomerizations'))
        responsePlottingStyle = 'steps';
        yticks = 0:5:50;
    elseif (strcmp(signalSource, 'photocurrents'))
        responsePlottingStyle = 'lines';
        yticks = -100:5:0;
    end
    
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 600 345], 'Color', [1 1 1]);
    ax = subplot('Position', [0.15 0.25 0.83 0.70]);
    
    for coneTypeIndex = 1:3
        switch (coneTypeIndex)
            case 1
                edgeColor = [1 0 0];
                faceColor = [1 0.5 0.5];
            case 2
                edgeColor = [0 0.8 0];
                faceColor = [0.5 1 0.5];
            case 3
                edgeColor = [0 0 1];
                faceColor = [0.5 0.8 1];
        end
        instances = squeeze(compositeResponseInstances(coneTypeIndex,:,:));
        meanResponse = squeeze(compositeNoiseFreeResponses(coneTypeIndex,:));
     
        rRange = prctile(instances, instancesPercentRange, 1);
        rLow = squeeze(rRange(1,:));
        rHigh = squeeze(rRange(2,:));

        
        renderResponseRangeAreaPlot(ax,compositeResponseTimeAxis,rLow, rHigh, meanResponse, edgeColor, faceColor, responsePlottingStyle);
        hold on;
        plot(ax, stimOnsetTime*[1 1], responseRange, 'k-', 'LineWidth', 1.5);
    end       
    
    
    set(ax, 'XLim', timeLimits, 'XTick', timeTicks, 'XTickLabels', timeTickLabels, ...
        'YLim', responseRange, 'YTick', yticks, ...
        'FontSize', 28, 'LineWidth', 1.0);
    
    if (~isempty(yAxisLabel))
        ylabel(sprintf('\\it %s',yAxisLabel), 'FontWeight', 'normal', 'FontSize', 36)
    end
    xlabel('\it time (sec)', 'FontWeight', 'normal', 'FontSize', 36);
    grid on;
    
    drawnow
    
end


function renderResponseRangeAreaPlot(ax,x, yLow, yHigh, yMean, edgeColor, faceColor, plottingStyle)

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
        'LineWidth',1.0)
    hold(ax, 'on');
    plot(ax, xMeanTrace, yMeanTrace, 'k-', 'Color', edgeColor, 'LineWidth', 3);
    
end

