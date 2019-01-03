function hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
    timeAxis, noStimResponseInstances, stimResponseInstances, ...
    noStimNoiseFreeResponse, stimNoiseFreeResponse, responseRange, ...
    signalSource, yAxisLabel, figNo)

    stimOnsetTime = timeAxis(end);
    dt = timeAxis(2)-timeAxis(1);
    compositeResponseTimeAxis = cat(2, timeAxis, timeAxis-timeAxis(1)+dt+timeAxis(end)); 
    stimOnsetTime = stimOnsetTime - compositeResponseTimeAxis(1)+dt;
    compositeResponseTimeAxis = compositeResponseTimeAxis - compositeResponseTimeAxis(1);
    
    compositeNoiseFreeResponses = cat(2, noStimNoiseFreeResponse, stimNoiseFreeResponse);
    compositeResponseInstances = cat(3, noStimResponseInstances, stimResponseInstances);
    pause
    
    timeTicks = [-500:50:500];
    if (strcmp(signalSource, 'isomerizations'))
        areaPlotMode = 'steps';
        yticks = 0:5:50;
    elseif (strcmp(signalSource, 'photocurrents'))
         areaPlotMode = 'lines';
        yticks = -100:5:0;
    end
    
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 600 345], 'Color', [1 1 1]);
    ax = subplot('Position', [0.15 0.25 0.83 0.71]);
    
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
        renderAreaPlot(ax,compositeResponseTimeAxis,min(instances,[],1),max(instances,[],1), meanResponse, edgeColor, faceColor, areaPlotMode);
        hold on;
        plot(ax, stimOnsetTime*[1 1], responseRange, 'k-', 'LineWidth', 1.5);
    end       
    
    
            
    set(ax, 'XLim', [compositeResponseTimeAxis(1)-dt compositeResponseTimeAxis(end)+dt], 'YLim', responseRange, ...
        'FontSize', 28, 'XTick', timeTicks, ...
        'YTick', yticks, 'LineWidth', 1.0);
    if (~isempty(yAxisLabel))
        ylabel(sprintf('\\it %s',yAxisLabel), 'FontWeight', 'normal', 'FontSize', 36)
    end
    xlabel('\it time (msec)', 'FontWeight', 'normal', 'FontSize', 36);
    grid on;
    
    drawnow
    
end

function renderAreaPlot(ax,x, yLow, yHigh, yMean, edgeColor, faceColor, mode)

    v = [x(1) yLow(1)];
    yMeanTrace = yMean(1);
    xMeanTrace = x(1);
    dt = x(2)-x(1);
    
    for k = 1:(numel(yHigh)-1)
        if (strcmp(mode, 'lines'))
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
    
    if (strcmp(mode, 'lines'))
        v = cat(1 ,v, [x(k+1) yHigh(k+1)]);
        yMeanTrace = cat(2, yMeanTrace, yMean(k+1));
        xMeanTrace = cat(2, xMeanTrace, x(k+1));
    else
        v = cat(1 ,v, [x(k+1) yHigh(k+1); x(k+1)+dt yHigh(k+1); x(k+1)+dt yLow(k+1); x(k+1) yLow(k+1)]);
        yMeanTrace = cat(2, yMeanTrace, [yMean(k+1) yMean(k+1)]);
        xMeanTrace = cat(2, xMeanTrace, [x(k+1) x(k+1)+dt]);
    end
    
    for k = (numel(yLow)):-1:2
        if (strcmp(mode, 'lines'))
            newV = [x(k) yLow(k)];
        else
            newV = [x(k) yLow(k); x(k-1) yLow(k)];
        end
        v = cat(1, v, newV);
    end
    
    v = cat(1,v, [x(1) yLow(1)]);
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

