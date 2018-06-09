function timeAxisLimits = renderNullTestComboResponse(ax1, ax2, signalSource, noStimResponses, stimResponses, noStimResponsesNoiseFree, stimResponsesNoiseFree, responseLevels, timeAxis, plotType, row, col, rows)
    
    comboResponses = cat(2, noStimResponses, stimResponses);
    if ((isempty(noStimResponsesNoiseFree)) || (isempty(stimResponsesNoiseFree)))
        meanResponse = mean(comboResponses,1);
    else
        meanResponse = cat(2, noStimResponsesNoiseFree, stimResponsesNoiseFree);
    end
    responseLevels = double(responseLevels);
    dt = timeAxis(2)-timeAxis(1);
    timeAxisComboSecondResponseStart = timeAxis(end)+dt;
    timeAxisCombo = double(cat(2, timeAxis, timeAxis(end)+timeAxis))+dt;
    tFirstHalfIndices  = find(timeAxisCombo <= timeAxis(end));
    tSecondHalfIndices = find(timeAxisCombo > timeAxis(end));
    
     
    timeAxisLimits = [timeAxisCombo(1) timeAxisCombo(end)];
    if (timeAxisCombo(1) == timeAxisCombo(end))
       timeAxisLimits = timeAxisCombo(1) + [-0.5 0.5];
    end
       
    if (strcmp(plotType, 'density'))
        responseDistribution = compute2Dhistogram(comboResponses, responseLevels);
        pcolor(ax1, timeAxisCombo, responseLevels(1:end-1), responseDistribution);
        hold(ax1, 'on');
    elseif (strcmp(plotType, 'line'))
        hold(ax1, 'on');
        for instanceIndex = 1:size(responseInstancesNum,1)
            y = squeeze(comboResponses(instanceIndex, :));
            xflip = [timeAxisCombo(1 : end - 1) fliplr(timeAxisCombo)];
            yflip = [y(1 : end - 1) fliplr(y)];
            patch(ax1, xflip, yflip, [0.1 0.1 0.1], 'EdgeAlpha', 0.1, 'FaceColor', 'none');
        end
    end
    plot(ax1, timeAxisComboSecondResponseStart*[1 1], [responseLevels(1) responseLevels(end)], 'k-', 'LineWidth', 2.0);
    dy = (responseLevels(2)-responseLevels(1))/2;
    plot(ax1, timeAxisCombo, meanResponse-dy, 'k-', 'LineWidth', 4.0);
    plot(ax1, timeAxisCombo, meanResponse-dy, 'g-', 'LineWidth', 3.0);
    
    hold(ax1, 'off');
    grid(ax1, 'on'); box(ax1, 'on');
   
    if (row == rows) && (col == 1)
        xlabel(ax1,'time (msec)', 'FontWeight', 'bold');
        if (strcmp(signalSource, 'isomerizations'))
            ylabel(ax1,'isomerizations (R*/cone)', 'FontWeight', 'bold');
        else
            ylabel(ax1,'photocurrent (pAmps)', 'FontWeight', 'bold');
        end
    else
        set(ax1, 'XTickLabel', {}, 'YTickLabel', {});
    end
    set(ax1, 'XLim', [timeAxisLimits(1) timeAxisLimits(2)]);
    set(ax1, 'YLim', [responseLevels(1) responseLevels(end-1)]); 
     
    maxProbability = max(responseDistribution(:));
    set(ax1, 'CLim', [0 maxProbability]);
    
    xTicks = 0:25:1000;
    set(ax1, 'XTick', xTicks, 'XTickLabel', xTicks);
    
    if (~isempty(ax2))
        textOffset = (timeAxisCombo(end)-timeAxisCombo(1))/6;
        x = timeAxisCombo(tFirstHalfIndices(1)) + textOffset;
        y = 1.01*responseLevels(end);
        t = text(ax1, x,y, 'null stimulus');
        t(1).FontSize = 18;
        x = timeAxisCombo(tSecondHalfIndices(1)) + textOffset;
        t = text(ax1, x, y, 'test stimulus');
        t(1).FontSize = 18;
        t(1).Color = 'r';
        
        meanResponseDistributionNullStimulus = mean(responseDistribution(:,tFirstHalfIndices),2);
        meanResponseDistributionTestStimulus = mean(responseDistribution(:,tSecondHalfIndices),2);
     
        plot(ax2, meanResponseDistributionNullStimulus, responseLevels(1:end-1), ...
            'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.6 0.6 0.6], 'LineWidth', 1.5);
        hold(ax2, 'on');
        plot(ax2, meanResponseDistributionTestStimulus, responseLevels(1:end-1), ...
            'ro-', 'MarkerSize', 10, 'MarkerFaceColor', [1.0 0.6 0.6], 'LineWidth', 1.5);

        if (row == rows) && (col == 1)
            %xlabel(ax2,'count', 'FontWeight', 'bold');
        else
            set(ax2, 'XTickLabel', {}, 'YTickLabel', {});
        end
        maxAll = max([ ...
            max(meanResponseDistributionNullStimulus) ...
            max(meanResponseDistributionTestStimulus) ...
            ]);
        set(ax2, 'XLim', [0 maxAll]);
        set(ax2, 'YLim', [responseLevels(1) responseLevels(end-1)]); 
        set(ax2, 'XTick', [0 maxAll], 'XTickLabel', {'0', sprintf('%2.2f', maxAll)}, 'YTickLabel', {});
        xlabel('probability', 'FontSize', 18);
    end
    
    colormap(1-gray(1024));
    drawnow; 
end

function responseDistribution = compute2Dhistogram(responses, responseLevels)
    nTrials = size(responses,1);
    responseDistribution = zeros(numel(responseLevels)-1, size(responses,2));
    for tBin = 1:size(responses,2)
        N = histcounts(squeeze(responses(:,tBin)),responseLevels);
        responseDistribution(:,tBin) = N';
    end
    % Force the sum of responseDistribution at each time bin to equal to 100
    responseDistribution = responseDistribution/nTrials;
end