function plotTemporalResponse(timeAxis, signal, lineColor, signalName, plotType, labelXaxis, labelYaxis)
    if (strcmp(plotType, 'stem'))
        stem(timeAxis*1000, signal, 'Color', lineColor, 'MarkerFaceColor', lineColor, 'MarkerSize', 6, 'LineWidth', 1.0);
    elseif (strcmp(plotType, 'line'))
        plot(timeAxis*1000, signal, '-', 'Color', lineColor, 'LineWidth', 1.5);
    elseif (strcmp(plotType, 'dashed line'))
        plot(timeAxis*1000, signal, '--', 'Color', lineColor, 'LineWidth', 1.5);
    elseif (strcmp(plotType, 'dotted line'))
        plot(timeAxis*1000, signal, ':', 'Color', lineColor, 'LineWidth', 1.5);
    else
        error('Uknown plotType ''%s'' in plotTemporalResponse', plotType);
    end
    if (labelYaxis)
        ylabel(sprintf('\\it %s',signalName));
    end
    
    xTicks = [0 100 200 300 400 500 600 700 800 900 1000];
    if (labelXaxis)
        xlabel('\it time (msec)');
        xTickLabels = {'0', '', '200', '', '400', '', '600', '', '800', '', '1000'};
    else
        xTickLabels = {};
    end
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)]*1000, 'XTick', xTicks, 'XTickLabel', xTickLabels, 'FontSize', 14);
end
