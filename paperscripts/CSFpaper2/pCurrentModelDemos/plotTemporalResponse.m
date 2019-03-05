function plotTemporalResponse(timeAxis, signal, lineColor, signalName, plotType, labelXaxis, labelYaxis, yLims, yTicks, yTickLabels)
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
    
    xTicks = [0:100:8000];
    if (labelXaxis)
        xlabel('\it time (msec)');
        if (timeAxis(end) > 1000)
            xTickLabels = {'0', '', '200', '', '400', '', '600', '', '800', '', ...
                    '1000', '', '1200', '', '1400', '', '1600', '', '1800', '', ...
                    '2000', '', '2200', '', '2400', '', '2600', '', '2800', '', ...
                    '3000', '', '3200', '', '3400', '', '3600', '', '3800', '', ...
                    '4000', '', '4200', '', '4400', '', '4600', '', '4800', '', ...
                    '5000', '', '5200', '', '5400', '', '5600', '', '5800', '', ...
                    '6000', '', '6200', '', '6400', '', '6600', '', '6800', '', ...
                    '7000', '', '7200', '', '7400', '', '7600', '', '7800', ''};
        else
            xTickLabels = sprintf('%2.0f\n', xTicks);
        end
    else
        xTickLabels = {};
    end
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)]*1000, 'XTick', xTicks, 'XTickLabel', xTickLabels, 'FontSize', 14);
    if (~isempty(yLims))
        set(gca, 'YLim', yLims, 'YTick', yTicks, 'YTickLabel', yTickLabels);
    else
        set(gca, 'YTick', [], 'XTick', []);
        %set(gca, 'XColor', 'none', 'YColor', 'none');
    end
    grid off; box on;
end
