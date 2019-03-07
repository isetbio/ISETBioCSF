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
    
    if (timeAxis(end)>= 0.99)
        xTicks = [0:250:8000];
        xTickLabels = {'0', '', '500', '',  '1000', '', ...
                    '1500', '', '2000', '', '2500', '', '3000', '', '3500', '', ...
                    '4000', '', '4500', '', '5000', '', '5500', '', '6000', '', ...
                    '6500', '', '7000', '', '7500', '', '8000'};
    elseif (timeAxis(end)>0.4)
        xTicks = [0:100:8000];
        xTickLabels = {'0', '', '200', '', '400', '', '600', '', '800', '', ...
                    '1000', '', '1200', '', '1400', '', '1600', '', '1800', '', ...
                    '2000', '', '2200', '', '2400', '', '2600', '', '2800', '', ...
                    '3000', '', '3200', '', '3400', '', '3600', '', '3800', '', ...
                    '4000', '', '4200', '', '4400', '', '4600', '', '4800', '', ...
                    '5000', '', '5200', '', '5400', '', '5600', '', '5800', '', ...
                    '6000', '', '6200', '', '6400', '', '6600', '', '6800', '', ...
                    '7000', '', '7200', '', '7400', '', '7600', '', '7800', ''};
    else
        xTicks = [0:50:4000];
        xTickLabels = {'0', '', '100', '', '200', '', '300', '', '400', '', ...
                    '500', '', '600', '', '700', '', '800', '', '900', '', ...
                    '1000', '', '1100', '', '1200', '', '1300', '', '1400', '', ...
                    '1500', '', '1600', '', '1700', '', '1800', '', '1900', '', ...
                    '2000', '', '2100', '', '2200', '', '2300', '', '2400', '', ...
                    '2500', '', '2600', '', '2700', '', '2800', '', '2900', '', ...
                    '3000', '', '3100', '', '3200', '', '3300', '', '3400', '', ...
                    '3500', '', '3600', '', '3700', '', '3800', '', '3900', ''};
    end
        
    if (labelXaxis)
        xlabel('\it time (msec)');
    end
    
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)]*1000, 'XTick', xTicks, 'XTickLabel', xTickLabels, 'FontSize', 14);
    if (~isempty(yLims))
        set(gca, 'YLim', yLims, 'YTick', yTicks, 'YTickLabel', yTickLabels);
    else
        if (max(signal) ~= min(signal))
            yTicks = linspace(min(signal), max(signal), 4);
        else
            yTicks = [];
        end
        set(gca, 'YTick', yTicks, 'YTickLabel', {});
        %set(gca, 'XColor', 'none', 'YColor', 'none');
    end
    grid on; box on;
end
