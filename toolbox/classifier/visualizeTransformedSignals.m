function hFig = visualizeTransformedSignals(timeAxis, noStimResponses, stimResponses, signalSource, stimContrast, spatialFilterName)

    hFig = figure(1234); clf;
    set(hFig, 'Position', [10 10 850 800], 'Color', [1 1 1]);

    nTrials = size(noStimResponses, 1);
    responseRange = [...
            min([min(stimResponses(:)) min(noStimResponses(:))]) ...
            max([max(stimResponses(:)) max(noStimResponses(:))])
        ];
    
    q = 25/100;
    quantileResponses = quantile(noStimResponses, [q 1-q], 1);
    x = [timeAxis(1) timeAxis fliplr(timeAxis)];
    f = 1:size(x,2);
    y = [quantileResponses(1,1) quantileResponses(1,:) fliplr(quantileResponses(2,:))];
    maxY = max(y);
    minY = min(y);
    noStimResponsesRange = [x(:) y(:)];
    quantileResponses = quantile(stimResponses, [q 1-q], 1);
    y = [quantileResponses(1,1) quantileResponses(1,:) fliplr(quantileResponses(2,:))];
    maxY = max([maxY max(y)]);
    minY = min([minY min(y)]);
    stimResponsesRange = [x(:) y(:)];
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 2, ...
           'heightMargin',   0.09, ...
           'widthMargin',    0.06, ...
           'leftMargin',     0.08, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.03);

    subplot('Position', subplotPosVectors(1,1).v);
    plot(timeAxis, noStimResponses, 'k-');
    set(gca, 'YLim', responseRange,  'XLim', [timeAxis(1) timeAxis(end)]);
    box off; grid on
    title(sprintf('NO-STIM, nTrials = %d', nTrials));
    ylabel(sprintf('%s output\n(%s-based)', spatialFilterName, signalSource), 'FontWeight', 'bold');
    set(gca, 'FontSize', 14);

    subplot('Position', subplotPosVectors(1,2).v);
    plot(timeAxis, stimResponses, 'k-');
    set(gca, 'YLim', responseRange,  'XLim', [timeAxis(1) timeAxis(end)]);
    box off; grid on
    title(sprintf('C = %2.3f%%, nTrials = %d', stimContrast, nTrials));
    set(gca, 'FontSize', 14);
     
    subplot('Position', subplotPosVectors(2,1).v);
    patch('Faces', f, 'Vertices', noStimResponsesRange, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
    set(gca, 'YLim', [minY-0.1*abs(minY) maxY+0.1*abs(maxY)],  'XLim', [timeAxis(1) timeAxis(end)]);
    box off; grid on
    ylabel(sprintf('%s output\n(%s-based)', spatialFilterName, signalSource), 'FontWeight', 'bold');
    title(sprintf('quantile=%2.2f', q));
    set(gca, 'FontSize', 14);
    xlabel('time (ms)', 'FontWeight', 'bold');
    
    subplot('Position', subplotPosVectors(2,2).v);
    patch('Faces', f, 'Vertices', stimResponsesRange, 'FaceColor', [1 0.2 0.3], 'FaceAlpha', 0.5);
    hold off;
    set(gca, 'YLim', [minY-0.1*abs(minY) maxY+0.1*abs(maxY)],  'XLim', [timeAxis(1) timeAxis(end)]);
    box off; grid on
    title(sprintf('quantile=%2.2f', q));
    set(gca, 'FontSize', 14);
    xlabel('time (ms)', 'FontWeight', 'bold');
    drawnow;
end

