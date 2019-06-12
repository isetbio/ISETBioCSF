function plotPerformance(performance, stimDescriptor, eyePosition, signalName, figNo)

    % Plot the raw and fitted psychometric function
    hFig = figure(figNo+1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 400]);
    
    
    % The smooth (fitted) psychometric function
    plot(performance.smoothPsychometricFunction.contrast*100, performance.smoothPsychometricFunction.performance, 'r-', 'LineWidth', 1.5); hold on;
    % The raw (measured) psychometric function
    
    % The contrast threshold
    plot(performance.contrastThreshold*[1 1]*100, [0 performance.performanceCriterion]*100, 'b-', 'LineWidth', 1.5);
    plot([0.0001 performance.contrastThreshold]*100, performance.performanceCriterion*[1 1]*100, 'b-', 'LineWidth', 1.5);
    plot(performance.contrastThreshold*100, performance.performanceCriterion*100, 'bx');
     
    plot(performance.rawPsychometricFunction.contrast*100, performance.rawPsychometricFunction.performance, 'ko', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.8 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'LineWidth', 1.0);
    
    contrastLims = [min(performance.rawPsychometricFunction.contrast) max(performance.rawPsychometricFunction.contrast)]*100;
    axis 'square'
    set(gca, 'XScale', 'log', 'XLim', contrastLims, 'YLim', [0.4 1.01]*100, 'XScale', 'log', 'FontSize', 16);
    set(gca, 'XTick', [0.03 0.1 0.3 1 3]);
    grid on;
    
    t = text(performance.contrastThreshold*100, performance.performanceCriterion*100, ...
        sprintf('%s C_t = %2.2f%%', '\leftarrow', 100*performance.contrastThreshold), ...
        'Color', 'b', 'FontSize',14);
    %t.fontSize = 30;
    xlabel('\it contrast (%)');
    ylabel('\it classification accuracy (%)');
    title(sprintf('%s, %s\n signal: %s', stimDescriptor, eyePosition, signalName), 'FontWeight', 'normal');
    
end
