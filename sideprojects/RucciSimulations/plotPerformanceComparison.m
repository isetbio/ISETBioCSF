function plotPerformanceComparison(performance, performanceSTABILIZED, stimDescriptor, signalName, figNo)

    % Plot the raw and fitted psychometric function
    hFig = figure(figNo+1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 400]);
    
    
    % The smooth (fitted) psychometric function
    plot(performance.smoothPsychometricFunction.contrast*100, performance.smoothPsychometricFunction.performance, 'r-', 'LineWidth', 1.5); hold on;
    plot(performanceSTABILIZED.smoothPsychometricFunction.contrast*100, performanceSTABILIZED.smoothPsychometricFunction.performance, 'k-', 'LineWidth', 1.5); hold on;
    % The raw (measured) psychometric function
    
    % The contrast threshold
    %plot(performance.contrastThreshold*[1 1]*100, [0 performance.performanceCriterion]*100, 'b-', 'LineWidth', 1.5);
    %plot([0.0001 performance.contrastThreshold]*100, performance.performanceCriterion*[1 1]*100, 'b-', 'LineWidth', 1.5);
    plot(performance.contrastThreshold*100, performance.performanceCriterion*100, 'rx');
    plot(performanceSTABILIZED.contrastThreshold*100, performanceSTABILIZED.performanceCriterion*100, 'kx');
     
    plot(performance.rawPsychometricFunction.contrast*100, performance.rawPsychometricFunction.performance, 'ko', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.8 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'LineWidth', 1.0);
    
    plot(performanceSTABILIZED.rawPsychometricFunction.contrast*100, performanceSTABILIZED.rawPsychometricFunction.performance, 'ko', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.0);
    
    contrastLims = [min(performance.rawPsychometricFunction.contrast) max(performance.rawPsychometricFunction.contrast)]*100;
    axis 'square'
    set(gca, 'XScale', 'log', 'XLim', contrastLims, 'YLim', [0.4 1.01]*100, 'XScale', 'log', 'FontSize', 16);
    set(gca, 'XTick', [0.03 0.1 0.3 1 3]);
    grid on;
    legend({...
        sprintf('dynamic,    cT=%2.2f%%', 100*performance.contrastThreshold), ...
        sprintf('stabilized, cT=%2.2f%%', 100*performanceSTABILIZED.contrastThreshold)'}, 'Location', 'SouthEast');
    
    xlabel('\it contrast (%)');
    ylabel('\it classification accuracy (%)');
    title(sprintf('%s\n signal: %s', stimDescriptor, signalName), 'FontWeight', 'normal');
    
end
