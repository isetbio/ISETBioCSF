function makeFixationalEyeMovementFigure()
    
    % Generate 512 3-second emPaths with a sample time of 1 msec
    emDurationSeconds = 3;
    sampleTimeSeconds = 1 / 1000;
    nTrials = 1024;
    
    microSaccadeType = 'none'; % , 'heatmap/fixation based', 'none'
    fixEMobj = fixationalEM();
    fixEMobj.randomSeed = 1;
    fixEMobj.microSaccadeType = microSaccadeType;
       
    % Compute the emPaths
    computeVelocity = true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
            nTrials, computeVelocity, 'useParfor', true);
        
    maxDurationSecondsForFixationSpan = 0.15;
    emPosRange = 5;
    plotData(fixEMobj, maxDurationSecondsForFixationSpan, emPosRange);
end


function plotData(fixEMobj, maxDurationSecondsForFixationSpan, emPosRange )

    % Compute fixationMap
    emPosDelta = 0.6;
    [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice] = ...
        fixEMobj.computeFixationMap(fixEMobj.timeAxis, ...
        fixEMobj.emPosArcMin, emPosRange*[-1 1], emPosDelta, ...
        'maxDurationSeconds', maxDurationSecondsForFixationSpan);
    
    % Compute dispacements
    displacementMode = 'D2';
    nTrials = size(fixEMobj.emPosArcMin,1);
    for trialNo = 1:nTrials
        xPosArcMin = squeeze(fixEMobj.emPosArcMin(trialNo,:,1));
        yPosArcMin = squeeze(fixEMobj.emPosArcMin(trialNo,:,2));
        [meanD(trialNo,:), ~, meanDscrambled(trialNo,:), ...
             DrandomWalk, timeLagsMilliseconds] = ...
        fixEMobj.performDisplacementAnalysis(yPosArcMin, fixEMobj.timeAxis, 'mode', displacementMode);
    end
    meanD = mean(meanD,1);
    meanDscrambled = mean(meanDscrambled,1);
    
    % Make figure
    hFig = figure(1);
    clf;
    set(hFig, 'Position', [0 0 840 500], 'Color', [1 1 1]);
    
    % The trajectories
    ax = subplot('Position', [0.065 0.6 0.60 0.39]);
    plotTrajectories(ax, fixEMobj, emPosRange, maxDurationSecondsForFixationSpan);
    
    % The velocities
    ax = subplot('Position', [0.065 0.1 0.60 0.39]);
    plotVelocityDistributionPlot(ax, fixEMobj, maxDurationSecondsForFixationSpan)
    
    % The fixation map
    ax = subplot('Position', [0.725 0.6 0.28 0.39]);
    tBins = find(fixEMobj.timeAxis<=maxDurationSecondsForFixationSpan);
    firstPath = squeeze(fixEMobj.emPosArcMin(1, tBins, :));
    plotFixationMap(ax, fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice, firstPath, emPosRange);
    
    % The displacement dynamics
    ax = subplot('Position', [0.725 0.1 0.28 0.39]);
    plotDisplacementPlot(ax,timeLagsMilliseconds, meanD, meanDscrambled, displacementMode);
    
    NicePlot.exportFigToPDF('FixationalEMmodelDemo.pdf', hFig, 300)
end

function plotVelocityDistributionPlot(ax, fixEMobj, maxDurationSecondsForFixationSpan)

    subjectDR.mean = 30;
    subjectDR.std = 17;
    subjectDG.mean = 89;
    subjectDG.std = 63;
    
    allInstancesVelocities = fixEMobj.velocityArcMin(:);
    [NN, edges] = histcounts(allInstancesVelocities, [0:10:250]);
    bar(edges(1:(end-1)),NN, 1, 'FaceAlpha', 0.4, 'FaceColor', [0.8 0.8 0.8]);
    hold on;
    plot(subjectDR.mean*[1 1], [0 max(NN)], 'k-', 'LineWidth', 1.5, 'Color', [1 0.3 0.4]);
    plot((subjectDR.mean-subjectDR.std)*[1 1], [0 max(NN)], 'k--', 'LineWidth', 1.5, 'Color', [1 0.3 0.4]);
    plot(subjectDG.mean*[1 1], [0 max(NN)], 'k-','LineWidth', 1.5, 'Color', [0.2 0.4 1]);
    plot((subjectDG.mean+subjectDG.std)*[1 1], [0 max(NN)], 'k--','LineWidth', 1.5, 'Color', [0.2 0.4 1]);
    
    timeTicks = 0:50:250;
    set(gca, 'XLim', [0 250], ...
        'XTick', timeTicks, ...
        'XTickLabel', sprintf('%2.0f\n', timeTicks), ...
        'YTick', [], ...
        'FontSize', 16);
    xlabel('\it velocity (arc min / sec)');
    box on; grid on
end

function plotDisplacementPlot(ax,timeLagsMilliseconds, meanD, meanDscrambled, displacementMode)
    plot(timeLagsMilliseconds, meanD, 'k-', 'LineWidth', 1.5); hold on;
    plot(timeLagsMilliseconds, meanDscrambled, 'k-.', 'LineWidth', 1.5); hold off
    axis 'square';
    box on;
    if (strcmp(displacementMode, 'D1'))
        yTicks = [0.2 0.5 1 2];
        yLims = [0 2];
    else
        yTicks = [0.1 0.3 1 3 10 30];
        yLims = [0.06  30];
    end
    yTickLabels = sprintf('%2.1f\n', yTicks);
    set(gca, 'XScale', 'log', 'Yscale', 'log');
    set(gca, 'XLim', [1.7 800], 'XTick', [1 3 10 30 100 300 1000], 'YLim', yLims, ...
        'YTick', yTicks, 'YTickLabel', yTickLabels, ...
        'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
    grid on; box on
    xlabel('\it time lag (ms)'); 
    if (strcmp(displacementMode, 'D1'))
        ylabel('\it displacement (arc min)');
    else
        ylabel('\it displacement (arc min^2)');
    end
end

function plotFixationMap(ax, fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice, firstPath, emPosRange)
    
    contourf(fixationMapSupportX, fixationMapSupportY, ...
        fixationMap, 0:0.05:1, 'LineColor', [0.5 0.5 0.5]);
    hold on;
    plot(firstPath(:,1), firstPath(:,2), '-', 'Color', [0 0 0], 'LineWidth', 1.5);
    
    plot([0 0], emPosRange, 'k-');
    plot(emPosRange, [0 0], 'k-');
    plot(fixationMapSupportX, ...
        -emPosRange + fixationMapXSlice * emPosRange * 0.9, '-', ...
        'Color', [1 0 0], 'LineWidth', 1.5);
    plot(emPosRange - fixationMapYSlice * emPosRange * 0.9, ...
        fixationMapSupportY, '-', 'Color', [0 0 1], 'LineWidth', 1.5);
    axis 'square';
    box on;
    grid on
    set(ax, 'YLim', emPosRange*[-1 1], 'XLim', emPosRange*[-1 1], ...
        'XTick', [-20:2:20], 'YTick', [-20:2:20], 'FontSize', 16);
    xlabel('\it position (arc min)');
    ylabel('\it position (arc min)');
    cmap = brewermap(1024, 'Greys');
    colormap(cmap);
    
end

function plotTrajectories(ax, fixEMobj, emPosRange, maxDurationSecondsForFixationSpan)

    TLims = [0 maxDurationSecondsForFixationSpan]*1000;
    allXpos = squeeze(fixEMobj.emPosArcMin(:, :, 1));
    allYpos = squeeze(fixEMobj.emPosArcMin(:, :, 2));
    plot(fixEMobj.timeAxis*1000, allYpos, '-', 'LineWidth', 1.5, 'Color', [0.4 0.4 0.4 0.1]);
    
        
    hold on;
    xPosArcMin = squeeze(fixEMobj.emPosArcMin(1, :, 1));
    yPosArcMin = squeeze(fixEMobj.emPosArcMin(1, :, 2));
    plot(fixEMobj.timeAxis*1000, xPosArcMin, '-', 'Color', [1 0.5 0.5], 'LineWidth', 3);
    plot(fixEMobj.timeAxis*1000, xPosArcMin, '-', 'Color', [1 0 0], 'LineWidth', 1.5);
    plot(fixEMobj.timeAxis*1000, yPosArcMin, '-', 'Color', [0.5 0.5 1], 'LineWidth', 3);
    plot(fixEMobj.timeAxis*1000, yPosArcMin, '-', 'Color', [0 0 1], 'LineWidth', 1.5);
    plot(fixEMobj.timeAxis*1000, yPosArcMin * 0, 'k-');
    hold off;
    grid on;
    box on
    timeTicks = 0:20:300;
    set(gca, 'XLim', TLims , 'XTick', timeTicks, ...
        'YLim', emPosRange*[-1 1], 'YTick', [-20:2:20], 'FontSize', 16);
    xlabel('\it time (msec)');
    ylabel('\it position (arc min)'); 
end
