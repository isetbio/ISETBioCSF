function [emPosArcMin, timeAxis] = generateFEM(nTrials, emDurationSeconds)
    sampleTimeSeconds = 1.0 / 1000;

    microSaccadeType = 'none'; % , 'heatmap/fixation based', 'none'
    fixEMobj = fixationalEM();
    fixEMobj.randomSeed = 1;
    fixEMobj.microSaccadeType = microSaccadeType;
    
    computeVelocity = ~true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
            nTrials, computeVelocity, 'useParfor', true);
        
    timeAxis = fixEMobj.timeAxis;
    emPosArcMin = fixEMobj.emPosArcMin;
    
    % Center at (0,0)
    meanEMPos = mean(emPosArcMin,2);
    for iTrial = 1:nTrials
        emPosArcMin(iTrial,:,1) = emPosArcMin(iTrial,:,1) - meanEMPos(iTrial,1,1);
        emPosArcMin(iTrial,:,2) = emPosArcMin(iTrial,:,2) - meanEMPos(iTrial,1,2);
    end

    % Plot the eye movements
    plotFEM(emPosArcMin, timeAxis);

end

function plotFEM(emPosArcMin, timeAxis)
    trialNo = 1;
    xPosArcMin = squeeze(emPosArcMin(trialNo,:,1));
    yPosArcMin = squeeze(emPosArcMin(trialNo,:,2));
    
    maxPos = 20; % max(abs(emPosArcMin(:)));
    if (maxPos == 0)
        maxPos = 1;
    end
    posRange = maxPos*[-1 1];
    
    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 900 400], 'Color', [1 1 1]);
    subplot(1,2,1);
    plot(timeAxis*1000, xPosArcMin, 'r-', 'LineWidth', 1.5); hold on
    plot(timeAxis*1000, yPosArcMin, 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', posRange, 'FontSize', 14, 'XTick', 0:200:1000, 'YTick', -20:5:20);
    grid on
    xlabel('time (msec)');
    ylabel('space (arc min)');
    axis 'square';
    
    subplot(1,2,2);
    plot(xPosArcMin, yPosArcMin, 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', posRange, 'YLim', posRange, 'FontSize', 14, 'XTick', -20:5:20, 'YTick', -20:5:20);
    grid on
    axis 'square';
    xlabel('space (arc min)');
end