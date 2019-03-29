function examineSNR()

    pulseDurations = [50 100 200 400]/1000;
    
    if (1==1)
        for pulseDurationIndex = 1:numel(pulseDurations)
            pulseDurationSeconds = pulseDurations(pulseDurationIndex);
            dataFileName = sprintf('results_%dmsec.mat', pulseDurationSeconds*1000);
            load(dataFileName, 'd', 'contrastLevels', 'adaptationLevels');
            for iAdaptationIndex = 1:numel(adaptationLevels)
            for iContrastIndex = 1:numel(contrastLevels)
                theConeExcitationSNR(pulseDurationIndex,iAdaptationIndex, iContrastIndex) = d{iAdaptationIndex, iContrastIndex}.theConeExcitationSNR;
                thePhotoCurrentSNR(pulseDurationIndex,iAdaptationIndex, iContrastIndex) = d{iAdaptationIndex, iContrastIndex}.thePhotoCurrentSNR;
            end
            end
            clear 'd'
        end

        save('theSNRs', 'theConeExcitationSNR', 'thePhotoCurrentSNR', 'contrastLevels', 'adaptationLevels', 'pulseDurations');
    else
        load('theSNRs');
    end
    
    plotSNRsAsAFunctionOfPulseDuration(thePhotoCurrentSNR, 'pCurrent', pulseDurations, contrastLevels, adaptationLevels, 100)
    plotSNRsAsAFunctionOfPulseDuration(theConeExcitationSNR, 'cone exc.', pulseDurations, contrastLevels, adaptationLevels, 200);
    
end

function plotSNRsAsAFunctionOfPulseDuration(theSNRs, signalName, pulseDurations, contrastLevels, adaptationLevels, figNo)
    
    SNRTicks = [0.03 0.1 0.3 1 3 10 30 100 300];
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1200 950]);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', numel(contrastLevels)/2, ...
       'rowsNum', numel(adaptationLevels), ...
       'heightMargin',   0.02, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
   
    markerSize = 10;
    for iContrastIndex = 1:(numel(contrastLevels)/2)
        iContrastIndex2 = find(contrastLevels == -contrastLevels(iContrastIndex));
        contrast = contrastLevels(iContrastIndex2)*100;
        
        if (contrast < 3)
            SNRLims = [0.03 3];
        elseif (contrast < 6)
            SNRLims = [0.05 5];
        elseif (contrast < 12)
            SNRLims = [0.1 10];
        elseif (contrast < 24)
            SNRLims = [0.2 20];
        elseif (contrast < 48)
            SNRLims = [0.4 40];
        else
            SNRLims = [1 100];
        end
        
        for iAdaptationIndex = 1:numel(adaptationLevels)    
            
            background = adaptationLevels(iAdaptationIndex);
            
            subplot('Position', subplotPosVectors(iContrastIndex, iAdaptationIndex).v);
            color = 0.4*[1 1 1];
            plot(pulseDurations, squeeze(theSNRs(:, iAdaptationIndex, iContrastIndex)), 'ko-', ...
                'Color', 0.5*color, 'MarkerFaceColor', color, ...
                'MarkerSize', markerSize, 'LineWidth', 1.5);
            hold on;

            color = 0.8*[1 1 1];
            plot(pulseDurations, theSNRs(:, iAdaptationIndex, iContrastIndex2), 'ko-', ...
                'Color', 0.5*color, 'MarkerFaceColor', color, ...
                'MarkerSize', markerSize, 'LineWidth', 1.5);
            grid on
            set(gca, 'FontSize', 14, 'YLim', SNRLims, 'YTick', SNRTicks, 'YScale', 'log', ...
                'XTick', [100 200 300 400]/1000, 'XTickLabel', [100 200 300 400]/1000, 'XLim', [20 420]/1000);
            if (iContrastIndex == numel(contrastLevels)/2)  
                if (iAdaptationIndex == 1)
                    xlabel('pulse duration (sec)');
                end
            else
                set(gca, 'XTickLabel', {});
            end
            
            if (iAdaptationIndex == 1)
                ylabel(sprintf('SNR (%s)',signalName));
                legend({sprintf('%2.0f%% decr.', contrast), sprintf('%2.0f%% incr.', contrast)}, 'Location', 'NorthWest');
            else
                set(gca, 'YTickLabel', {});
            end
            grid on; box on;
            
            if (iContrastIndex == 1)
                title(sprintf('%2.0f R*/c/s', background), 'FontWeight', 'normal');
            end
        end
    end
    
    
end
