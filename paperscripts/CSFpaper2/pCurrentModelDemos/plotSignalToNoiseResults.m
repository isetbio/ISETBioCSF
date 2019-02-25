function plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstances, ...
    timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstances, ...
    photocurrentSNR, coneExcitationSNR, transformDecibelsToRatios, adaptationPhotonRates, stepWeberContrasts, legends, figNo)
% Plot the signal-to-noise results
%
% Syntax:
%   plotSignalToNoiseResults(timeAxis, stepResponses, adaptationPhotonRates, stepWeberContrasts, legends, figNo)
%
% Description:
%    Plot the photocurrent responses for the different step stimuli which
%    vary in step amplitude and polarity and in the adaptation level.
%    The plots are arranged in a [mStepStrengths x nAdaptationLevels] array of panels
%
% Inputs:
%    timeAxis              - time axis of the response
%    stepResponses         - matrix of step responses
%    adaptationPhotonRates - vector with the examined adaptation levels
%    stepWeberContrasts    - vector with the examined step Weber contrasts 
%    legends               - cell array containing strings with the examined
%                            conditions
%    figNo                 - figure number
%
% Output:
%    None.
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019
    
    
    % Plot photocurrent response traces
    timeLimits = [0 300]; timeTicks = 0:50:1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotTraces(stepWeberContrasts, adaptationPhotonRates, ...
        timeAxis, photoCurrents, noisyPhotoCurrentsInstances, 'photocurrent', ...
        timeLimits, timeTicks);
   
    % Plot cone excitations response traces
    figNo = figNo + 1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotTraces(stepWeberContrasts, adaptationPhotonRates, ...
        timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstances, 'cone excitations', ...
        timeLimits, timeTicks);
   
    % Plot the SNRs across adaptation levels and contrasts
    figNo = figNo + 1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotSNRs(stepWeberContrasts, adaptationPhotonRates, coneExcitationSNR, photocurrentSNR, transformDecibelsToRatios);
    
    figNo = figNo + 1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotSNRsAcrossAllConditions(coneExcitationSNR, photocurrentSNR,  transformDecibelsToRatios)
end

function plotSNRsAcrossAllConditions(coneExcitationSNR, photocurrentSNR,  transformDecibelsToRatios)
    SNRLims = [-37 27];
    SNRTicks = -50:5:50;
    
    if (transformDecibelsToRatios)
       coneExcitationSNR = 10.^(coneExcitationSNR/10);
       photocurrentSNR = 10.^(photocurrentSNR/10);
       SNRLims = 10.^(SNRLims/10);
       SNRTicks = 10.^(SNRTicks/10);
    end
   
    plot(coneExcitationSNR(:), photocurrentSNR(:), 'rp');
    hold on;
    plot(SNRLims, SNRLims, 'k-');
    set(gca, 'FontSize', 14, 'XLim', SNRLims, 'YLim', SNRLims, 'XTick', SNRTicks, 'YTick', SNRTicks);
    
    grid on; box on;
    axis 'square';
    
    if (transformDecibelsToRatios)
        xAxisLabel = 'SNR (cone excitation)';
        yAxisLabel = 'SNR (pCurrent)';
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
    else
        xAxisLabel = 'SNR (cone excitation) (dB)';
        yAxisLabel = 'SNR (pCurrent) (dB)';
    end
    
    xlabel(sprintf('\\it %s', xAxisLabel));
    ylabel(sprintf('\\it %s', yAxisLabel));
    
end

function plotSNRs(stepWeberContrasts, adaptationPhotonRates, coneExcitationSNR, photocurrentSNR, transformDecibelsToRatios)

    nAdaptationLevels = size(coneExcitationSNR,1);
    mStepStrengths = size(coneExcitationSNR,2);
    
    for pulseStrengthIndex = 1:mStepStrengths
       
        weberContrastIndex = floor((pulseStrengthIndex-1)/2)+1;
        if (mod(pulseStrengthIndex,2) == 1) 
            weberSign = 1;
        else
            weberSign = -1;
        end
        signedWeberContrasts(pulseStrengthIndex) = weberSign * stepWeberContrasts(weberContrastIndex);
    end
    
   
    decrementIndices = find(signedWeberContrasts<0);
    incrementIndices = find(signedWeberContrasts>0);
   
    WeberContrastLims = [1e-2 1.1]*100;
    WeberContrastTicks = [0.03 0.1 3 10 30 100];
    SNRLims = [-40 30];
    SNRTicks = -30:10:100;
    SNRDiffLims = [-15 -5];
    SNRDiffTicks = (-30:1:0);
    adaptationRateLims = [30 12000];
    adaptationRateTicks = [30 100 300 1000 3000 10000];
    diffSNR = photocurrentSNR - coneExcitationSNR;
    
    subplot(2,4,1);
    plotSNR(-signedWeberContrasts(decrementIndices)*100, ...
            squeeze(coneExcitationSNR(:, decrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', WeberContrastLims , SNRLims, WeberContrastTicks,...
            SNRTicks, transformDecibelsToRatios, 'contrast (%)', 'SNR (cone excitation)', 'decrement', 'NorthWest');
    
    subplot(2,4,2);
    plotSNR(signedWeberContrasts(incrementIndices)*100, ...
            squeeze(coneExcitationSNR(:, incrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', WeberContrastLims , SNRLims, WeberContrastTicks,...
            SNRTicks, transformDecibelsToRatios, 'contrast (%)', 'SNR (cone excitation)', 'increment', 'NorthWest');
    
    subplot(2,4,5);
    plotSNR(-signedWeberContrasts(decrementIndices)*100, ...
            squeeze(photocurrentSNR(:, decrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', WeberContrastLims , SNRLims, WeberContrastTicks,...
            SNRTicks, transformDecibelsToRatios, 'contrast (%)', 'SNR (pCurrent)', 'decrement', 'NorthWest');
    
    subplot(2,4,6);
    plotSNR(signedWeberContrasts(incrementIndices)*100, ...
            squeeze(photocurrentSNR(:, incrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', WeberContrastLims , SNRLims, WeberContrastTicks,...
            SNRTicks, transformDecibelsToRatios, 'contrast (%)', 'SNR (pCurrent)', 'increment', 'NorthWest');
        
    
    subplot(2,4,3);
    if (transformDecibelsToRatios)
        yAxisLabel = 'SNR(pCurrent) / SNR(cone excitation)';
    else
        yAxisLabel = 'SNR(pCurrent) - SNR(cone excitation) (dB)';
    end
    
    plotSNR(-signedWeberContrasts(decrementIndices)*100, squeeze(diffSNR (:, decrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', WeberContrastLims , SNRDiffLims, WeberContrastTicks,...
            SNRDiffTicks, transformDecibelsToRatios, 'contrast (%)', yAxisLabel, 'decrement', 'SouthWest');   
   
    subplot(2,4,4);
    plotSNR(signedWeberContrasts(incrementIndices)*100, squeeze(diffSNR (:, incrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', WeberContrastLims , SNRDiffLims, WeberContrastTicks,...
            SNRDiffTicks, transformDecibelsToRatios, 'contrast (%)', yAxisLabel, 'increment', 'SouthWest');
        
    
    subplot(2,4,7);
    plotSNR(adaptationPhotonRates, (diffSNR(:,decrementIndices))', ...
            -signedWeberContrasts(decrementIndices)*100, '% contrast', adaptationRateLims , SNRDiffLims, adaptationRateTicks, ...
            SNRDiffTicks, transformDecibelsToRatios, 'adaptation level (R*/c/sec)', yAxisLabel, 'decrement', 'SouthWest');
    
    subplot(2,4,8);
    plotSNR(adaptationPhotonRates, (diffSNR(:,incrementIndices))', ...
            signedWeberContrasts(incrementIndices)*100, '% contrast', adaptationRateLims, SNRDiffLims, adaptationRateTicks, ...
            SNRDiffTicks, transformDecibelsToRatios, 'adaptation level (R*/c/sec)', yAxisLabel, 'increment', 'SouthWest');
        
end

function plotSNR(x, SNR, variedPropertyValue, variedPropertyUnits, XLims, Ylims, XTicks, YTicks, transformDecibelsToRatios, theXLabel, theYLabel,  theTitle, theLegendLocation)

   if (transformDecibelsToRatios)
       SNR = 10.^(SNR/10);
       Ylims = 10.^(Ylims/10);
       YTicks = 10.^(YTicks/10);
   end
   
   legends = cell(1,numel(variedPropertyValue));
   cmap = brewermap(numel(variedPropertyValue), '*spectral');
   cmapSaturation = 0.5;
   
   for kIndex = 1:numel(variedPropertyValue)
       color = squeeze(cmap(kIndex,:));
       plot(x, squeeze(SNR(kIndex,:)), 's-', 'LineWidth', 1.5, 'MarkerSize', 12, ...
           'Color', cmapSaturation*color, 'MarkerFaceColor', cmapSaturation*color + (1-cmapSaturation)*[1 1 1], ...
           'MarkerEdgeColor', cmapSaturation*color);
       if (kIndex == 1)
           hold on;
       end
       legends{kIndex} = sprintf('%2.0f %s',variedPropertyValue(kIndex), variedPropertyUnits);
   end
   
   if (transformDecibelsToRatios)
       YScale = 'log';
   else
       YScale = 'linear';
   end
   
   YTickLabels = {};
   for k = 1:numel(YTicks)
       if (YTicks(k) >= 10)
            YTickLabels{k} = sprintf('%2.0f', YTicks(k));
       elseif (YTicks(k) >= 1)
            YTickLabels{k} = sprintf('%2.1f', YTicks(k));
       elseif (YTicks(k) >= 0.1)
           YTickLabels{k} = sprintf('%2.2f', YTicks(k));
       elseif (YTicks(k) >= 0.01)
           YTickLabels{k} = sprintf('%2.3f', YTicks(k));
       elseif (YTicks(k) >= 0.001)
           YTickLabels{k} = sprintf('%2.4f', YTicks(k));
       else
           YTickLabels{k} = sprintf('%f', YTicks(k));
       end
   end
   
   set(gca, 'FontSize', 14, 'XLim', XLims, 'YLim', Ylims, ...
       'XTick', XTicks, 'YTick', YTicks, 'YTickLabel', YTickLabels, ...
       'XScale', 'log', 'YScale', YScale);
   xlabel(sprintf('\\it %s', theXLabel));
   ylabel(sprintf('\\it %s', theYLabel));
   
   grid on; box on;
   if (numel(legends)<5)
        legend(legends, 'Location', theLegendLocation);
   else
        cmapSaturation = 0.7;
        cmap = bsxfun(@plus, cmapSaturation*cmap, (1-cmapSaturation)*[1 1 1]);
        colormap(cmap);
        normalizedTicks = linspace(0,1,numel(variedPropertyValue));
        cbarHandle = colorbar('Ticks', normalizedTicks, 'TickLabels', sprintf('%2.0f\n', variedPropertyValue));
        cbarHandle.Label.String = sprintf('%s',variedPropertyUnits);
   end
   
   title(theTitle)
end

function plotTraces(stepWeberContrasts, adaptationPhotonRates, timeAxis, meanTraces, noisyTraces, signalName, timeLimits, timeTicks)

    nAdaptationLevels = size(meanTraces,1);
    mStepStrengths = size(meanTraces,2);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nAdaptationLevels, ...
       'rowsNum', mStepStrengths, ...
       'heightMargin',   0.02, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.03);
   
    
   for pulseStrengthIndex = 1:mStepStrengths
        
        weberContrastIndex = floor((pulseStrengthIndex-1)/2)+1;
        
        if (mod(pulseStrengthIndex,2) == 1) 
            signedWeberContrast = stepWeberContrasts(weberContrastIndex);
        else
            signedWeberContrast =-stepWeberContrasts(weberContrastIndex);
        end
        
        labelXaxis = false;
        if (pulseStrengthIndex == mStepStrengths)
           labelXaxis = true;
        end
           
        for adaptationIndex = 1:nAdaptationLevels
           labelYaxis = false;
           if (adaptationIndex == 1)
               labelYaxis = true;
           end
           
           % The mean response
           responseMean = squeeze(meanTraces(adaptationIndex, pulseStrengthIndex,:));
           
           % A noisy response instance
           responseInstance = squeeze(noisyTraces(adaptationIndex, pulseStrengthIndex,:));

           % Plot mean response and the instance
           subplot('Position', subplotPosVectors(pulseStrengthIndex, adaptationIndex).v);
           plotTemporalResponse(timeAxis, responseInstance, 'k', '', 'line', labelXaxis, labelYaxis);
           hold on;
           plotTemporalResponse(timeAxis, responseMean, 'r', '', 'line', labelXaxis, labelYaxis);
           
           if (strcmp(signalName, 'photocurrent')) 
               yLims = [-90 0];
               yTicks = -90:10:0;
               set(gca, 'YLim',  yLims, 'YTick', yTicks);
           else
               yLims = [0 20000];
               yTicks = 0:5000:20000;
               set(gca, 'YLim',  yLims, 'YTick', yTicks, 'YTickLabel', {'0', '', '10k', '', '20k', '', '30k', '', '40k', ''});
           end
           
           set(gca, 'XLim', timeLimits, 'XTick', timeTicks);
           
           stepPhotonRate = signedWeberContrast*adaptationPhotonRates(adaptationIndex);
           if (pulseStrengthIndex == 1)
               title(sprintf('bkgnd:%2.0f photons/cone/sec\n stepDiff: %2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex), stepPhotonRate), 'FontSize', 12);
           else
               title(sprintf('stepDiff: %2.0f photons/cone/sec', stepPhotonRate), 'FontSize', 12);
           end
               
           if (adaptationIndex > 1)
               set(gca, 'YTickLabel', {});
           else
               if (strcmp(signalName, 'photocurrent')) 
                   ylabel('photocurrent (pAmps)')
               else
                   ylabel('cone excitation rate (R*/c/s)')
               end
           end
           grid on; box on 
       end
   end
end



