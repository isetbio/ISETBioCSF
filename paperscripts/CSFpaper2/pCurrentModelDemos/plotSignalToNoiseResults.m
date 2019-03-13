function plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstances, ...
    timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstances, ...
    photocurrentSNR, coneExcitationSNR, transformDecibelsToRatios, adaptationPhotonRates, stepWeberContrasts, SNRLims, SNRTicks, legends, figNo)
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
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 1025]);
    plotSNRs(stepWeberContrasts, adaptationPhotonRates, coneExcitationSNR, photocurrentSNR, transformDecibelsToRatios, SNRLims, SNRTicks);
    
    figNo = figNo + 1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 950]);
    plotSNRsAcrossAllConditions(coneExcitationSNR, photocurrentSNR,  transformDecibelsToRatios, SNRLims, SNRTicks)
end

function plotSNRsAcrossAllConditions(coneExcitationSNR, photocurrentSNR,  transformDecibelsToRatios, SNRLims, SNRTicks)
    
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

function plotSNRs(stepWeberContrasts, adaptationPhotonRates, coneExcitationSNR, photocurrentSNR, transformDecibelsToRatios, SNRLims, SNRTicks)

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
   
    WeberContrastLims = [0.02 1.0]*100;
    WeberContrastTicks = [0.03 0.1 3 10 30 100];
    SNRDiffLims = [-18 6];
    SNRDiffTicks = (-30:2:20);
    adaptationRateLims = [40 8000];
    adaptationRateTicks = [30 100 300 1000 3000 10000];
    diffSNR = photocurrentSNR - coneExcitationSNR;
    cMap = brewermap(numel(decrementIndices), '*spectral');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 1, ...
       'rowsNum', 3, ...
       'heightMargin',   0.05, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.12, ...
       'rightMargin',    0.02, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.00);
   
    subplot('Position', subplotPosVectors(1,1).v);
    plotSNR(signedWeberContrasts(incrementIndices)*100, ...
            squeeze(coneExcitationSNR(:, decrementIndices)), ...
            squeeze(coneExcitationSNR(:, incrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', cMap, WeberContrastLims , SNRLims, WeberContrastTicks,...
            SNRTicks, transformDecibelsToRatios, 'contrast (%)', 'SNR (cone excitation)',  'NorthWest');
    
    subplot('Position', subplotPosVectors(2,1).v);
    plotSNR(signedWeberContrasts(incrementIndices)*100, ...
            squeeze(photocurrentSNR(:, decrementIndices)), ...
            squeeze(photocurrentSNR(:, incrementIndices)), ...
            adaptationPhotonRates, 'R*/c/sec', cMap, WeberContrastLims , SNRLims, WeberContrastTicks,...
            SNRTicks, transformDecibelsToRatios, 'contrast (%)',  'SNR (pCurrent)',  'NorthWest');
    
    
    if (transformDecibelsToRatios)
        yAxisLabel = 'SNR(pCurrent) / SNR(cone excitation)';
    else
        yAxisLabel = 'SNR(pCurrent) - SNR(cone excitation) (dB)';
    end  
    cMap = brewermap(numel(decrementIndices), 'greys');

    
    subplot('Position', subplotPosVectors(3,1).v);
    decrementsDiffSNR = (diffSNR(:,decrementIndices))';
    incrementsDiffSNR = (diffSNR(:,incrementIndices))';
    contrasts = signedWeberContrasts(incrementIndices)*100;
    
    % average over contrast
    averageOverContrasts = true;
    if (averageOverContrasts)
        decrementsDiffSNR = mean(decrementsDiffSNR,1);
        dincrementsDiffSNR = mean(incrementsDiffSNR,1);
        contrasts = 100;
    end
    % Decrements
    
    plotSNR(adaptationPhotonRates, decrementsDiffSNR, 'decrements', ...
            contrasts, '% contrast', cMap, adaptationRateLims , SNRDiffLims, adaptationRateTicks, ...
            SNRDiffTicks, transformDecibelsToRatios, 'adaptation level (R*/c/sec)', yAxisLabel,  'SouthWest');
    
    % Increments
    plotSNR(adaptationPhotonRates, 'increments', incrementsDiffSNR, ...
            contrasts, '% contrast', cMap, adaptationRateLims , SNRDiffLims, adaptationRateTicks, ...
            SNRDiffTicks, transformDecibelsToRatios, 'adaptation level (R*/c/sec)', yAxisLabel,  'SouthWest');
    
        
end

function plotSNR(x, decSNR, incSNR, variedPropertyValue, variedPropertyUnits, cMap, XLims, Ylims, XTicks, YTicks, transformDecibelsToRatios, theXLabel, theYLabel,  theLegendLocation)

   if (transformDecibelsToRatios)
       if (~ischar(decSNR))
            decSNR = 10.^(decSNR/10);
       end
       if (~ischar(incSNR))
            incSNR = 10.^(incSNR/10);
       end
       Ylims = 10.^(Ylims/10);
       YTicks = 10.^(YTicks/10);
   end
   
   legends = cell(1,numel(variedPropertyValue));
   
   hold on;
   if (~ischar(decSNR))
       for kIndex = 1:numel(variedPropertyValue)
           color = squeeze(cMap(kIndex,:));
           plot(x, squeeze(decSNR(kIndex,:)), 's-', 'LineWidth', 1.5, 'MarkerSize', 12, ...
               'Color', color*0.5, 'MarkerFaceColor', color, ...
               'MarkerEdgeColor', color*0.5);
           legends{kIndex} = sprintf('%2.0f %s (dec)',variedPropertyValue(kIndex), variedPropertyUnits);
       end
   end
   
   if (~ischar(incSNR))
       no = numel(legends);
       for kIndex = 1:numel(variedPropertyValue)
           color = squeeze(cMap(kIndex,:));
           plot(x, squeeze(incSNR(kIndex,:)), 'o-', 'LineWidth', 1.5, 'MarkerSize', 12, ...
               'Color', color*0.5, 'MarkerFaceColor', color, ...
               'MarkerEdgeColor', color*0.5);
           legends{kIndex+no} = sprintf('%2.0f %s (inc)',variedPropertyValue(kIndex), variedPropertyUnits);
       end
   end
   
   if (transformDecibelsToRatios)
       YScale = 'log';
   else
       YScale = 'linear';
   end
   
   YTickLabels = {};
   if (transformDecibelsToRatios)
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
   else
       for k = 1:numel(YTicks)
            YTickLabels{k} = sprintf('%2.0f', YTicks(k));
       end
   end
   
   set(gca, 'FontSize', 14, 'XLim', XLims, 'YLim', Ylims, ...
       'XTick', XTicks, 'YTick', YTicks, 'YTickLabel', YTickLabels, ...
       'XScale', 'log', 'YScale', YScale);
   xlabel(sprintf('\\it %s', theXLabel));
   ylabel(sprintf('\\it %s', theYLabel));
   
   grid on; box on;
   if (strcmp(variedPropertyUnits, '% contrast'))
%         colormap(cMap);
%         normalizedTicks = linspace(0,1,numel(variedPropertyValue));
%         cbarHandle = colorbar('Location', 'North', 'Ticks', normalizedTicks, 'TickLabels', sprintf('%2.0f\n', variedPropertyValue));
%         cbarHandle.Label.String = sprintf('%s',variedPropertyUnits);
   else
        legend(legends, 'Location', theLegendLocation);  
   end

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
   
   yLims = [-93 -25];
   yTicks = [-100:5:0];
   yTickLabels = sprintf('%2.0f\n', yTicks);
       
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
           plotTemporalResponse(timeAxis, responseInstance, 'k', '', 'line', labelXaxis, labelYaxis, yLims, yTicks, yTickLabels);
           hold on;
           plotTemporalResponse(timeAxis, responseMean, 'r', '', 'line', labelXaxis, labelYaxis, yLims, yTicks, yTickLabels);
           
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



