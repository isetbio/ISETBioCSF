function plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstances, ...
    timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstances, ...
    photocurrentSNR, coneExcitationSNR, adaptationPhotonRates, stepWeberContrasts, legends, figNo)
% Plot the signal-to-noise results
%
% Syntax:
%   plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, stepWeberContrasts, legends, figNo)
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
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotTraces(stepWeberContrasts, adaptationPhotonRates, timeAxis, photoCurrents, noisyPhotoCurrentsInstances, 'photocurrent');
   
    % Plot cone excitations response traces
    figNo = figNo + 1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotTraces(stepWeberContrasts, adaptationPhotonRates, timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstances, 'cone excitations');
   
    % Plot the SNRs across adaptation levels and contrasts
    figNo = figNo + 1000;
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    plotSNRs(stepWeberContrasts, adaptationPhotonRates, coneExcitationSNR, photocurrentSNR);
    
end

function plotSNRs(stepWeberContrasts, adaptationPhotonRates, coneExcitationSNR, photocurrentSNR)

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
    
    cmap = brewermap(nAdaptationLevels, '*spectral');
   
    decrementIndices = find(signedWeberContrasts<0);
    incrementIndices = find(signedWeberContrasts>0);
   
    subplot(2,3,1);
    plotSNR(-signedWeberContrasts(decrementIndices), ...
            squeeze(coneExcitationSNR(:, decrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(cone excitation) (dB)', 'decrement');
   
    subplot(2,3,4);
    plotSNR(signedWeberContrasts(incrementIndices), ...
            squeeze(coneExcitationSNR(:, incrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(cone excitation) (dB)', 'increment');
        
    subplot(2,3,2);
    plotSNR(-signedWeberContrasts(decrementIndices), ...
            squeeze(photocurrentSNR(:, decrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(pCurrent) (dB)', 'decrement');
        
    subplot(2,3,5);
    plotSNR(signedWeberContrasts(incrementIndices), ...
            squeeze(photocurrentSNR(:, incrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(pCurrent) (dB)', 'increment');
        
    subplot(2,3,3);
    SNRRatioDecrements = squeeze(photocurrentSNR(:, decrementIndices)) ./ ...
              squeeze(coneExcitationSNR(:, decrementIndices));
    plotSNR(-signedWeberContrasts(decrementIndices), ...
            SNRRatioDecrements, ...
            cmap, adaptationPhotonRates, 'SNR(pCurrent)/SNR(cone excitation)', 'decrement');
        
   
    subplot(2,3,6);
    SNRRatioIncrements = squeeze(photocurrentSNR(:, incrementIndices)) ./ ...
              squeeze(coneExcitationSNR(:, incrementIndices));
          
    plotSNR(signedWeberContrasts(incrementIndices), ...
            SNRRatioIncrements, ...
            cmap, adaptationPhotonRates, 'SNR(pCurrent)/SNR(cone excitation)', 'increment');

end

function plotTraces(stepWeberContrasts, adaptationPhotonRates, timeAxis, meanTraces, noisyTraces, signalName)

    nAdaptationLevels = size(meanTraces,1);
    mStepStrengths = size(meanTraces,2);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nAdaptationLevels, ...
       'rowsNum', mStepStrengths, ...
       'heightMargin',   0.02, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);
   
    
   for pulseStrengthIndex = 1:mStepStrengths
        
        weberContrastIndex = floor((pulseStrengthIndex-1)/2)+1;
        
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
               title(sprintf('%2.0f photons/cone/sec', stepWeberContrasts(weberContrastIndex)*adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
               yLims = [-90 0];
               yTicks = -90:10:0;
               set(gca, 'YLim',  yLims, 'YTick', yTicks);
           else
               title(sprintf('%2.0f photons/cone/sec', stepWeberContrasts(weberContrastIndex)*adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
               yLims = [0 2000];
               yTicks = 0:500:2000;
               set(gca, 'YLim',  yLims, 'YTick', yTicks);
           end
           
           if (adaptationIndex > 1)
               set(gca, 'YTickLabel', {});
           end
           grid on; box on 
           if (pulseStrengthIndex == 1)
               title(sprintf('%2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
           end
       end
   end
end

function plotSNR(weberC, SNR, cmap, adaptationPhotonRates, theYLabel, theTitle)

   legends = cell(1,numel(adaptationPhotonRates));

   for adaptationIndex = 1:numel(adaptationPhotonRates)
       color = squeeze(cmap(adaptationIndex,:));
       plot(weberC, squeeze(SNR(adaptationIndex,:)), 's-', 'LineWidth', 1.5, 'MarkerSize', 12, ...
           'Color', 0.5*color, 'MarkerFaceColor', 0.5*color+[0.5 0.5 0.5], ...
           'MarkerEdgeColor', 0.5*color);
       if (adaptationIndex == 1)
           hold on;
       end
       legends{adaptationIndex} = sprintf('%2.0f R*/c/sec',adaptationPhotonRates(adaptationIndex));
   end
   
   set(gca, 'FontSize', 14, 'XLim', [1e-2 1],  'XScale', 'log', 'YScale', 'linear');
   xlabel('\it weber contrast');
   ylabel(sprintf('\\it %s', theYLabel));
   grid on; box on;
   legend(legends, 'Location', 'NorthWest');
   title(theTitle)
end
