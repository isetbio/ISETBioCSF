function plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstance, photocurrentSNR, coneExcitationSNR, adaptationPhotonRates, stepWeberContrasts, legends, figNo)
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

    nAdaptationLevels = size(photoCurrents,1);
    mStepStrengths = size(photoCurrents,2);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    
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
        if (mod(pulseStrengthIndex,2) == 1) 
            weberSign = 1;
        else
            weberSign = -1;
        end
        signedWeberContrasts(pulseStrengthIndex) = weberSign * stepWeberContrasts(weberContrastIndex);
             
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
           responseMean = squeeze(photoCurrents(adaptationIndex, pulseStrengthIndex,:));
           
           % A noisy response instance
           responseInstance = squeeze(noisyPhotoCurrentsInstance(adaptationIndex, pulseStrengthIndex,:));

           % Plot mean response and the instance
           subplot('Position', subplotPosVectors(pulseStrengthIndex, adaptationIndex).v);
           plotTemporalResponse(timeAxis, responseMean, 'r', '', 'line', labelXaxis, labelYaxis);
           hold on;
           plotTemporalResponse(timeAxis, responseInstance, 'k', '', 'line', labelXaxis, labelYaxis);
           
           text(150, -5, sprintf('%2.0f photons/cone/sec', stepWeberContrasts(weberContrastIndex)*adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
           yLim = [-90 0];
           set(gca, 'YLim',  yLim, 'YTick', [-90:10:0]);
           if (adaptationIndex > 1)
               set(gca, 'YTickLabel', {});
           end
           grid on; box on 
           if (pulseStrengthIndex == 1)
               title(sprintf('%2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
           end
       end
   end
   
   
   
   % Plot the SNRs across adaptation levels and contrasts
   figNo = figNo + 2000;
   hFig = figure(figNo); clf;
   set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
   cmap = brewermap(nAdaptationLevels, '*spectral');
   
   decrementIndices = find(signedWeberContrasts<0);
   incrementIndices = find(signedWeberContrasts>0);
   
   subplot(2,3,1);
   plotSNR(-signedWeberContrasts(decrementIndices), ...
            squeeze(coneExcitationSNR(:, decrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(cone excitation)', 'decrement');
   
   
   subplot(2,3,4);
   plotSNR(signedWeberContrasts(incrementIndices), ...
            squeeze(coneExcitationSNR(:, incrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(cone excitation)', 'increment');
        
   
   subplot(2,3,2);
   plotSNR(-signedWeberContrasts(decrementIndices), ...
            squeeze(photocurrentSNR(:, decrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(pCurrent)', 'decrement');
        
   
   subplot(2,3,5);
   plotSNR(signedWeberContrasts(incrementIndices), ...
            squeeze(photocurrentSNR(:, incrementIndices)), ...
            cmap, adaptationPhotonRates, 'SNR(pCurrent)', 'increment');
        

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
   
   set(gca, 'FontSize', 14, 'XLim', [1e-2 1], 'YLim', [0.01 1100], 'XScale', 'log', 'YScale', 'log');
   xlabel('\it weber contrast');
   ylabel(sprintf('\\it %s', theYLabel));
   grid on; box on;
   legend(legends, 'Location', 'NorthWest');
   title(theTitle)
end
