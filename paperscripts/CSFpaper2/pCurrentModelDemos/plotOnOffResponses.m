function plotOnOffResponses(timeAxis, stepResponses, modelResponses, adaptationPhotonRates, stepWeberConstants, legends, figNo)
% Plot the photocurrent responses for the examined On/OFF stimuli 
%
% Syntax:
%   plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, stepWeberConstants, legends, figNo)
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
%    stepWeberConstants    - vector with the examined step Weber contrasts 
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

    nAdaptationLevels = size(stepResponses,1);
    mStepStrengths = size(stepResponses,2);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    
    % Color scheme to use for the different impulse resposes
    cmap = brewermap(mStepStrengths/2 , 'reds');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nAdaptationLevels, ...
       'rowsNum', 1, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.1, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.10, ...
       'topMargin',      0.02);
   
   for adaptationIndex = 1:nAdaptationLevels
        subplot('Position', subplotPosVectors(1, adaptationIndex).v);
        
        labelYaxis = false;
        if (adaptationIndex == 1)
          labelYaxis = true;
        end
        
        legends = {};
        for pulseStrengthIndex = 1:2:mStepStrengths
           labelXaxis = true;
           
           weberContrastIndex = floor(pulseStrengthIndex/2)+1;
           lineColor = squeeze(cmap(weberContrastIndex,:));
           response = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex,:));

           % Plot
           yLims = [-90 -25];
           yTicks = [-100:5:0];
           
           yTickLabels = sprintf('%2.0f\n', yTicks);
           if (~isempty(modelResponses{adaptationIndex,pulseStrengthIndex}.noisyMembraneCurrents))
                noisyResponse = squeeze(modelResponses{adaptationIndex,pulseStrengthIndex}.noisyMembraneCurrents(1,:));
                plotTemporalResponse(timeAxis, noisyResponse, lineColor, ...
               '', 'line', labelXaxis, labelYaxis, ...
               yLims, yTicks, yTickLabels);  
           end
           hold on;
           plotTemporalResponse(timeAxis, response, lineColor, ...
               '', 'line', labelXaxis, labelYaxis, ...
               yLims, yTicks, yTickLabels);
           
           if (adaptationIndex > 1)
               set(gca, 'YTickLabel', {});
           else
              ylabel('photocurrent (pAmps)'); 
           end
           legends{numel(legends)+1} = sprintf('c = %2.1f%%', stepWeberConstants(weberContrastIndex)*100);
           drawnow
        end
        legend(legends, 'Location', 'SouthWest');
        
        title(sprintf('bkgnd: %2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
   end
end
