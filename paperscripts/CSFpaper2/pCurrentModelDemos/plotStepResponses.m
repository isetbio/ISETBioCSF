function plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, stepWeberConstants, legends, figNo)
% Plot the photocurrent responses for the examined step stimuli 
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
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nAdaptationLevels, ...
       'rowsNum', mStepStrengths/2, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
   
   for pulseStrengthIndex = 1:2:mStepStrengths
        labelXaxis = false;
        if (pulseStrengthIndex == mStepStrengths-1)
           labelXaxis = true;
        end
           
        for adaptationIndex = 1:nAdaptationLevels
           labelYaxis = false;
           if (adaptationIndex == 1)
               labelYaxis = true;
           end
           weberContrastIndex = floor(pulseStrengthIndex/2)+1;
           % step increment response
           responseIncrement = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex,:));
           % step decrement response
           responseDecrement = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex+1,:));
           % Plot
           subplot('Position', subplotPosVectors(weberContrastIndex, adaptationIndex).v);
           plotTemporalResponse(timeAxis, responseIncrement, 'r', '', 'line', labelXaxis, labelYaxis);
           hold on;
           plotTemporalResponse(timeAxis, responseDecrement, 'b', '', 'line', labelXaxis, labelYaxis);
           plotTemporalResponse(timeAxis, responseDecrement(1)-(responseDecrement-responseDecrement(1)), 'b', 'current (pAmps)', 'dotted line', labelXaxis, labelYaxis);
           text(150, -5, sprintf('%2.0f photons/cone/sec', stepWeberConstants(weberContrastIndex)*adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
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
end
