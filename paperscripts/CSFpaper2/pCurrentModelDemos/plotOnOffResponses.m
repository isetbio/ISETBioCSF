function plotOnOffResponses(timeAxis, stepResponses, adaptationPhotonRates, stepWeberConstants, legends, figNo)
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
    cmap = brewermap(mStepStrengths/2 , 'spectral');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nAdaptationLevels, ...
       'rowsNum', 1, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
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
           % step increment response
           responseIncrement = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex,:));
           % Plot
           plotTemporalResponse(timeAxis, responseIncrement, lineColor, '', 'line', labelXaxis, labelYaxis);
           if (pulseStrengthIndex==1)
               hold on;
           end
           
           yLim = [-100 0];
           set(gca, 'YLim',  yLim, 'YTick', [-100:10:0]);
           if (adaptationIndex > 1)
               set(gca, 'YTickLabel', {});
           end
           grid on; box on 
           legends{numel(legends)+1} = sprintf('c = %2.1f%%', stepWeberConstants(weberContrastIndex)*100);
        end
        legend(legends);
        
        title(sprintf('adaptation: %2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
   end
end
