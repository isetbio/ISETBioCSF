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
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 725]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 1, ...
       'rowsNum', 1, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.15, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.09, ...
       'topMargin',      0.01);
   
   for adaptationIndex = 1:nAdaptationLevels
       dResponses = squeeze(stepResponses(adaptationIndex, :,:));
       response = squeeze(stepResponses(adaptationIndex, 1,:));
       dResponse(adaptationIndex) = ceil(max(abs(dResponses(:)-response(1)))/5)*5;
   end

   legends = {};
   for adaptationIndex = 1:nAdaptationLevels
        subplot('Position', subplotPosVectors(1, 1).v);
        hold on
        if (adaptationIndex == 1)
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'greys');
        elseif (adaptationIndex == 2)
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'blues');
        else 
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'reds');
        end
        
        labelYaxis = false;
        if (adaptationIndex == 1)
          labelYaxis = true;
        end
        
        
        %response = squeeze(stepResponses(adaptationIndex, 1,:));
        %yLims = response(1) + dResponse(adaptationIndex) * [-1 1];
        yLims = [-93 -25]; %[floor(min(stepResponses(:))/5)*5 ceil(max(stepResponses(:))/5)*5]; 
        yTicks = [-100:5:0];
        yTickLabels = sprintf('%2.0f\n', yTicks);
        
        for pulseStrengthIndex = 1:2:mStepStrengths
           labelXaxis = true;
           
           weberContrastIndex = floor(pulseStrengthIndex/2)+1;
           lineColor = squeeze(cmap(weberContrastIndex,:));
           response = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex,:));

           % Plot
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
               %set(gca, 'YTickLabel', {});
           else
              ylabel('\it photocurrent (pA)'); 
           end
           
           set(gca, 'FontSize', 24);
           
           legends{numel(legends)+1} = sprintf('bkgnd: %2.0f R*/c/s; c = %2.1f%%', adaptationPhotonRates(adaptationIndex), stepWeberConstants(weberContrastIndex)*100);
           drawnow
        end
        %legend(legends, 'Location', 'NorthEast');
        
        %title(sprintf('bkgnd: %2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
   end
end
