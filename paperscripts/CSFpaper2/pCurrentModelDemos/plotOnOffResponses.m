function plotOnOffResponses(timeAxis, stepResponses, modelResponses, adaptationPhotonRates, stepWeberConstants, negativePulseDelay, figNo)
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
    mWeberContrasts = mStepStrengths /2;

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 2, ...
       'rowsNum', 1, ...
       'heightMargin',   0.02, ...
       'widthMargin',    0.09, ...
       'leftMargin',     0.08, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.12, ...
       'topMargin',      0.01);
   
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1000 725]);
    subplot('Position', subplotPosVectors(1, 1).v);
    
   for adaptationIndex = 1:nAdaptationLevels
       dResponses = squeeze(stepResponses(adaptationIndex, :,:));
       response = squeeze(stepResponses(adaptationIndex, 1,:));
       dResponse(adaptationIndex) = ceil(max(abs(dResponses(:)-response(1)))/5)*5;
   end

   legends = {};
   for adaptationIndex = 1:nAdaptationLevels
        
        hold on
        if (adaptationIndex == 1)
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'greys');
        elseif (adaptationIndex == 2)
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'blues');
        elseif (adaptationIndex == 3)
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'greens');
        elseif (adaptationIndex == 4)
            % Color scheme to use for the different impulse resposes
            cmap = brewermap(mStepStrengths/2 , 'reds');
        end
        
        labelYaxis = false;
        if (adaptationIndex == 1)
          labelYaxis = true;
        end
        
        
        %response = squeeze(stepResponses(adaptationIndex, 1,:));
        %yLims = response(1) + dResponse(adaptationIndex) * [-1 1];
        yLims = [-87 -5]; %[floor(min(stepResponses(:))/5)*5 ceil(max(stepResponses(:))/5)*5]; 
        yTicks = [-100:10:0];
        yTickLabels = sprintf('%2.0f\n', yTicks);
        
        for pulseStrengthIndex = 1:2:mStepStrengths
           labelXaxis = true;
           
           weberContrastIndex = floor(pulseStrengthIndex/2)+1;
           lineColor = squeeze(cmap(weberContrastIndex,:));
           response = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex,:));
           
           [posResponseAmplitude(adaptationIndex,weberContrastIndex), ...
            posResponseTimeToPeak(adaptationIndex,weberContrastIndex), ...
            negResponseAmplitude(adaptationIndex,weberContrastIndex), ...
            negResponseTimeToPeak(adaptationIndex,weberContrastIndex)] = ...
               analyzePhotoCurrentWaveform(response, timeAxis, negativePulseDelay);
               
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
           
           if (pulseStrengthIndex == 1)
                legends{numel(legends)+1} = sprintf('bkgnd: %2.0f R*/c/s', adaptationPhotonRates(adaptationIndex));
           end
           
           drawnow
        end
        
        ylabel('\it photocurrent (pA)'); 
        set(gca, 'FontSize', 20);
   end
   
   subplot('Position', subplotPosVectors(1, 2).v);
   hold on;
   
   
   if (mWeberContrasts < nAdaptationLevels)
       legends = {};
       cmap = brewermap(mWeberContrasts*2 , '*RdYlBu');
       
       for weberContrastIndex = 1:mWeberContrasts
           lineColor = squeeze(cmap(weberContrastIndex,:));
           faceColor = lineColor/2+0.5*[1 1 1];
           plot(adaptationPhotonRates, squeeze(negResponseAmplitude(:,weberContrastIndex)), '--o', ...
              'LineWidth', 1.5, 'Color', lineColor, ...
               'MarkerFaceColor', faceColor, 'MarkerSize', 12);
           legends{numel(legends)+1} = sprintf('c: %2.2f', -stepWeberConstants(weberContrastIndex));
       end
       
       for weberContrastIndex = 1:mWeberContrasts
           lineColor = squeeze(cmap(mWeberContrasts+weberContrastIndex,:));
           faceColor = lineColor/2+0.5*[1 1 1];
           plot(adaptationPhotonRates, squeeze(posResponseAmplitude(:,weberContrastIndex)), '-o', ...
              'LineWidth', 1.5, 'Color', lineColor, ...
               'MarkerFaceColor', faceColor, 'MarkerSize', 12);
           legends{numel(legends)+1} = sprintf('c: %2.2f', stepWeberConstants(weberContrastIndex));
       end
       
       xlabel('\it background (R*/c/s)');
       ylabel('\it photocurrent differential (pA)');
       set(gca, 'FontSize', 20, 'XLim', [adaptationPhotonRates(1)*0.9 adaptationPhotonRates(end)*1.1], ...
           'YTick', [-100:5:100], 'YLim', [-50 20], 'XScale', 'log');
       legend(legends, 'Location', 'NorthWest');
       grid on; box on;
   else 
       for adaptationIndex = 1:nAdaptationLevels

            if (adaptationIndex == 1)
                lineColor = [0 0 0];
                faceColor = [0.8 0.8 0.8];
            elseif (adaptationIndex == 2)
                lineColor = [0 0 1];
                faceColor = [0.5 0.5 1.0];
            elseif (adaptationIndex == 3)
                lineColor = [0 0.8 0];
                faceColor = [0.5 0.8 0.5];
            elseif (adaptationIndex == 4)
                lineColor = [1 0 0];
                faceColor = [1.0 0.5 0.5];
            end
            plot(stepWeberConstants, posResponseAmplitude(adaptationIndex,:), ...
                'o-', 'LineWidth', 1.5, 'Color', lineColor, ...
                'MarkerFaceColor', faceColor, 'MarkerSize', 12);
       end

       for adaptationIndex = 1:nAdaptationLevels
            if (adaptationIndex == 1)
                lineColor = [0 0 0];
                faceColor = [0.5 0.5 0.5];
            elseif (adaptationIndex == 2)
                lineColor = [0 0 1];
                faceColor = [0.5 0.5 1.0];
            elseif (adaptationIndex == 3)
                lineColor = [0 0.8 0];
                faceColor = [0.5 0.8 0.5];
            elseif (adaptationIndex == 4)
                lineColor = [1 0 0];
                faceColor = [1.0 0.5 0.5];
            end
            plot(stepWeberConstants, negResponseAmplitude(adaptationIndex,:), '--o', ...
                'LineWidth', 1.5, 'Color', lineColor, ...
                'MarkerFaceColor', faceColor, 'MarkerSize', 12);
       end
       xlabel('\it contrast');
       ylabel('\it photocurrent differential (pA)');
       set(gca, 'FontSize', 20, 'XLim', [0.045 1.0], 'XTick', [0.05 0.1 0.2 0.4 0.8], 'YTick', [-100:5:100], 'YLim', [-50 20], 'XScale', 'log');
       legend(legends, 'Location', 'NorthWest');
       grid on; box on;
   end
   
   
end

function [posResponseAmplitude, posResponseTimeToPeak, ...
            negResponseAmplitude, negResponseTimeToPeak] = ...
               analyzePhotoCurrentWaveform(response, timeAxis, negativePulseDelay)
           
    [posResponseAmplitude, idx] = max(response - response(1));
    posResponseTimeToPeak = timeAxis(idx);
    
    [negResponseAmplitude, idx] = min(response - response(1));
    negResponseAmplitude = negResponseAmplitude;
    negResponseTimeToPeak = timeAxis(idx)-negativePulseDelay;
           
end

           