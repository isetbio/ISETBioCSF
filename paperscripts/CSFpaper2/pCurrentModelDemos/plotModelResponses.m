function plotModelResponses(modelResponses, legends, figNo)
% Plot the photocurrent model component responses for a number of models
%
% Syntax:
%   plotModelResponses(modelResponses, legends, figNo)
%
% Description:
%    Plot the 6 components of the photocurrent model for a number of models
%    The plots are arranged in a [6 rows x nModels] array of panels
%
% Inputs:
%    modelResponses        - cell array containing the model component
%                            responses for the examined conditions
%    legends               - cell array containing strings with the examined
%                            conditions
% Output:
%    None.
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019

    nModels = length(modelResponses);
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 1000]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nModels, ...
       'rowsNum', 7, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.11, ...
       'leftMargin',     0.18, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);


    for iModel = 1:nModels
        model = modelResponses{iModel};

        if (iModel == 1)
            labelYaxis = true;
        else
            labelYaxis = false;
        end
        
        
        subplot('Position', subplotPosVectors(1,iModel).v);
        yLims = model.pRate(1)+[0 200];
        yTicks = 0:200:100000;
        yTickLabels = sprintf('%4.0f\n',yTicks);
        
        plotTemporalResponse(model.timeAxis, model.pRate, 'k', ...
            sprintf('absorbed photon rate\n(photons/cone/sec)'), 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        title(legends{iModel}, 'FontSize', 12);
        
        subplot('Position', subplotPosVectors(2,iModel).v);
        yLims = []; % model.opsin(1)+[0 13];
        yTicks = 0:2:20000;
        yTickLabels = sprintf('%4.0f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.opsin, 'k', ...
            'opsin activation', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        
        
        subplot('Position', subplotPosVectors(3,iModel).v);
        yLims = []; % model.pde(1)+[-0.05 0.4];
        yTicks = 30:0.1:600;
        yTickLabels = sprintf('%2.1f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.pde, 'k', 'PDE', 'line', false, labelYaxis,...
            yLims, yTicks, yTickLabels);
        
        subplot('Position', subplotPosVectors(4,iModel).v);
        yLims = []; % model.ca(1)+[-0.003 0.001];
        yTicks = 0.004:0.001:1;
        yTickLabels = sprintf('%1.3f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.ca, 'k', '', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        hold on;
        plotTemporalResponse(model.timeAxis, model.caSlow, 'b', 'Ca & slowCa', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        
        subplot('Position', subplotPosVectors(5,iModel).v);
        yLims = []; % model.gC(1)+[-0.003 7];
        yTicks = 0:1:10000;
        yTickLabels = sprintf('%1.3f\n',yTicks);
        if (~isempty(yLims))&&(yLims(2)==yLims(1))
            yLims = [-1 1];
        end
        plotTemporalResponse(model.timeAxis, model.gC, 'k', 'GC', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        
        subplot('Position', subplotPosVectors(6,iModel).v);
        yLims = []; % model.cGMP(1)+[-0.08 0.02];
        yTicks = 13:0.02:21;
        yTickLabels = sprintf('%1.3f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.cGMP, 'k', 'cGMP', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        
        subplot('Position', subplotPosVectors(7,iModel).v);
        yTicks = -100:1:0; % -100:0.05:0;
        yTickLabels = sprintf('%2.2f\n', yTicks);
        yLims = [-45 -40]; % model.membraneCurrent(1)+[-0.3 1]/10;
        plotTemporalResponse(model.timeAxis, model.membraneCurrent, 'k', sprintf('photocurrent\n(pAmps)'), 'line', true, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        
        drawnow
        pause(1)
    end
end