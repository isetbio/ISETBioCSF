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
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1220 845]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nModels, ...
       'rowsNum', 7, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.06, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.01, ...
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
        plotTemporalResponse(model.timeAxis, model.pRate, 'r', sprintf('absorbed photon rate\n(photons/cone/sec)'), 'line', false, labelYaxis);
        yLim = model.pRate(1)+[0 1100];
        yTicks = 0:200:100000;
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%4.0f\n',yTicks));
        grid on; box on;
        title(legends{iModel}, 'FontSize', 12);
        
        subplot('Position', subplotPosVectors(2,iModel).v);
        plotTemporalResponse(model.timeAxis, model.opsin, 'r', 'opsin activation', 'line', false, labelYaxis);
        yLim = model.opsin(1)+[0 13];
        yTicks = 0:2:20000;
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%4.0f\n',yTicks));
        grid on; box on;
        
        subplot('Position', subplotPosVectors(3,iModel).v);
        plotTemporalResponse(model.timeAxis, model.pde, 'r', 'PDE', 'line', false, labelYaxis);
        yLim = model.pde(1)+[-0.05 0.4];
        yTicks = 30:0.1:600;
        set(gca, 'YLim', yLim, 'YTick', yTicks);
        grid on; box on;
        
        subplot('Position', subplotPosVectors(4,iModel).v);
        plotTemporalResponse(model.timeAxis, model.ca, 'r', '', 'line', false, labelYaxis);
        hold on;
        plotTemporalResponse(model.timeAxis, model.caSlow, 'b', 'Ca & slowCa', 'line', false, labelYaxis);
        yLim = model.ca(1)+[-0.003 0.001];
        yTicks = 0.004:0.001:1;
        set(gca, 'YLim', yLim, 'YTick', yTicks);
        grid on; box on;
        
        subplot('Position', subplotPosVectors(5,iModel).v);
        plotTemporalResponse(model.timeAxis, model.gC, 'r', 'GC', 'line', false, labelYaxis);
        yLim = [min(model.gC(1:end-10)) max(model.gC(1:end-10))];
        yTicks = 0:1:10000;
        set(gca, 'YLim', yLim, 'YTick', yTicks);
        grid on; box on;
        
        subplot('Position', subplotPosVectors(6,iModel).v);
        plotTemporalResponse(model.timeAxis, model.cGMP, 'r', 'cGMP', 'line', false, labelYaxis);
        yTicks = 13:0.02:21;
        yLim = model.cGMP(1)+[-0.08 0.02];
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%2.2f\n', yTicks));
        grid on; box on;
        
        subplot('Position', subplotPosVectors(7,iModel).v);
        plotTemporalResponse(model.timeAxis, model.membraneCurrent, 'r', sprintf('photocurrent\n(pAmps)'), 'line', true, labelYaxis);
        yTicks = -90:0.2:0;
        yLim = model.membraneCurrent(1)+[-0.3 1];
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%2.1f\n', yTicks));
        grid on; box on;
    end
end