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
    
    if (~isempty(figNo))
        hFig = figure(figNo); clf;
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 1000]);
    
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', nModels, ...
           'rowsNum', 8, ...
           'heightMargin',   0.03, ...
           'widthMargin',    0.11, ...
           'leftMargin',     0.18, ...
           'rightMargin',    0.03, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.02);
    end
    

    panelFigurePosition = [10 10 200 200];
    
    for iModel = 1:nModels
        model = modelResponses{iModel};

        if (iModel == 1)
            labelYaxis = true;
        else
            labelYaxis = false;
        end
        
        %
        % The pRate response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(1,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            figPos = panelFigurePosition;
            figPos(3) = 260;
            figPos(4) = 220;
            set(hFig, 'Color', [1 1 1], 'Position', figPos);
        end
        
        yLims = [max([0 (floor(min(model.pRate)/1000)-1)*1000]) ceil(max(model.pRate)/1000+1)*1000];
        yRange = (yLims(2)-yLims(1));
        if (yRange > 10000)
            deltaT = 2000;
        elseif (yRange > 5000)
            deltaT = 1000;
        elseif (yRange > 2500)
            deltaT = 500;
        elseif (yRange > 1000)
            deltaT = 500;
        elseif (yRange > 1000)
            deltaT = 200;
        elseif (yRange > 500)
            deltaT = 40;
        else
            deltaT = 20;
        end

        yTicks = 0:deltaT:100000;
        yTickLabels = sprintf('%4.0f\n',yTicks);
        
        plotTemporalResponse(model.timeAxis, model.pRate, 'r', ...
            sprintf('excitation rate\n(photons/cone/sec)'), 'line', true, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        %title(legends{iModel}, 'FontSize', 12);
        if (isempty(figNo))
            NicePlot.exportFigToPDF('pRate_response.pdf', hFig, 300);
        end
        
        %
        % The Opsin response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(2,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            set(hFig, 'Color', [1 1 1], 'Position', panelFigurePosition);
        end
        
        yLims = []; % model.opsin(1)+[0 13];
        yTicks = 0:2:20000;
        yTickLabels = sprintf('%4.0f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.opsin, 'r', ...
            'opsin activation', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        if (isempty(figNo))
            NicePlot.exportFigToPDF('Opsin_response.pdf', hFig, 300);
        end
        
        %
        % The PDE response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(3,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            set(hFig, 'Color', [1 1 1], 'Position', panelFigurePosition);
        end
        
        yLims = []; % model.pde(1)+[-0.05 0.4];
        yTicks = 30:0.1:600;
        yTickLabels = sprintf('%2.1f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.pde, 'r', 'PDE', 'line', false, labelYaxis,...
            yLims, yTicks, yTickLabels);
        if (isempty(figNo))
            NicePlot.exportFigToPDF('PDE_response.pdf', hFig, 300);
        end
        
        %
        % The Ca response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(4,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            set(hFig, 'Color', [1 1 1], 'Position', panelFigurePosition);
        end
        
        yLims = [];
        yTicks = -1:0.05:1;
        yTickLabels = sprintf('%1.3f\n',yTicks);
        
        plotTemporalResponse(model.timeAxis, model.caSlow, 'b', 'Ca', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        hold on;
        plotTemporalResponse(model.timeAxis, model.ca, 'r', 'Ca', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        
        if (isempty(figNo))
            NicePlot.exportFigToPDF('Ca_response.pdf', hFig, 300);
        end
        
        %
        % The GC response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(5,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            set(hFig, 'Color', [1 1 1], 'Position', panelFigurePosition);
        end
        
        yLims = []; % model.gC(1)+[-0.003 7];
        yTicks = 0:1:10000;
        yTickLabels = sprintf('%1.3f\n',yTicks);
        if (~isempty(yLims))&&(yLims(2)==yLims(1))
            yLims = [-1 1];
        end
        plotTemporalResponse(model.timeAxis, model.gC, 'r', 'GC', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        if (isempty(figNo))
            NicePlot.exportFigToPDF('GC_response.pdf', hFig, 300);
        end
        
        
        %
        % The cGMP response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(6,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            set(hFig, 'Color', [1 1 1], 'Position', panelFigurePosition);
        end
        yLims = []; % model.cGMP(1)+[-0.08 0.02];
        yTicks = 13:0.02:21;
        yTickLabels = sprintf('%1.3f\n',yTicks);
        plotTemporalResponse(model.timeAxis, model.cGMP, 'r', 'cGMP', 'line', false, labelYaxis, ...
            yLims, yTicks, yTickLabels);
        if (isempty(figNo))
            NicePlot.exportFigToPDF('cGMP_response.pdf', hFig, 300);
        end
        
        %
        % The photocurrent response
        %
        if (~isempty(figNo))
            subplot('Position', subplotPosVectors(7,iModel).v);
        else
            hFig = figure(); clf;
            figNo = figNo + 1;
            set(hFig, 'Color', [1 1 1], 'Position', panelFigurePosition);
        end
        
        yRange = max(model.membraneCurrent)-min(model.membraneCurrent);
        yRange = max(model.noisyMembraneCurrents(:))-min(model.noisyMembraneCurrents(:));
        
        if yRange > 50
            deltaTick = 10;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.0f\n', yTicks);
        elseif yRange > 20
            deltaTick = 4;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.0f\n', yTicks);
        elseif yRange > 10
            deltaTick = 2;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.0f\n', yTicks);
        elseif yRange > 5
            deltaTick = 1;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.0f\n', yTicks);
        elseif yRange > 2
            deltaTick = 0.4;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.1f\n', yTicks);
        elseif yRange > 1
            deltaTick = 0.2;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.1f\n', yTicks);
        elseif yRange > 0.5
            deltaTick = 0.1;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.1f\n', yTicks);
        else
            deltaTick = 0.02;
            yTicks = -100:deltaTick :0;
            yTickLabels = sprintf('%.2f\n', yTicks);
        end
        
        yLims = [floor(min(model.membraneCurrent(:))/deltaTick)*deltaTick ceil(max(model.membraneCurrent(:))/deltaTick)*deltaTick];
        yLims = [floor(min(model.noisyMembraneCurrents(:))/deltaTick)*deltaTick ceil(max(model.noisyMembraneCurrents(:))/deltaTick)*deltaTick];
        
        plotTemporalResponse(model.timeAxis, model.membraneCurrent, 'r', sprintf('photocurrent (pA)'), 'line', false, labelYaxis, ...
            yLims, yTicks, {});
        if (isempty(figNo))
            NicePlot.exportFigToPDF('pCurrent_response.pdf', hFig, 300);
        end
        
        %
        % The noisy photocurrent
        %
        if (~isempty(model.noisyMembraneCurrents))
            if (~isempty(figNo))
                subplot('Position', subplotPosVectors(8,iModel).v);
            else
                hFig = figure(); clf;
                figNo = figNo + 1;
                figPos = panelFigurePosition;
                figPos(3) = 240;
                figPos(4) = 220;
                set(hFig, 'Color', [1 1 1], 'Position', figPos);
            end
            plotTemporalResponse(model.timeAxis, squeeze(model.noisyMembraneCurrents), 'r', sprintf('photocurrent\n(pAmps)'), 'line', true, labelYaxis, ...
                yLims, yTicks, yTickLabels);
            if (isempty(figNo))
                NicePlot.exportFigToPDF('pCurrentNoisy_response.pdf', hFig, 300);
            end
        end
        
    end
    
end