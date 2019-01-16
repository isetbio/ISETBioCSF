function varargout = formatFigureForPaper(hFig, varargin)
    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.addParameter('plotRatiosOfOtherConditionsToFirst', false, @islogical);
    p.addParameter('theAxes', []);
    p.addParameter('theAxes2', []);
    p.addParameter('theRatioAxes', []);
    p.addParameter('theRatioLims', []);
    p.addParameter('theRatioTicks', []);
    p.addParameter('theLegend', []);
    p.addParameter('theLegend2', []);
    p.addParameter('theText', []);
    p.addParameter('theText2', []);
    p.addParameter('theTextFontSize', [], @isnumeric);
    p.addParameter('theFigureTitle', '', @ischar);
    p.addParameter('theLegendPosition', [], @(x)(isnumeric(x)||ischar(x)));
    p.addParameter('figureHasFinalSize', false, @islogical);
    
    p.parse(varargin{:});
    
    varargout = {};
    varargout{1} = [];
    varargout{2} = [];
    
    figureType = p.Results.figureType;
    plotRatiosOfOtherConditionsToFirst = p.Results.plotRatiosOfOtherConditionsToFirst;
    
    theAxes =  p.Results.theAxes;
    theAxes2 =  p.Results.theAxes2;
    theRatioAxes = p.Results.theRatioAxes;
    theRatioLims = p.Results.theRatioLims;
    theRatioTicks = p.Results.theRatioTicks;
    theLegend = p.Results.theLegend;
    theLegend2 = p.Results.theLegend2;
    theText = p.Results.theText;
    theText2 = p.Results.theText2;
    theTextFontSize = p.Results.theTextFontSize;
    theFigureTitle = p.Results.theFigureTitle;
    theLegendPosition = p.Results.theLegendPosition;
    figureHasFinalSize = p.Results.figureHasFinalSize;
    
    switch (figureType)
        case 'CONE_APERTURE'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);
            else
                axis(theAxes, 'square');
                set(theAxes, 'FontSize', 22, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                ytickformat(theAxes, '%.1f');
                xtickformat(theAxes, '%.0f');
                if (~isempty(theText))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
            end
        case 'ABERRATION_MAP_PSF_COMBO'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 1050 500], 'Color', [1 1 1]);
            else
                axis(theAxes, 'square');
                axis(theAxes, 'xy');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                ytickformat(theAxes, '%.1f');
                xtickformat(theAxes, '%.1f');
                if (~isempty(theText))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
            end
            
        case 'RESPONSE_INSTANCE_SINGLE_CONDIITION'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 900 500], 'Color', [1 1 1]);
            else
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75, 'Ydir', 'Normal');
                box(theAxes, 'on');
                grid(theAxes, 'on');
                
                if (~isempty(theLegend))
                    set(theLegend,  'Location', 'NorthEast');
                end
            end
            
        case 'RESPONSE_INSTANCE_TWO_CONDIITIONS'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 900 860], 'Color', [1 1 1]);
            else
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75, 'Ydir', 'Normal');
                box(theAxes, 'on');
                grid(theAxes, 'on');
            end
            
        case 'RESPONSE_INSTANCE_THREE_CONDIITIONS'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 900 1500], 'Color', [1 1 1]);
            else
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75, 'Ydir', 'Normal');
                box(theAxes, 'on');
                grid(theAxes, 'on');
                xtickformat(theAxes, '%.0f');
                
                if (~isempty(theLegend))
                    set(theLegend,  'Location', 'NorthEast');
                end
                if (~isempty(theText))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
                
            end
            
        case {'PIGMENT_QUANTAL_EFFICIENCY', 'PIGMENT_TRANSMITTANCE'}
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);
            else
                axis(theAxes, 'square');
                set(theAxes, 'FontSize', 22, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                ytickformat(theAxes, '%.1f');
                xtickformat(theAxes, '%.0f');
                if (~isempty(theLegend))
                    set(theLegend,  'Location', 'SouthEast');
                end
                if (~isempty(theText))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
            end
            
        case 'DISPLAY_PROPERTIES'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 1300 950], 'Color', [1 1 1]);
            else
                axis(theAxes, 'square');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                ytickformat(theAxes, '%.1f');
                xtickformat(theAxes, '%.1f');
                if (~isempty(theLegend))
                    set(theLegend,  'Location', 'SouthEast');
                end
                if (~isempty(theText))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
            
            end
            
        case 'PSYCHOMETRIC_FUNCTIONS_2_CLASSIFIERS'
            if (isempty(theAxes)) 
                set(hFig, 'Position', [10 10 890 450], 'Color', [1 1 1]);
                varargout{1} = subplot('Position', [0.08 0.13 0.42 0.86]);
                varargout{2} = subplot('Position', [0.56 0.13 0.42 0.86]);     
                
            else
                %axis(theAxes, 'square');
                axis(theAxes, 'xy');
                set(theAxes, 'XScale', 'log', 'FontSize', 22, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                xtickformat(theAxes, '%.3f'); ytickformat(theAxes, '%.1f');
                
                %axis(theAxes2, 'square');
                axis(theAxes2, 'xy');
                set(theAxes2, 'XScale', 'log', 'FontSize', 22, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes2, 'on');
                grid(theAxes2, 'on');
                xtickformat(theAxes2, '%.3f'); ytickformat(theAxes2, '%.1f');
                set(theAxes2, 'YTickLabel', {});
                
                if (~isempty(theLegend))
                    set(theLegend,  'Location', 'NorthWest');
                end
            
                if (~isempty(theLegend2))
                    set(theLegend2,  'Location', 'NorthWest');
                end
                
                if (~isempty(theText))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
                
                if (~isempty(theText2))
                    if (isempty(theTextFontSize))
                        theTextFontSize = 30;
                    end
                    set(theText2, 'FontSize', theTextFontSize, ...
                        'FontWeight', 'Bold', ...
                        'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [ 0 0 0], ...
                        'LineWidth', 1.0);
                end
            end
            
        case 'PSYCHOMETRIC_FUNCTIONS'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 1000 450], 'Color', [1 1 1]);
                varargout{1} = subplot('Position', [0.06 0.13 0.44 0.78]);
                varargout{2} = subplot('Position', [0.55 0.13 0.44 0.78]);     
            else
                axis(theAxes, 'square');
                axis(theAxes, 'xy');
                set(theAxes, 'XScale', 'log', 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                xtickformat(theAxes, '%.3f'); ytickformat(theAxes, '%.2f');
                
                if (~isempty(theLegend))
                    set(theLegend,  'Location', 'SouthEast');
                end
            
                axis(theRatioAxes, 'square');
                axis(theRatioAxes, 'xy');
                set(theRatioAxes, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theRatioAxes, 'on');
                grid(theRatioAxes, 'on');
                xtickformat(theRatioAxes, '%3.0g'); ytickformat(theRatioAxes, '%.0f'); 
            end
            
        case 'CONE_SPATIAL_POOLING'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 1400 500], 'Color', [1 1 1]);
            else
                axis(theAxes, 'equal');
                axis(theAxes, 'xy');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                xtickformat(theAxes, '%.1f'); ytickformat(theAxes, '%.1f');
            end
        case {'STIMULUS_AND_OPTICAL_IMAGE', 'STIMULUS_OPTICAL_IMAGE_ISOMERIZATIONS_PHOTOCURRENTS'}
            if (isempty(theAxes))
                if (strcmp(figureType, 'STIMULUS_AND_OPTICAL_IMAGE'))
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 900 400]); 
                else
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1400 450]);
                end
            else
                axis(theAxes, 'square');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
                xtickformat(theAxes, '%.2f'); ytickformat(theAxes, '%.2f');
            end
            
        case {'MOSAICS', 'MOSAICS2', 'PSFS', 'PSFS2', 'OTFS'}
            if (isempty(theAxes))
                if (strcmp(figureType, 'MOSAICS'))
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 840 900]);
                elseif (strcmp(figureType, 'MOSAICS2'))
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 900 500]);
                elseif (strcmp(figureType, 'PSFS'))
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 925]);
                elseif (strcmp(figureType, 'OTFS'))
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 925]);
                elseif (strcmp(figureType, 'PSFS2'))
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 450 830]);
                end
            else
                if (strcmp(figureType, 'MOSAICS')) || (strcmp(figureType, 'MOSAICS2'))
                    %axis(theAxes, 'equal');
                    xtickformat('%0.1f');
                    ytickformat('%0.1f');
                else
                    axis(theAxes, 'square');
                end
                
                set(theAxes, 'FontSize', 22, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes, 'on');
                grid(theAxes, 'on');
            end

            if (~isempty(theText))
                if (isempty(theTextFontSize))
                    theTextFontSize = 30;
                end
                set(theText, 'FontSize', theTextFontSize, ...
                    'FontWeight', 'Bold', ...
                    'BackgroundColor', [1 1 1], ...
                    'EdgeColor', [ 0 0 0], ...
                    'LineWidth', 1.0);
            end
            
        case 'CSF'
            csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
            csLims = [1.5 15000];
            sfTicks = [1 2 5 10 20 50 100];
            sfLims  = [1.5 80];
            dx1 = 0.2; dx2 = 10; dy1 = 0.5; dy2 = 1000;
                
            if (isempty(theAxes)) && (isempty(theLegend))
                if (plotRatiosOfOtherConditionsToFirst)
                    if (figureHasFinalSize)
                        set(hFig, 'Color', [1 1 1], 'Position', [10 10 190 340]);
                        varargout{1} = subplot('Position', [0.215 0.31 0.78 0.68]);
                        varargout{2} = subplot('Position', [0.215 0.08 0.78 0.22]);
                    else
                        set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 900]);
                        varargout{1} = subplot('Position', [0.19 0.30 0.80 0.69]);
                        varargout{2} = subplot('Position', [0.19 0.07 0.80 0.22]);
                    end
                    
                else
                    if (figureHasFinalSize)
                        set(hFig, 'Color', [1 1 1], 'Position', [10 10 190 265]);
                    else
                        set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 700]);
                    end
                    varargout{1} = subplot('Position', [0.19 0.09 0.80 0.90]);
                end
            end

            if (figureHasFinalSize)
                fontSize = 9;
                legendFontSize = 7;
                textFontSize = 12;
                labelFontSize = 10;
                axesLineWidth = 0.25;
            else
                fontSize = 22;
                legendFontSize = 18;
                textFontSize = 30;
                labelFontSize = 24;
                axesLineWidth = 0.25;
            end
                    
            if (~isempty(theAxes))
                set(theAxes, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, ...
                    'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', [csLims(1)-dy1 csLims(2)+dy2]);
                box(theAxes, 'on'); grid(theAxes, 'on');
                set(theAxes, 'FontSize', fontSize, 'TickLength',[0.02, 0.02], 'LineWidth', axesLineWidth);
                xlabel(theAxes,'\it spatial frequency (c/deg)', 'FontWeight', 'normal', 'FontSize', labelFontSize);
                ylabel(theAxes,'\it contrast sensitivity', 'FontWeight', 'normal', 'FontSize', labelFontSize);
            end
            
            if (~isempty(theLegend))
                if (~isempty(theLegendPosition))
                    if (ischar(theLegendPosition))
                        set(theLegend, 'FontSize', legendFontSize, 'Location', theLegendPosition);
                    else
                        set(theLegend, 'FontSize', legendFontSize, 'Position', theLegendPosition, 'Units', 'normalized');
                    end
                else
                    set(theLegend, 'FontSize', legendFontSize, 'Location', 'NorthEast');
                end
            end
            
            if (~isempty(theText))
                if (isempty(theTextFontSize))
                    theTextFontSize = textFontSize;
                end
                set(theText, 'FontSize', theTextFontSize, ...
                    'FontWeight', 'Bold', ...
                    'BackgroundColor', [1 1 1], ...
                    'EdgeColor', [ 0 0 0], ...
                    'LineWidth', 1.0);
            end
            
            if (plotRatiosOfOtherConditionsToFirst)
                if (~isempty(theRatioAxes))
                    set(theRatioAxes, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, ...
                        'XLim', [sfLims(1)-dx1 sfLims(2)+dx2]);
                    
                    if (~isempty(theRatioLims))
                        set(theRatioAxes, 'YLim', theRatioLims);
                    end
                    if (~isempty(theRatioTicks))
                        set(theRatioAxes, 'YTick', theRatioTicks);
                    end
                    ytickformat(theRatioAxes,'%.2f')
                    box(theRatioAxes, 'on'); grid(theRatioAxes, 'on');
                    set(theRatioAxes, 'FontSize', fontSize, 'TickLength',[0.02, 0.02] , 'LineWidth', axesLineWidth);
                    xlabel(theRatioAxes,'\it spatial frequency (c/deg)', 'FontWeight', 'normal', 'FontSize', labelFontSize);
                    ylabel(theRatioAxes,'\it sensitivity ratio', 'FontWeight', 'normal', 'FontSize', labelFontSize);
                    set(theAxes, 'XTickLabel', {});
                end
            end
            
        otherwise
            error('Unknown figure type: ''%s''.', figureType);
    end % switch
    
    if (figureHasFinalSize)
        set(theAxes,'GridLineStyle','-');
        set(theAxes,'GridAlpha', 0.15);
        set(theAxes,'MinorGridLineStyle','-')
        set(theAxes,'MinorGridAlpha', 0.1);
        
        if (plotRatiosOfOtherConditionsToFirst)
            if (~isempty(theRatioAxes))
                set(theRatioAxes,'GridLineStyle','-');
                set(theRatioAxes,'GridAlpha', 0.15);
                set(theRatioAxes,'MinorGridLineStyle','-')
                set(theRatioAxes,'MinorGridAlpha', 0.1);
            end
        end
    end
    if (~isempty(theFigureTitle))
        title(theAxes, theFigureTitle);
    end
                
end

