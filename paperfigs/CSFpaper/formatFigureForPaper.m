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
    p.addParameter('theTextFontSize', [], @isnumeric);
    p.addParameter('theFigureTitle', '', @ischar);
    
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
    theTextFontSize = p.Results.theTextFontSize;
    theFigureTitle = p.Results.theFigureTitle;
    
    switch (figureType)
        case 'CONE_APERTURE'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);
            else
                axis(theAxes, 'square');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
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
                xtickformat(theAxes, '%.0f');
                
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
            
        case 'PIGMENT_TRANSMITTANCE'
            if (isempty(theAxes))
                set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);
            else
                axis(theAxes, 'square');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
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
                
                axis(theAxes2, 'square');
                axis(theAxes2, 'xy');
                set(theAxes2, 'XScale', 'log', 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                box(theAxes2, 'on');
                grid(theAxes2, 'on');
                xtickformat(theAxes2, '%.3f'); ytickformat(theAxes2, '%.2f');
                
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
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 900]);
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
            
        case {'MOSAICS', 'PSFS'}
            if (isempty(theAxes))
                set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 900]);
            else
                if (strcmp(figureType, 'MOSAICS'))
                    axis(theAxes, 'equal');
                else
                    axis(theAxes, 'square');
                end
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
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
            sfLims  = [1.2 100];
            dx1 = 0.2; dx2 = 20; dy1 = 0.4; dy2 = 1000;
                
            if (isempty(theAxes)) && (isempty(theLegend))
                if (plotRatiosOfOtherConditionsToFirst)
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 900]);
                    varargout{1} = subplot('Position', [0.16 0.30 0.83 0.69]);
                    varargout{2} = subplot('Position', [0.16 0.07 0.83 0.22]);
                else
                    set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 700]);
                    varargout{1} = subplot('Position', [0.16 0.09 0.83 0.90]);
                end
            end

            if (~isempty(theAxes))
                set(theAxes, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, ...
                    'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', [csLims(1)-dy1 csLims(2)+dy2]);
                box(theAxes, 'on'); grid(theAxes, 'on');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                xlabel(theAxes,'spatial frequency (c/deg)', 'FontWeight', 'bold', 'FontSize', 20);
                ylabel(theAxes,'contrast sensitivity', 'FontWeight', 'bold', 'FontSize', 20);
            end
            
            if (~isempty(theLegend))
                set(theLegend,  'FontSize', 16);
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
                    set(theRatioAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02], 'LineWidth', 0.75);
                    xlabel(theRatioAxes,'spatial frequency (c/deg)', 'FontWeight', 'bold', 'FontSize', 20);
                    ylabel(theRatioAxes,'sensitivity ratio', 'FontWeight', 'bold', 'FontSize', 20);
                    set(theAxes, 'XTickLabel', {});
                end
            end
    end % switch
    
    if (~isempty(theFigureTitle))
        title(theAxes, theFigureTitle);
    end
                
end

