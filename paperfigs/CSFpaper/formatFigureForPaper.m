function varargout = formatFigureForPaper(hFig, varargin)
    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.addParameter('plotRatiosOfOtherConditionsToFirst', false, @islogical);
    p.addParameter('theAxes', []);
    p.addParameter('theRatioAxes', []);
    p.addParameter('theRatioLims', []);
    p.addParameter('theRatioTicks', []);
    p.addParameter('theLegend', []);
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
    theRatioAxes = p.Results.theRatioAxes;
    theRatioLims = p.Results.theRatioLims;
    theRatioTicks = p.Results.theRatioTicks;
    theLegend = p.Results.theLegend;
    theText = p.Results.theText;
    theTextFontSize = p.Results.theTextFontSize;
    theFigureTitle = p.Results.theFigureTitle;
    
    switch (figureType)
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
            csLims = [2 12000];
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
                    if (~isempty(theRatioTicks))
                        set(theRatioAxes, 'YTick', theRatioTicks);
                    end
                    if (~isempty(theRatioLims))
                        set(theRatioAxes, 'YLim', theRatioLims);
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

