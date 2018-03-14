function formatFigureForPaper(hFig, varargin)
    p = inputParser;
    p.addParameter('theAxes', []);
    p.addParameter('theLegend', []);
    p.addParameter('figureType', 'CSF', @ischar);
    p.parse(varargin{:});
    
    figureType = p.Results.figureType;
    theAxes =  p.Results.theAxes;
    theLegend = p.Results.theLegend;
    
    switch (figureType)
        case 'CSF'
            if (isempty(theAxes)) && (isempty(theLegend))
                set(hFig, 'Color', [1 1 1], 'Position', [10 10 500 800]);
                subplot('Position', [0.16 0.07 0.83 0.92]);
            end
    
            if (~isempty(theAxes))
                csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
                csLims = [2 12000];
                sfTicks = [1 2 5 10 20 50 100];
                sfLims  = [1.2 100];
                dx1 = 0.2; dx2 = 20; dy1 = 0.4; dy2 = 1000;

                set(theAxes, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, ...
                    'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', [csLims(1)-dy1 csLims(2)+dy2]);
                box(theAxes, 'on'); grid(theAxes, 'on');
                set(theAxes, 'FontSize', 18, 'TickLength',[0.02, 0.02]);
                xlabel(theAxes,'spatial frequency (c/deg)', 'FontWeight', 'bold', 'FontSize', 20);
                ylabel(theAxes,'contrast sensitivity', 'FontWeight', 'bold', 'FontSize', 20);
            end
            
            if (~isempty(theLegend))
                set(theLegend, 'FontName', 'Monaco', 'FontSize', 14);
            end
    end % switch
    
end

