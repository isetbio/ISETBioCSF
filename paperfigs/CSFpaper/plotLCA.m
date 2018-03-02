function plotLCA
    lambda = 400:850;
    lambdaFocus = 550;
    defocus = 633.46 * (1/(lambdaFocus - 214.1) - 1./(lambda-214.1));
    
    xLims = [lambda(1) lambda(end)];
    yLims = [-2 1];
    xTicks = 400:50:850;
    yTicks = -2:0.2:1;
    export = struct(...
        'format','PDF', ...
        'name', 'LCA.pdf' ...
        );
    renderXYplot(lambda,defocus, xLims, yLims, xTicks, yTicks, 'wavelength (nm)', 'defocus (D)', export);
end

function renderXYplot(x,y, xLims, yLims, xTicks, yTicks, xLabel, yLabel, export)
    % Render figure
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);
    
    outline.x = [xLims(1) xLims(1) xLims(2) xLims(2) xLims(1)];
    outline.y = [yLims(1) yLims(2) yLims(2) yLims(1) yLims(1)];
    plot(x, y, 'b-', 'LineWidth', 2);
    hold on;
    dy = 0.02*(yLims(2)-yLims(1));
    dx = 0.02*(xLims(2)-xLims(1));
    plot(outline.x, outline.y, 'k-', 'LineWidth', 0.5);
    set(gca, 'XLim', [xLims(1)-dx xLims(2)+dx], 'YLim', [yLims(1)-dy yLims(2)+dy],...
        'XTick', xTicks, 'YTick', yTicks, 'LineWidth', 0.75);
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel(xLabel, 'FontWeight', 'bold');
    ylabel(yLabel, 'FontWeight', 'bold');
    grid on; box on;
    axis 'square';
    if strcmp(export.format, 'PDF')
        NicePlot.exportFigToPDF(export.name, hFig, 300);
    end
end

