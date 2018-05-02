function plotMacularPigmentTransmittance

    lens = Lens();
    lambda = lens.wave;
    
    macular = Macular('wave', lambda);
    macularTransmittance = macular.transmittance;

    
    xLims = [lambda(1) lambda(end)];
    yLims = [0 1];
    xTicks = 400:50:850;
    yTicks = 0:0.1:1;
    
    localDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    
    export = struct(...
        'format','PDF', ...
        'name', fullfile(localDir, 'MacularPigmentTransmittance.pdf') ...
        );
    renderXYplot(lambda, macularTransmittance, xLims, yLims, xTicks, yTicks, 'wavelength (nm)', 'transmittance', export);
end

function renderXYplot(lambda,macularTransmittance, xLims, yLims, xTicks, yTicks, xLabel, yLabel, export)
    % Render figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'PIGMENT_TRANSMITTANCE');
    
    outline.x = [xLims(1) xLims(1) xLims(2) xLims(2) xLims(1)];
    outline.y = [yLims(1) yLims(2) yLims(2) yLims(1) yLims(1)];
    
    colors = brewermap(3,'Set2');
    markerEdgeColor = [0 0 0];
    
    % Lens transmittance on the top plot
    subplot('Position', [0.12 0.12 0.86 0.88]);
%     plot(lambda, lensTransmittance, '-o', 'Color', squeeze(colors(1,:)), 'LineWidth', 1.5, ...
%         'MarkerFaceColor', squeeze(colors(1,:)), 'MarkerEdgeColor', markerEdgeColor, 'MarkerSize', 9);
%     hold on;
    plot(lambda, macularTransmittance, '-o', 'Color', squeeze(colors(1,:)), 'LineWidth', 1.5, ...
        'MarkerFaceColor', squeeze(colors(1,:)), 'MarkerEdgeColor', markerEdgeColor, 'MarkerSize', 9);
%     plot(lambda, macularTransmittance.*lensTransmittance, '-o', 'Color', squeeze(colors(3,:)), 'LineWidth', 1.5, ...
%         'MarkerFaceColor', squeeze(colors(3,:)), 'MarkerEdgeColor', markerEdgeColor, 'MarkerSize', 9);
    dy = 0.02*(yLims(2)-yLims(1));
    dx = 0.02*(xLims(2)-xLims(1));
    %plot(outline.x, outline.y, 'k-', 'LineWidth', 0.5);
    set(gca, 'XLim', [xLims(1)-dx xLims(2)+dx], 'YLim', [yLims(1)-dy yLims(2)+dy],...
        'XTick', xTicks, 'YTick', yTicks, 'LineWidth', 0.75);
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel(xLabel, 'FontWeight', 'bold');
    ylabel(yLabel, 'FontWeight', 'bold');
    grid on; box on;
    axis 'square';
    hL = legend({'lens', 'macular pigment', 'total'}, 'Location', 'SouthEast');
    
    t = text(410, 0.95, ' A ');
    formatFigureForPaper(hFig, 'figureType', 'PIGMENT_TRANSMITTANCE', 'theAxes', gca, 'theLegend', hL, 'theText', t);
    
    if strcmp(export.format, 'PDF')
        NicePlot.exportFigToPDF(export.name, hFig, 300);
    end
end