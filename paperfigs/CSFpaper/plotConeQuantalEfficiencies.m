function plotConeQuantalEfficiencies

    c = coneMosaic;
    lambda = c.wave;
    quantalEfficiencies(:,1) = c.qe(:,1) ./ c.macular.transmittance;
    quantalEfficiencies(:,2) = c.qe(:,2) ./ c.macular.transmittance;
    quantalEfficiencies(:,3) = c.qe(:,3) ./ c.macular.transmittance;
    

    xLims = [lambda(1) lambda(end)];
    yLims = [0 0.5];
    xTicks = 400:50:850;
    yTicks = 0:0.1:1;
    
    localDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    
    export = struct(...
        'format','PDF', ...
        'name', fullfile(localDir, 'ConeQuantalEfficiencies.pdf') ...
        );
    renderXYplot(lambda, quantalEfficiencies', xLims, yLims, xTicks, yTicks, 'wavelength (nm)', 'transmittance', export);
end

function renderXYplot(lambda,quantalEfficiencies, xLims, yLims, xTicks, yTicks, xLabel, yLabel, export)
    % Render figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'PIGMENT_TRANSMITTANCE');
    

    dy = 0.02*(yLims(2)-yLims(1));
    dx = 0.02*(xLims(2)-xLims(1));
    colors = brewermap(3,'Set2');
    markerEdgeColor = [0 0 0];
    
    % Lens transmittance on the top plot
    subplot('Position', [0.12 0.12 0.86 0.88]);
    hold on
    for coneIndex = 1:3
        plot(lambda, squeeze(quantalEfficiencies(coneIndex,:)), '-o', 'Color', squeeze(colors(coneIndex,:)), 'LineWidth', 1.5, ...
        'MarkerFaceColor', squeeze(colors(coneIndex,:)), 'MarkerEdgeColor', markerEdgeColor, 'MarkerSize', 9);
    end
    set(gca, 'XLim', [xLims(1)-dx xLims(2)+dx], 'YLim', [yLims(1)-dy yLims(2)+dy],...
        'XTick', xTicks, 'YTick', yTicks, 'LineWidth', 0.75);
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel(xLabel, 'FontWeight', 'bold');
    ylabel(yLabel, 'FontWeight', 'bold');
    grid on; box on;
    axis 'square';
    hL = legend({'L-cones', 'M-cones', 'S-cones'});
    
    t = text(410, 0.475, ' B ');
    formatFigureForPaper(hFig, 'figureType', 'PIGMENT_TRANSMITTANCE', 'theAxes', gca, 'theLegend', hL, 'theText', t);
    set(hL, 'Location', 'NorthEast');
    if strcmp(export.format, 'PDF')
        NicePlot.exportFigToPDF(export.name, hFig, 300);
    end
end

