function loadSummationData()

    % Load wvf data
    employedOptics = 'None';
    datafile = sprintf('SummationDataExtendedRange_%s.mat',employedOptics);
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    spatialSummationDataWvfOptics = spatialSummationData;
    
    % Load no optics data
    employedOptics = 'None';
    datafile = sprintf('SummationDataExtendedRange_%s.mat',employedOptics);
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    spatialSummationDataNoOptics = spatialSummationData;
    
    % Plot the data
    plotTheData(spatialPoolingSigmaArcMinList, spatialSummationDataWvfOptics, spatialSummationDataNoOptics)
end

function plotTheData(spatialPoolingSigmaArcMinList, spatialSummationDataWvfOptics, spatialSummationDataNoOptics)
    % Figure setup
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', numel(spatialPoolingSigmaArcMinList), ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.1, ...
           'topMargin',      0.01);
       
    hFig = figure(1); clf;
    set(hFig, 'Position', [1 500 1670 600], 'Color', [1 1 1]);
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        % Get spatial pooling sigma
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k);
        
        % Get threshold data struct associated with this spatial pooling sigma
        plotDataWvfOptics = spatialSummationDataWvfOptics(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin));
        plotDataNoOptics = spatialSummationDataNoOptics(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin));
       
        % Extract selected data
        ISETbioSpotAreasArcMinSquared = plotDataWvfOptics.spotArea;
        ISETbioSummationAreaArcMinSquared = plotDataWvfOptics.summationArea;
        ISETbioThresholdContrastsWvfOptics = plotDataWvfOptics.thresholdContrasts;
        ISETbioThresholdContrastsNoOptics = plotDataNoOptics.thresholdContrasts;
        
        % Render the plot
        thresholdRange = [0.5*1e-4 1e-1];
        subplot('Position', subplotPosVectors(1,k).v);
        plot(ISETbioSpotAreasArcMinSquared, ISETbioThresholdContrastsWvfOptics, 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [0.9 0.6 0.7], 'LineWidth', 1.5);
        hold on;
        plot(ISETbioSpotAreasArcMinSquared, ISETbioThresholdContrastsNoOptics, 'bo-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.6 0.9], 'LineWidth', 1.5);
        plot(ISETbioSummationAreaArcMinSquared*[1 1], thresholdRange, 'kv-', 'LineWidth', 1.5);
        hold off
        set(gca, 'XLim', [1e-1 2*1e3], 'YLim', [0.5*1e-4 0.5*1e-1], 'XTick', 10.^[-2:1:4], 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
        grid on;
        xlabel('spot area (arc min^2)', 'FontSize', 18, 'FontWeight', 'bold');
        if (k==1)
            ylabel('threshold contrast', 'FontSize', 18, 'FontWeight', 'bold');
        else
            set(gca, 'YTickLabel', {});
        end
        legend('wvf optics', 'no optics', 'Location', 'NorthOutside', 'Orientation', 'Horizontal')
    end
    
    NicePlot.exportFigToPDF('summation.pdf', hFig, 300);
end