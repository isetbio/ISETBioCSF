function plotSummaryData
    load('summaryData');
    
    for classifierIndex = 1:numel(classifiersList)
        classifierType = classifiersList{classifierIndex};
        s = summaryData(classifierType);
        thresholdsEnergy = [];
        thresholdContrasts = [];
        summationAreasList = [];
        for opticsModelIndex = 1:numel(s)
            m = s{opticsModelIndex};
            summationSigmasList = keys(m.spatialSummationData);
            for summationSigmaIndex = 1:numel(summationSigmasList)
                opticsLabels{opticsModelIndex} = m.optics;
                summationData = m.spatialSummationData(summationSigmasList{summationSigmaIndex});
                examinedSpotAreas = summationData.spotArea;
                summationAreasList(summationSigmaIndex) = summationData.summationArea;
                thresholdsEnergy(summationSigmaIndex,opticsModelIndex,:) = summationData.thresholdsEnergy;
                thresholdContrasts(summationSigmaIndex,opticsModelIndex,:) = summationData.thresholdContrasts;
            end
        end
        hFig = addPlot(classifierIndex, classifierType, summationAreasList, opticsLabels, examinedSpotAreas, thresholdsEnergy);
    end % classifierIndex
    
    NicePlot.exportFigToPDF('summaryData.pdf', hFig, 300);
end

function hFig = addPlot(classifierIndex, classifierType, summationAreasList, opticsLabels, examinedSpotAreas, thresholdsEnergy)
    
    hFig = figure(1);
    if (classifierIndex == 1) 
        clf;
        set(hFig, 'Position', [1 1 2200 805]);
    end
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 7, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
       
    coneApertureDiamMicrons = 3.0;
    ConeApertureAreaArcMinSquared = pi * (coneApertureDiamMicrons/2 * (60/300))^2;
    thresholdEnergyRange = [0.005 10];
    
    for k = 1:numel(summationAreasList)
        if (numel(summationAreasList) == 1)
            subplot('Position', subplotPosVectors(1,classifierIndex).v);
        else
            subplot('Position', subplotPosVectors(2,k).v);
        end
        plot(examinedSpotAreas, squeeze(thresholdsEnergy(k,1,:)), ...
            'ro-',  'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.7 0.7]);
        hold on;
        
        plot(examinedSpotAreas, squeeze(thresholdsEnergy(k,2,:)), ...
            'bo-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 1]);
        
        plot(ConeApertureAreaArcMinSquared*[1 1], thresholdEnergyRange, 'kv:', 'LineWidth', 1.5);
       
        plot(summationAreasList(k)*[1 1], thresholdEnergyRange, 'kv-', 'LineWidth', 1.5);
        hold off;
        
        set(gca, 'XLim', [1e-1 1e2], 'YLim', thresholdEnergyRange, ...
            'XTick', 10.^[-2:1:4], 'XScale', 'log', 'YScale', 'log', ...
            'FontSize', 16);
        axis 'square'
        grid on;
        
        if (k == 1) && (numel(summationAreasList)>1)
            legend({opticsLabels{1}, opticsLabels{2}, 'cone area', 'pooling area'}, 'Location', 'SouthEast', 'Orientation', 'Vertical');
            xlabel('spot area (arc min^2)', 'FontSize', 18, 'FontWeight', 'bold');
            ylabel('threshold energy', 'FontSize', 18, 'FontWeight', 'bold');
        else
            set(gca, 'YTickLabel', {});
        end
        if (~isnan(summationAreasList(k)))
            title(sprintf('classifier: %s\nsum. area: %2.2f arc min^2', classifierType, summationAreasList(k)));
        else
            title(sprintf('classifier:%s\nsum. area: none',  classifierType));
        end
    end

end

