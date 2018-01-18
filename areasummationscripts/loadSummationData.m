function loadSummationData()

    classifiersList = {'mlpt', 'svm', 'svmGaussianRF'};
    
    summaryData = containers.Map();
    
    for classifierIndex = 1:numel(classifiersList)
        
        classifierType = classifiersList{classifierIndex};
        
        employedOptics1 = 'AOoptics';
        employedOptics1legend = 'AOoptics';
        datafile = sprintf('SummationData_%s_%s.mat', classifierType, employedOptics1);
        load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
        spatialSummationData1 = spatialSummationData;
        spatialPoolingSigmaArcMinList1 = spatialPoolingSigmaArcMinList;
        
        employedOptics2 = 'WvfHumanMeanOTFmagMeanOTFphase';
        employedOptics2legend = 'wvfHuman';
        datafile = sprintf('SummationData_%s_%s.mat',  classifierType, employedOptics2);
        load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
        spatialSummationData2 = spatialSummationData;
        spatialPoolingSigmaArcMinList2 = spatialPoolingSigmaArcMinList;
        
        summaryData(classifierType) = {
            struct( ...
                'optics', employedOptics1legend, ...
                'spatialSummationData', spatialSummationData1 ...  
            ), ...
            struct( ...
                'optics', employedOptics2legend, ...
                'spatialSummationData', spatialSummationData2 ...
            )
        };
    
        spatialSummationData3 = [];
        employedOptics3 = '';
        employedOptics3legend = '';

        % Plot the data
        hFig = plotTheData(classifierType, spatialPoolingSigmaArcMinList, spatialSummationData1, spatialSummationData2, spatialSummationData3, ...
            employedOptics1legend, employedOptics2legend, employedOptics3legend);
    
        NicePlot.exportFigToPDF(sprintf('%s_summary.pdf', classifierType), hFig, 300);
    end
    
    save('summaryData.mat', 'summaryData', 'classifiersList');
end

function hFig = plotTheData(classifierType, spatialPoolingSigmaArcMinList, spatialSummationData1, spatialSummationData2, spatialSummationData3, employedOptics1, employedOptics2, employedOptics3)
    % Figure setup
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', numel(spatialPoolingSigmaArcMinList), ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.006, ...
           'bottomMargin',   0.15, ...
           'topMargin',      0.05);
       
    if strcmp(classifierType, 'svmGaussianRF')
       figSize = [2400 500];
    else
       figSize = [500 500];
    end
        
    hFig = figure(1); clf;
    set(hFig, 'Position', [1 500 figSize(1) figSize(2)], 'Color', [1 1 1]);
    
    coneApertureDiamMicrons = 3.0;
    ConeApertureAreaArcMinSquared = pi * (coneApertureDiamMicrons/2 * (60/300))^2;
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        % Get spatial pooling sigma
        
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k);
        if ~isnan(spatialPoolingSigmaArcMin)
            % Get threshold data struct associated with this spatial pooling sigma
            plotData1 = spatialSummationData1(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin));
            plotData2 = spatialSummationData2(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin));
            if (~isempty(spatialSummationData3))
                plotData3 = spatialSummationData3(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin));
            end
        else
            plotData1 = spatialSummationData1('none');
            plotData2 = spatialSummationData2('none');
            if (~isempty(spatialSummationData3))
                plotData3 = spatialSummationData3('none');
            end
        end
        
        % Render the plot
        thresholdRange = [0.5*1e-4 1e-1];
        thresholdEnergyRange = [0.005 10];
        subplot('Position', subplotPosVectors(1,k).v);
        plot(plotData1.spotArea, plotData1.thresholdsEnergy, 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [0.9 0.6 0.7], 'LineWidth', 1.5);
        hold on;
        plot(plotData2.spotArea, plotData2.thresholdsEnergy, 'bo-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.6 0.9], 'LineWidth', 1.5);
        if (~isempty(spatialSummationData3))
            plot(plotData3.spotArea, plotData3.thresholdsEnergy, 'ko-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.6 0.9], 'LineWidth', 1.5);
        end
        plot(plotData1.summationArea*[1 1], thresholdEnergyRange, 'kv-', 'LineWidth', 1.5);
        plot(ConeApertureAreaArcMinSquared*[1 1], thresholdEnergyRange, 'mv--', 'LineWidth', 1.5);
        hold off
        set(gca, 'XLim', [1e-1 1e2], 'YLim', thresholdEnergyRange, 'XTick', 10.^[-2:1:4], 'XScale', 'log', 'YScale', 'log', 'FontSize', 16);
        axis 'square'
        grid on;
        xlabel('spot area (arc min^2)', 'FontSize', 18, 'FontWeight', 'bold');
        if (k==1)
            ylabel('threshold energy', 'FontSize', 18, 'FontWeight', 'bold');
        end
        if (k == 1)
            if (~isempty(spatialSummationData3))
                legend({employedOptics1, employedOptics2, employedOptics3, 'summation area', 'cone area'}, 'Location', 'SouthEast', 'Orientation', 'Vertical');
            else
                legend({employedOptics1, employedOptics2, 'summation area', 'cone area'}, 'Location', 'SouthEast', 'Orientation', 'Vertical');
            end
        end
        
        if (~isnan(plotData1.summationArea))
            title(sprintf('sum. area: %2.2f arc min^2\nclassifier:%s', plotData1.summationArea, classifierType));
        else
            title(sprintf('classifier:%s', classifierType));
        end
    end
    
    
end