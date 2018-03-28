function hFigs = visualizeTransformedEnsembleSignals(V1filterEnsemble, timeAxis, noStimResponseInstances, stimResponseInstances, signalSource, stimContrast, spatialFilterName)

    responseQuantizationLevelsNum = 100;
    
    if (strcmp(signalSource, 'isomerizations'))
        plotType = 'density';
    else
        plotType = 'line';
    end
            
    responseRange = [...
        min([min(noStimResponseInstances(:)) min(stimResponseInstances(:))]) ...
        max([max(noStimResponseInstances(:)) max(stimResponseInstances(:))]) ...
        ];
    
    hFigs = visualizeEnsembleResponses(V1filterEnsemble, ...
        timeAxis, noStimResponseInstances, stimResponseInstances,  ...
        responseRange, responseQuantizationLevelsNum, plotType, ...
        sprintf('%s (%s)', spatialFilterName, signalSource), 5003+sum(signalSource-'a'));

end

function hFigs = visualizeEnsembleResponses(V1filterEnsemble, timeAxis, noStimResponseInstances, stimResponseInstances, responseRange, responseLevelsNum, plotType, yAxisLabel, figNo)


    responseLevels = linspace(responseRange(1), responseRange(2), responseLevelsNum);
    
    
    responseInstancesNum = size(noStimResponseInstances,1);
    unitsNums = size(noStimResponseInstances,2);
    
    rows = abs(V1filterEnsemble{1}.rowColPosition(1))*2+1;
    cols = abs(V1filterEnsemble{1}.rowColPosition(2))*2+1;
    halfRows = (rows-1)/2;
    halfCols = (cols-1)/2;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', rows, ...
       'colsNum', cols, ...
       'heightMargin',   0.02, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.04, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.01);
        
   timeAxis = timeAxis*1000;
   hFigs = [];
   
   for unitIndex = 1:unitsNums
       bandwidthIndex = V1filterEnsemble{unitIndex}.bandwidthIndex;
       orientationIndex = V1filterEnsemble{unitIndex}.orientationIndex;
       ft2DIndex = V1filterEnsemble{unitIndex}.ft2Dindex;
       
       if (numel(hFigs) < ft2Dindex)
            hFigs(ft2Dindex) = figure(figNo+ft2Dindex); clf;
            set(hFigs(ft2Dindex), 'Position', [10+ft2Dindex*100 10+ft2Dindex*50 1490 1340], 'Color', [1 1 1]);
            set(gcf,'renderer','opengl');
       end
       
       row = V1filterEnsemble{unitIndex}.rowColPosition(1) + halfRows+1;
       col = V1filterEnsemble{unitIndex}.rowColPosition(2) + halfCols+1;
       subplot('Position', subplotPosVectors(row,col).v);
       
       % Concatenate responses: null - test 
       unitResponses= squeeze(noStimResponseInstances(:,unitIndex,:));
       unitResponses = cat(2, unitResponses, squeeze(stimResponseInstances(:,unitIndex,:)));
       meanResponse = mean(unitResponses,1);
       
       timeAxisComboSecondResponseStart = timeAxis(end)+timeAxis(2)-timeAxis(1);
       timeAxisCombo = cat(2, timeAxis, timeAxis(end)+timeAxis);
       timeAxisCombo = timeAxisCombo+(timeAxisCombo(2)-timeAxisCombo(1));
       
       timeAxisLimits = [timeAxisCombo(1) timeAxisCombo(end)];
       if (timeAxisCombo(1) == timeAxisCombo(end))
            timeAxisLimits = timeAxisCombo(1) + [-0.5 0.5];
       end
    
       if (strcmp(plotType, 'density'))
            responseDistribution = compute2Dhistogram(unitResponses, responseLevels);
            imagesc(timeAxisCombo, responseLevels, responseDistribution);
            hold on;
       elseif (strcmp(plotType, 'line'))
           hold on;
           for instanceIndex = 1:size(responseInstancesNum,1)
               y = squeeze(unitResponses(instanceIndex, :));
               xflip = [timeAxisCombo(1 : end - 1) fliplr(timeAxisCombo)];
               yflip = [y(1 : end - 1) fliplr(y)];
               patch(xflip, yflip, [0.1 0.1 0.1], 'EdgeAlpha', 0.1, 'FaceColor', 'none');
           end 
       end
       
       plot(timeAxisCombo, meanResponse, 'w-', 'LineWidth', 6.0);
       plot(timeAxisCombo, meanResponse, 'r-', 'LineWidth', 2.0);
       plot(timeAxisComboSecondResponseStart*[1 1], [responseLevels(1) responseLevels(end)], 'k-');
       hold off;
       grid on; box on;
       set(gca, 'XTick', [0:25:200]);
       if (row == rows) && (col == 1)
           xlabel('time (msec)');
           ylabel('response');
       else
           set(gca, 'XTickLabel', {}, 'YTickLabel', {});
       end
       
       set(gca, 'CLim', [0 1]);
       set(gca, 'YLim', [responseLevels(1) responseLevels(end)], 'XLim', timeAxisLimits, 'FontSize', 14);
       axis 'xy';
       
       colormap(1-gray(1024));
       drawnow; 
   
   end % unitIndex
       
   
end

function responseDistribution = compute2Dhistogram(responses, responseLevels)

    responseDistribution = zeros(numel(responseLevels)-1, size(responses,2));
    for tBin = 1:size(responses,2)
        N = histcounts(squeeze(responses(:,tBin)),responseLevels);
        responseDistribution(:,tBin) = N';
    end
    responseDistribution = responseDistribution/max(responseDistribution(:));
end
