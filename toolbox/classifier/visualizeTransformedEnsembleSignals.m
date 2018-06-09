function hFigs = visualizeTransformedEnsembleSignals(V1filterEnsemble, timeAxis, noStimResponseInstances, stimResponseInstances, signalSource, stimContrast, spatialFilterName)

    if (strcmp(signalSource, 'isomerizations'))
        plotType = 'density';
    else
        plotType = 'line';
    end
            
    responseRange = [...
        min([min(noStimResponseInstances(:)) min(stimResponseInstances(:))]) ...
        max([max(noStimResponseInstances(:)) max(stimResponseInstances(:))]) ...
        ];
    
    responseQuantizationLevelsNum = 100;
    if (strcmp(signalSource, 'isomerizations'))
        responseRange(1) = 0;
        responseQuantizationLevelsNum = round(responseRange(2));
    end
    
    
    hFigs = visualizeEnsembleResponses(V1filterEnsemble, ...
        timeAxis, noStimResponseInstances, stimResponseInstances,  ...
        responseRange, responseQuantizationLevelsNum, plotType, signalSource, ...
        sprintf('%s (%s)', spatialFilterName, signalSource), 5003+sum(signalSource-'a'));

end

function hFigs = visualizeEnsembleResponses(V1filterEnsemble, timeAxis, noStimResponseInstances, stimResponseInstances, responseRange, responseLevelsNum, plotType, signalSource, figName, figNo)

    responseLevels = linspace(responseRange(1), responseRange(2), responseLevelsNum);
    responseInstancesNum = size(noStimResponseInstances,1);
    unitsNums = size(noStimResponseInstances,2);
    
    timeAxis = timeAxis*1000;
    hFigs = [];
   
    if (iscell(V1filterEnsemble))
        for unitIndex = 1:unitsNums
           bandwidthIndex = V1filterEnsemble{unitIndex}.bandwidthIndex;
           orientationIndex = V1filterEnsemble{unitIndex}.orientationIndex;
           ft2DIndex = V1filterEnsemble{unitIndex}.ft2Dindex;

           if (numel(hFigs) < ft2DIndex)
                rows = V1filterEnsemble{unitIndex}.rowsNum;
                cols = V1filterEnsemble{unitIndex}.colsNum;
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

                hFigs(ft2DIndex) = figure(figNo+ft2DIndex); clf;
                set(hFigs(ft2DIndex), 'Position', [10+ft2DIndex*100 10+ft2DIndex*50 1490 1340], 'Color', [1 1 1], 'Name', figureName);
                set(gcf,'renderer','opengl');
           end

           row = V1filterEnsemble{unitIndex}.rowColPosition(1) + halfRows+1;
           col = V1filterEnsemble{unitIndex}.rowColPosition(2) + halfCols+1;
           ax1 = subplot('Position', subplotPosVectors(rows+1-row,col).v);
           ax2 = [];
           
           noStimResponses = squeeze(noStimResponseInstances(:,unitIndex,:));
           stimResponses = squeeze(stimResponseInstances(:,unitIndex,:));
           renderNullTestComboResponse(ax1, ax2, ...
               signalSource, ...
               noStimResponses, stimResponses, [],[], ...
               responseLevels, timeAxis, plotType, row, col, rows);
       end % unitIndex
    else
        hFigs = figure(3333); clf;
        formatFigureForPaper(hFigs, ...
            'figureType','RESPONSE_INSTANCE_SINGLE_CONDIITION');
        
        ax1 = subplot('Position', [0.065 0.11 0.67 0.82]);
        ax2 = subplot('Position', [0.78  0.11 0.20 0.82]);
        
        timeAxisLimits = renderNullTestComboResponse(ax1, ax2, ...
            signalSource, ...
            noStimResponseInstances, stimResponseInstances, [], [], ...
            responseLevels, timeAxis, plotType, 1, 1, 1);
        
        formatFigureForPaper(hFigs, ...
            'figureType','RESPONSE_INSTANCE_SINGLE_CONDIITION', ...
            'theAxes', ax1);
        
        hL = legend(ax2, 'null stimulus', 'test stimulus');
        formatFigureForPaper(hFigs, ...
            'figureType','RESPONSE_INSTANCE_SINGLE_CONDIITION', ...
            'theAxes', ax2, ...
            'theLegend', hL);
        
        if (strcmp(signalSource, 'isomerizations'))
            ytickformat(ax1, '%.0f');
        else
            ytickformat(ax1, '%.1f');
        end
        
        xtickformat(ax2, '%0.2f');
        grid(ax2, 'on');
        drawnow
    end 
end

