function hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(timeAxis, noStimResponseInstances, stimResponseInstances, noStimNoiseFreeResponse, stimNoiseFreeResponse, responseRange, responseLevelsNum, plotType, signalSource, yAxisLabel, figNo)


    if (numel(timeAxis) > 1)
        if (ndims(noStimResponseInstances) == 3)
            responseNums = size(noStimResponseInstances,1);
            for k = 1:responseNums
                subplotPosVectors(responseNums-k+1,1).v = [0.065 0.05+(k-1)*0.33 0.67 0.25];
                subplotPosVectors(responseNums-k+1,2).v = [0.78 0.05+(k-1)*0.33 0.20 0.25];
            end
        else
            responseNums = 1;
            subplotPosVectors(1,1).v = [0.065 0.11 0.67 0.85];
            subplotPosVectors(1,2).v = [0.78 0.11 0.20 0.85];
        end
    else
       if (ndims(noStimResponseInstances) == 2)
            responseNums = size(noStimResponseInstances,1);
            for k = 1:responseNums
                subplotPosVectors(responseNums-k+1,1).v = [0.065 0.05+(k-1)*0.33 0.67 0.25];
                subplotPosVectors(responseNums-k+1,2).v = [0.78 0.05+(k-1)*0.33 0.20 0.25];
            end
        else
            responseNums = 1;
            subplotPosVectors(1,1).v = [0.065 0.11 0.67 0.85];
            subplotPosVectors(1,2).v = [0.78 0.11 0.20 0.85];
        end 
    end
    
    plotType = 'density';
    responseLevels = linspace(responseRange(1), responseRange(2), responseLevelsNum);
    
    hFig = figure(figNo); clf;
    formatFigureForPaper(hFig, ...
            'figureType','RESPONSE_INSTANCE_THREE_CONDIITIONS');
         
    for responseIndex = 1:responseNums
        ax1 = subplot('Position', subplotPosVectors(responseIndex,1).v);
        ax2 = subplot('Position', subplotPosVectors(responseIndex,2).v);
        
        if (responseNums == 1)
            timeAxisLimits = renderNullTestComboResponse(ax1, ax2, ...
                noStimResponseInstances, stimResponseInstances, signalSource, ...
                noStimNoiseFreeResponse, stimNoiseFreeResponse, ...
                responseLevels, timeAxis, plotType, 1, 1, 1);
        else
            timeAxisLimits = renderNullTestComboResponse(ax1, ax2, signalSource, ...
                squeeze(noStimResponseInstances(responseIndex,:,:)), squeeze(stimResponseInstances(responseIndex,:,:)), ...
                squeeze(noStimNoiseFreeResponse(responseIndex,:)), squeeze(stimNoiseFreeResponse(responseIndex,:)), ...
                responseLevels, timeAxis, plotType, 1, 1, 1);
        end
        
        t = [];
        if (responseNums > 1)
            if (responseIndex == 1)
               %t = text(ax1, 8, responseLevels(end)*0.87, ' A ');
               %xlabel(ax1, '');
               %xlabel(ax2, '');
            elseif (responseIndex == 2)
               %t = text(ax1, 8, responseLevels(end)*0.87, ' B ');
               %xlabel(ax1, '');
               %xlabel(ax2, '');
            elseif (responseIndex == 3)
               %t = text(ax1, 8, responseLevels(end)*0.87, ' C ');
            end
            
        end  
        
        formatFigureForPaper(hFig, ...
            'figureType','RESPONSE_INSTANCE_THREE_CONDIITIONS', ...
            'theAxes', ax1, ...
            'theText', t);
        
        hL = legend(ax2, 'null stimulus', 'test stimulus');
        formatFigureForPaper(hFig, ...
            'figureType','RESPONSE_INSTANCE_THREE_CONDIITIONS', ...
            'theAxes', ax2, ...
            'theLegend', hL);
        
        if (strcmp(signalSource, 'isomerizations'))
            ytickformat(ax1, '%.0f');
        else
            ytickformat(ax1, '%.1f');
        end
        xtickformat(ax2, '%0.2f');
        grid(ax2, 'on');
        
    end % responseIndex
end

