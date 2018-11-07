function hFig = visualizeBestRespondingLMSResponseInstancesAndNoiseFreeResponse(...
    timeAxis, noStimResponseInstances, stimResponseInstances, ...
    noStimNoiseFreeResponse, stimNoiseFreeResponse, responseRange, ...
    responseLevelsNum, coneTypeIndicesToVisualize, plotType, signalSource, yAxisLabel, figNo)

    % Select what we are visualizing
    noStimResponseInstances = noStimResponseInstances(coneTypeIndicesToVisualize,:,:);
    stimResponseInstances = stimResponseInstances(coneTypeIndicesToVisualize,:,:);
    noStimNoiseFreeResponse = noStimNoiseFreeResponse(coneTypeIndicesToVisualize,:);
    stimNoiseFreeResponse = stimNoiseFreeResponse(coneTypeIndicesToVisualize,:);
    
    % Arrange subplot regions
    responseNums = size(noStimResponseInstances,1);
    if (responseNums > 1)        
        for k = 1:responseNums
            subplotPosVectors(responseNums-k+1,1).v = [0.065 0.05+(k-1)*0.33 0.67 0.25];
            subplotPosVectors(responseNums-k+1,2).v = [0.78 0.05+(k-1)*0.33 0.19 0.25];
        end
    elseif (responseNums == 1)
        subplotPosVectors(1,1).v = [0.065 0.11 0.67 0.83];
        subplotPosVectors(1,2).v = [0.78 0.11 0.19 0.83];
    else
        error('responseNums in %s is: %d\n', mfilename(), responseNums)
    end


    
    plotType = 'density';
    responseLevels = linspace(responseRange(1), responseRange(2), responseLevelsNum);
    
    hFig = figure(figNo); clf;
    if (responseNums > 1)
        formatFigureForPaper(hFig, ...
            'figureType','RESPONSE_INSTANCE_THREE_CONDIITIONS');
    else
        formatFigureForPaper(hFig, ...
            'figureType','RESPONSE_INSTANCE_SINGLE_CONDIITION');
    end
    
    for responseIndex = 1:responseNums
        ax1 = subplot('Position', subplotPosVectors(responseIndex,1).v);
        ax2 = subplot('Position', subplotPosVectors(responseIndex,2).v);
        
        if (responseNums == 1)
            timeAxisLimits = renderNullTestComboResponse(ax1, ax2, signalSource, ...
                squeeze(noStimResponseInstances(responseIndex,:,:)), squeeze(stimResponseInstances(responseIndex,:,:)), ...
                squeeze(noStimNoiseFreeResponse(responseIndex,:)), squeeze(stimNoiseFreeResponse(responseIndex,:)), ...
                [min(responseLevels) max(responseLevels)], timeAxis, plotType, 1, 1, 1);
        else
            timeAxisLimits = renderNullTestComboResponse(ax1, ax2, signalSource, ...
                squeeze(noStimResponseInstances(responseIndex,:,:)), squeeze(stimResponseInstances(responseIndex,:,:)), ...
                squeeze(noStimNoiseFreeResponse(responseIndex,:)), squeeze(stimNoiseFreeResponse(responseIndex,:)), ...
                [min(responseLevels) max(responseLevels)], timeAxis, plotType, 1, 1, 1);
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
        
        if (responseNums > 1)
            formatFigureForPaper(hFig, ...
                'figureType','RESPONSE_INSTANCE_THREE_CONDIITIONS', ...
                'theAxes', ax1, ...
                'theText', t);

            hL = legend(ax2, 'null stimulus', 'test stimulus');
            formatFigureForPaper(hFig, ...
                'figureType','RESPONSE_INSTANCE_THREE_CONDIITIONS', ...
                'theAxes', ax2, ...
                'theLegend', hL);
        else
            formatFigureForPaper(hFig, ...
                'figureType','RESPONSE_INSTANCE_SINGLE_CONDIITION', ...
                'theAxes', ax1, ...
                'theText', t);

            hL = legend(ax2, 'null stimulus', 'test stimulus');
            formatFigureForPaper(hFig, ...
                'figureType','RESPONSE_INSTANCE_SINGLE_CONDIITION', ...
                'theAxes', ax2, ...
                'theLegend', hL);
        end
        
        
        if (strcmp(signalSource, 'isomerizations'))
            ytickformat(ax1, '%.0f');
        else
            ytickformat(ax1, '%.1f');
        end
        xtickformat(ax2, '%0.2f');
        grid(ax2, 'on');
        
    end % responseIndex
end

