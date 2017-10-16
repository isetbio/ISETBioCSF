function hFig = visualizeTransformedSignals(timeAxis, noStimResponseInstances, stimResponseInstances, signalSource, stimContrast, spatialFilterName)

    responseQuantizationLevelsNum = 100;
    p = 0;
    
    noStimNoiseFreeResponse = squeeze(mean(noStimResponseInstances,1));
    stimNoiseFreeResponse = squeeze(mean(stimResponseInstances,1));
    
    if (strcmp(signalSource, 'isomerizations'))
        plotType = 'density';
    else
        plotType = 'line';
    end
            
    responseRange = [...
        min([min(noStimResponseInstances(:)) min(stimResponseInstances(:))]) ...
        max([max(noStimResponseInstances(:)) max(stimResponseInstances(:))]) ...
        ];
    
    hFig = visualizeResponsesInstancesAndNoiseFreeResponsesAsDensityPlots(...
        timeAxis, noStimResponseInstances, stimResponseInstances, ...
        noStimNoiseFreeResponse, stimNoiseFreeResponse,  ...
        p, responseRange, responseQuantizationLevelsNum, plotType, sprintf('%s (%s)', spatialFilterName, signalSource), 5003+sum(signalSource-'a'));
end
