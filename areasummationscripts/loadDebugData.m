function loadDebugData

    employedOptics = 'AOoptics';
    datafile = fullfile('3levels',sprintf('SummationDataExtendedRange_%s.mat',employedOptics));
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k)
        data3levels = spatialSummationData(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin))
    end

    datafile = fullfile('1onlylevel',sprintf('SummationDataExtendedRange_%s.mat',employedOptics));
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k)
        data1level = spatialSummationData(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin))
    end
    
    minspotsizeRatios = data3levels.spotArea(1)/data1level.spotArea
    thresholdRatios = data3levels.thresholdContrasts(3)/data1level.thresholdContrasts(1)
    thresholdEnergyRatios = data3levels.thresholdsEnergy(3)/data1level.thresholdsEnergy(1)
    
    thresholdRatiosDividedByMinSpotSizeRatios = thresholdRatios / minspotsizeRatios
end

