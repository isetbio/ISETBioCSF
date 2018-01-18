function loadDebugData

    employedOptics = 'AOoptics';
    datafile = fullfile('5levels',sprintf('SummationDataExtendedRange_%s.mat',employedOptics));
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k)
        data5levels = spatialSummationData(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin))
    end

    datafile = fullfile('1levellast',sprintf('SummationDataExtendedRange_%s.mat',employedOptics));
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k);
        data1levellast = spatialSummationData(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin))
    end
    
    datafile = fullfile('1levelthird',sprintf('SummationDataExtendedRange_%s.mat',employedOptics));
    load(datafile, 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
    
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k);
        data1level = spatialSummationData(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin))
    end
    
%     minspotsizeRatios = data5levels.spotArea(1)/data1levellast.spotArea
%     thresholdRatios = data3levels.thresholdContrasts(3)/data1levellast.thresholdContrasts(1)
%     thresholdEnergyRatios = data3levels.thresholdsEnergy(3)/data1levellast.thresholdsEnergy(1)
%     
%     thresholdRatiosDividedByMinSpotSizeRatios = thresholdRatios / minspotsizeRatios
end

