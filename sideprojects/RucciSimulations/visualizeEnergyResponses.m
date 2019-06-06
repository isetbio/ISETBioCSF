function visualizeEnergyResponses(stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, resourcesDir, figNo)

    fname = fullfile(resourcesDir, sprintf('energyResponse_%s_instance_%1.0f_nTrials_%d.mat', stimDescriptor, analyzedNoiseInstance, nTrials));
    load(fname, 'energyResponse','emPathsDegs', 'timeAxis');
   
    nContrasts = numel(contrastLevels);
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 2500 350]);
    
    for theContrastLevel = 1:nContrasts
        standardOrientationMechanismResponses = squeeze(energyResponse.output(theContrastLevel,:,:));
        orthogonalOrientationMechanismResponses = squeeze(energyResponse.orthoOutput(theContrastLevel,:,:));
        subplot(1,nContrasts,theContrastLevel)
        hold on;
        for trialIndex = 1:nTrials
            plot(timeAxis*1000, standardOrientationMechanismResponses(trialIndex,:), 'r-');
            plot(timeAxis*1000, orthogonalOrientationMechanismResponses(trialIndex,:), 'b-');
        end
        xlabel('time (msec)');
        ylabel('energy response');
        legend({'standard ori', 'orthogonal ori'});
        title(sprintf('c = %2.3f (%s)', contrastLevels(theContrastLevel),stimDescriptor));
    end
    
end
